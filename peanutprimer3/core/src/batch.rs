//! Batch orchestration: parse input, find variants across all records, and
//! design each marker.
//!
//! Parsing/variant-scanning/candidate-generation/scoring are pure Rust and
//! safe to parallelize with `rayon`; the common-primer `design_primers()`
//! call is already internally mutex-guarded by the `primer3` crate itself
//! (see the doc comment on [`crate::common_primer`]), so no extra
//! synchronization is needed here -- rayon can call straight through
//! [`crate::common_primer::design_marker`] from multiple threads.

use rayon::prelude::*;
use std::sync::mpsc::Sender;

use crate::model::{MarkerResult, PresetParams};
use crate::parser;
use crate::variant::{self, Variant};

#[derive(Debug, Clone)]
pub enum ProgressEvent {
    Started { total_markers: usize },
    MarkerDone { index: usize, sequence_id: String },
    Complete,
}

struct Job {
    sequence_id: String,
    variant: Option<Variant>,
}

/// Runs primer design for every marker found in `input`. A record with no
/// detected variant yields one error-only [`MarkerResult`]; a record with
/// more than one variant yields one marker per variant (sharing the
/// sequence id), since each needs independent primer design.
pub fn run_batch(input: &str, params: &PresetParams, progress: Option<Sender<ProgressEvent>>) -> Vec<MarkerResult> {
    let records = parser::parse_records(input);

    let jobs: Vec<Job> = records
        .into_iter()
        .flat_map(|record| {
            let variants = variant::find_variants(&record.sequence);
            if variants.is_empty() {
                vec![Job { sequence_id: record.id, variant: None }]
            } else {
                variants.into_iter().map(|v| Job { sequence_id: record.id.clone(), variant: Some(v) }).collect()
            }
        })
        .collect();

    if let Some(tx) = &progress {
        let _ = tx.send(ProgressEvent::Started { total_markers: jobs.len() });
    }

    let results: Vec<MarkerResult> = jobs
        .into_par_iter()
        .enumerate()
        .map(|(index, job)| {
            let result = match &job.variant {
                Some(v) => crate::common_primer::design_marker(&job.sequence_id, v, params),
                None => MarkerResult {
                    sequence_id: job.sequence_id.clone(),
                    variant_description: "no variant found".to_string(),
                    flank_warning: None,
                    primary: None,
                    alternative: None,
                    errors: vec![
                        "No SNP/indel variant found in this sequence (expected an embedded IUPAC \
                         ambiguity code or bracket notation like [A/G])"
                            .to_string(),
                    ],
                },
            };
            if let Some(tx) = &progress {
                let _ = tx.send(ProgressEvent::MarkerDone { index, sequence_id: job.sequence_id.clone() });
            }
            result
        })
        .collect();

    if let Some(tx) = &progress {
        let _ = tx.send(ProgressEvent::Complete);
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    fn repeat_to_len(unit: &str, len: usize) -> String {
        unit.chars().cycle().take(len).collect()
    }

    #[test]
    fn run_batch_designs_multiple_markers() {
        let left = repeat_to_len("CTGACGATCGTAGGCATCGA", 150);
        let right = repeat_to_len("TAGCCTGAACGGTCATGCAT", 150);
        let input = format!(">marker_1\n{left}R{right}\n>marker_2\n{left}Y{right}\n");

        let results = run_batch(&input, &PresetParams::cotton_kasp_hulse_kemp_2015(), None);
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].sequence_id, "marker_1");
        assert_eq!(results[1].sequence_id, "marker_2");
        for r in &results {
            assert!(r.errors.is_empty(), "unexpected errors for {}: {:?}", r.sequence_id, r.errors);
            assert!(r.primary.is_some(), "expected a primary triplet for {}", r.sequence_id);
        }
    }

    #[test]
    fn run_batch_reports_error_for_record_without_a_variant() {
        let input = ">no_variant\nACGTACGTACGTACGTACGT\n";
        let results = run_batch(input, &PresetParams::cotton_kasp_hulse_kemp_2015(), None);
        assert_eq!(results.len(), 1);
        assert!(!results[0].errors.is_empty());
        assert!(results[0].primary.is_none());
    }

    #[test]
    fn run_batch_sends_progress_events() {
        let input = ">no_variant\nACGT\n";
        let (tx, rx) = std::sync::mpsc::channel();
        let _ = run_batch(input, &PresetParams::cotton_kasp_hulse_kemp_2015(), Some(tx));
        let events: Vec<_> = rx.try_iter().collect();
        assert!(matches!(events.first(), Some(ProgressEvent::Started { total_markers: 1 })));
        assert!(matches!(events.last(), Some(ProgressEvent::Complete)));
    }
}
