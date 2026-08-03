//! Common/flanking primer design, and the orientation-pairing logic that
//! automates away BatchPrimer3's manual "re-pair the flanking primer with
//! the matching allele-specific primers" step.
//!
//! There are exactly two valid triplet layouts per marker:
//! - **Option A**: allele-specific primers use the *forward*-anchored
//!   candidates (3' end at the variant) -> the common primer must be a
//!   *reverse* primer designed downstream of the variant.
//! - **Option B**: allele-specific primers use the *reverse-complement*-
//!   anchored candidates -> the common primer must be a *forward* primer
//!   designed upstream of the variant.
//!
//! [`design_marker`] builds both, auto-selects the better-scoring valid one
//! as `primary`, and keeps the other as `alternative`.

use std::ops::Range;

use primer3::{PrimerSettings, PrimerTask, SequenceArgs};

use crate::arms;
use crate::kasp;
use crate::model::{AlleleSpecificPrimer, CommonPrimer, MarkerResult, MarkerTriplet, Orientation, PresetParams};
use crate::scoring;
use crate::variant::anchor;
use crate::variant::Variant;

/// Minimum recommended flanking sequence length on each side of the
/// variant, per the working protocol (>=100bp -> >=101bp total).
const MIN_RECOMMENDED_FLANK: usize = 100;

/// A heterodimer is flagged as a warning when its structure's melting
/// temperature exceeds this threshold -- an arbitrary but reasonable bar
/// for "could plausibly interfere near typical PCR annealing conditions".
const HETERODIMER_WARNING_TM_THRESHOLD: f64 = 40.0;

struct DesignedCommonPrimer {
    primer: CommonPrimer,
    /// Position on the template used to design it (`variant.full_seq(0)`).
    position_on_allele0_template: Range<usize>,
    penalty: f64,
}

/// Designs candidate common/flanking primers on the given orientation,
/// sorted by primer3's own penalty (best first). The `PrimerTask::PickPrimerList`
/// task deliberately skips pairing (so it can return single left/right
/// primers rather than requiring both), but that also means primer3's
/// product-size-range logic -- which only applies during pairing -- never
/// kicks in here. Two things compensate for that: (1) the excluded region
/// covers the *entire* wrong side of the variant (not just the variant
/// itself), so a "left primer" can't be positioned somewhere irrelevant in
/// the right flank, and (2) [`build_triplet`] scans this candidate list
/// itself for the first one whose resulting product size actually falls in
/// `[product_min_size, product_max_size]`.
fn design_common_primer_candidates(
    variant: &Variant,
    common_orientation: Orientation,
    params: &PresetParams,
) -> Result<Vec<DesignedCommonPrimer>, String> {
    let template = variant.full_seq(0);
    let variant_start = variant.variant_start();
    let max_allele_len = variant.alleles.iter().map(String::len).max().unwrap_or(1).max(1);
    let variant_end = (variant_start + max_allele_len).min(template.len());

    // Exclude everything on the wrong side of the variant (plus the
    // variant itself), so the only region left to search is the correct
    // flank for this orientation.
    let excluded = match common_orientation {
        Orientation::Forward => (variant_start, template.len() - variant_start),
        Orientation::Reverse => (0, variant_end),
    };

    let seq_args = SequenceArgs::builder()
        .sequence(template)
        .excluded_region(excluded.0, excluded.1)
        .build()
        .map_err(|e| e.to_string())?;

    let settings = PrimerSettings::builder()
        .task(PrimerTask::PickPrimerList)
        .pick_left_primer(matches!(common_orientation, Orientation::Forward))
        .pick_right_primer(matches!(common_orientation, Orientation::Reverse))
        .primer_opt_size(params.primer_opt_size)
        .primer_min_size(params.primer_min_size)
        .primer_max_size(params.primer_max_size)
        .primer_opt_tm(params.primer_opt_tm)
        .primer_min_tm(params.primer_min_tm)
        .primer_max_tm(params.primer_max_tm)
        .primer_min_gc(params.primer_min_gc)
        .primer_max_gc(params.primer_max_gc)
        .primer_mv_conc(params.salt_conc_mm)
        .primer_dna_conc(params.dna_conc_nm)
        .product_size_range(params.product_min_size, params.product_max_size)
        .num_return(50)
        .build()
        .map_err(|e| e.to_string())?;

    let result = primer3::design_primers(&seq_args, &settings, None, None).map_err(|e| e.to_string())?;

    let candidates: &[primer3::PrimerRecord] = match common_orientation {
        Orientation::Forward => result.left_primers(),
        Orientation::Reverse => result.right_primers(),
    };

    let mut designed: Vec<DesignedCommonPrimer> = candidates
        .iter()
        .map(|rec| DesignedCommonPrimer {
            primer: CommonPrimer {
                orientation: common_orientation,
                sequence: rec.sequence().to_string(),
                tm: rec.tm(),
                gc_content: rec.gc_content(),
                self_any: rec.self_any(),
                self_end: rec.self_end(),
            },
            position_on_allele0_template: rec.position_on_template(),
            penalty: rec.penalty(),
        })
        .collect();
    designed.sort_by(|a, b| a.penalty.partial_cmp(&b.penalty).unwrap_or(std::cmp::Ordering::Equal));

    Ok(designed)
}

/// Computes the product size a given common-primer candidate would yield,
/// for the given allele-specific orientation, using allele 0's geometry
/// (see [`build_triplet`] for the indel caveat).
fn product_size_for_candidate(
    variant: &Variant,
    allele_specific_orientation: Orientation,
    allele0_primer: &AlleleSpecificPrimer,
    allele0_anchor: usize,
    common: &DesignedCommonPrimer,
) -> i64 {
    match allele_specific_orientation {
        Orientation::Forward => {
            let as_start = variant.variant_start() as i64 + allele0_primer.pos_relative_to_variant;
            common.position_on_allele0_template.end as i64 - as_start
        }
        Orientation::Reverse => {
            let as_end = allele0_anchor + allele0_primer.length;
            as_end as i64 - common.position_on_allele0_template.start as i64
        }
    }
}

/// For each allele, finds the single best-scoring candidate in the given
/// orientation. Returns `None` if any allele has no candidate that passes
/// its hard QC cutoffs -- a triplet needs every allele covered.
fn best_candidates_for_orientation(
    variant: &Variant,
    orientation: Orientation,
    params: &PresetParams,
) -> Option<Vec<AlleleSpecificPrimer>> {
    let mut out = Vec::with_capacity(variant.alleles.len());
    for allele_idx in 0..variant.alleles.len() {
        let candidates =
            anchor::generate_candidates(variant, allele_idx, params.primer_min_size, params.primer_max_size);
        let mut best: Option<AlleleSpecificPrimer> = None;

        for cand in candidates.into_iter().filter(|c| c.orientation == orientation) {
            let (final_seq, has_mismatch) = if params.second_mismatch_enabled {
                arms::maybe_inject(variant, &cand.sequence, params.second_mismatch_pos)
            } else {
                (cand.sequence.clone(), false)
            };

            let Some(scored) = scoring::score_candidate(&final_seq, params) else { continue };
            if scored.score <= 0.0 {
                continue;
            }
            if best.as_ref().is_some_and(|b| scored.score <= b.score) {
                continue;
            }

            let tailed_sequence = params.kasp_tails_enabled.then(|| match allele_idx {
                0 => Some(kasp::apply_tail(&final_seq, &params.kasp_fam_tail)),
                1 => Some(kasp::apply_tail(&final_seq, &params.kasp_hex_tail)),
                _ => None,
            }).flatten();

            best = Some(AlleleSpecificPrimer {
                allele: variant.alleles[allele_idx].clone(),
                orientation,
                length: final_seq.len(),
                sequence: final_seq,
                tailed_sequence,
                pos_relative_to_variant: cand.pos_relative_to_variant,
                tm: scored.tm,
                gc_content: scored.gc_content,
                self_any: scored.self_any,
                self_end: scored.self_end,
                score: scored.score,
                has_second_mismatch: has_mismatch,
            });
        }

        out.push(best?);
    }
    Some(out)
}

fn build_triplet(
    variant: &Variant,
    allele_specific_orientation: Orientation,
    params: &PresetParams,
) -> Option<MarkerTriplet> {
    let common_orientation = match allele_specific_orientation {
        Orientation::Forward => Orientation::Reverse,
        Orientation::Reverse => Orientation::Forward,
    };

    let allele_primers = best_candidates_for_orientation(variant, allele_specific_orientation, params)?;
    let common_candidates = design_common_primer_candidates(variant, common_orientation, params).ok()?;

    // Product size is computed using allele 0's geometry; for indels where
    // allele lengths differ, the true product size shifts by that length
    // difference per allele -- a documented approximation, since the
    // primary use case (biallelic SNPs) has identical allele lengths.
    let allele0_anchor = anchor::anchor_index(variant.variant_start(), variant.alleles[0].len());

    // Candidates are already sorted best-penalty-first; take the first one
    // whose resulting product size actually falls in the configured range
    // (PickPrimerList doesn't enforce this itself -- see the doc comment
    // on `design_common_primer_candidates`).
    let (designed_common, product_size) = common_candidates.into_iter().find_map(|c| {
        let size = product_size_for_candidate(variant, allele_specific_orientation, &allele_primers[0], allele0_anchor, &c);
        (size >= params.product_min_size as i64 && size <= params.product_max_size as i64)
            .then_some((c, size as usize))
    })?;

    let tms: Vec<f64> = allele_primers.iter().map(|p| p.tm).collect();
    let tm_difference_between_alleles =
        tms.iter().copied().fold(f64::MIN, f64::max) - tms.iter().copied().fold(f64::MAX, f64::min);
    let tm_balance_ok = tm_difference_between_alleles <= params.max_tm_difference;

    let mut heterodimer_warnings = Vec::new();
    for p in &allele_primers {
        let seq_to_check = p.tailed_sequence.as_deref().unwrap_or(&p.sequence);
        if let Ok(dimer) = primer3::calc_heterodimer(seq_to_check, &designed_common.primer.sequence)
            && dimer.structure_found() && dimer.tm() > HETERODIMER_WARNING_TM_THRESHOLD {
                heterodimer_warnings.push(format!(
                    "Possible heterodimer between the '{}' allele-specific primer and the common primer (structure Tm={:.1}C)",
                    p.allele, dimer.tm()
                ));
            }
    }

    let mut combined_score: f64 = allele_primers.iter().map(|p| p.score).sum();
    if !tm_balance_ok {
        combined_score -= 20.0;
    }
    combined_score -= 10.0 * heterodimer_warnings.len() as f64;

    Some(MarkerTriplet {
        common: designed_common.primer,
        allele_specific: allele_primers,
        product_size,
        tm_difference_between_alleles,
        tm_balance_ok,
        heterodimer_warnings,
        combined_score,
    })
}

fn flank_warning(variant: &Variant) -> Option<String> {
    if variant.left_flank.len() < MIN_RECOMMENDED_FLANK || variant.right_flank.len() < MIN_RECOMMENDED_FLANK {
        Some(format!(
            "Flanking sequence shorter than the recommended {}bp minimum on each side (left={}bp, right={}bp)",
            MIN_RECOMMENDED_FLANK,
            variant.left_flank.len(),
            variant.right_flank.len(),
        ))
    } else {
        None
    }
}

fn describe_variant(variant: &Variant) -> String {
    format!("{} (variant start {})", variant.alleles.join("/"), variant.variant_start())
}

/// Designs both orientation options for one marker, auto-selects the
/// better-scoring valid one as `primary`, and keeps the other (if also
/// valid) as `alternative` for the detail view.
pub fn design_marker(sequence_id: &str, variant: &Variant, params: &PresetParams) -> MarkerResult {
    let mut errors = Vec::new();

    if variant.alleles.len() < 2 {
        errors.push("Variant must have at least two alleles".to_string());
        return MarkerResult {
            sequence_id: sequence_id.to_string(),
            variant_description: describe_variant(variant),
            flank_warning: flank_warning(variant),
            primary: None,
            alternative: None,
            errors,
        };
    }

    let option_a = build_triplet(variant, Orientation::Forward, params);
    let option_b = build_triplet(variant, Orientation::Reverse, params);

    let (primary, alternative) = match (option_a, option_b) {
        (Some(a), Some(b)) if a.combined_score >= b.combined_score => (Some(a), Some(b)),
        (Some(a), Some(b)) => (Some(b), Some(a)),
        (Some(a), None) => (Some(a), None),
        (None, Some(b)) => (Some(b), None),
        (None, None) => {
            errors.push("No valid primer triplet found in either orientation".to_string());
            (None, None)
        }
    };

    MarkerResult {
        sequence_id: sequence_id.to_string(),
        variant_description: describe_variant(variant),
        flank_warning: flank_warning(variant),
        primary,
        alternative,
        errors,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::VariantSource;

    fn relaxed_params() -> PresetParams {
        let mut p = PresetParams::cotton_kasp_hulse_kemp_2015();
        // Widen constraints so the test isn't at the mercy of one specific
        // hand-written sequence happening to satisfy a narrow window.
        p.primer_min_tm = 45.0;
        p.primer_max_tm = 75.0;
        p.product_min_size = 50;
        p.product_opt_size = 100;
        p.product_max_size = 400;
        p
    }

    fn repeat_to_len(unit: &str, len: usize) -> String {
        unit.chars().cycle().take(len).collect()
    }

    #[test]
    fn design_marker_finds_a_valid_triplet_for_a_realistic_snp() {
        let left = repeat_to_len("CTGACGATCGTAGGCATCGA", 150);
        let right = repeat_to_len("TAGCCTGAACGGTCATGCAT", 150);
        let variant = Variant {
            left_flank: left,
            right_flank: right,
            alleles: vec!["G".to_string(), "A".to_string()],
            source: VariantSource::Iupac('R'),
        };

        let result = design_marker("test_marker", &variant, &relaxed_params());
        assert!(result.errors.is_empty(), "unexpected errors: {:?}", result.errors);
        let primary = result.primary.expect("expected a primary triplet to be found");
        assert_eq!(primary.allele_specific.len(), 2);
        // The two allele-specific primers must be in the SAME orientation
        // as each other, and the common primer in the OPPOSITE orientation.
        assert_eq!(primary.allele_specific[0].orientation, primary.allele_specific[1].orientation);
        assert_ne!(primary.allele_specific[0].orientation, primary.common.orientation);
        assert!(primary.product_size > 0);
    }

    #[test]
    fn design_marker_finds_a_valid_triplet_with_strict_hulse_kemp_defaults() {
        // Same shape as the relaxed test, but using the actual shipped
        // preset (Tm 55-60, product 50/50/100bp) with no relaxation, to
        // confirm the real working parameters aren't merely a theoretical
        // preset that never actually finds anything in practice.
        let left = repeat_to_len("CTGACGATCGTAGGCATCGA", 150);
        let right = repeat_to_len("TAGCCTGAACGGTCATGCAT", 150);
        let variant = Variant {
            left_flank: left,
            right_flank: right,
            alleles: vec!["G".to_string(), "A".to_string()],
            source: VariantSource::Iupac('R'),
        };

        let result = design_marker("test_marker", &variant, &PresetParams::cotton_kasp_hulse_kemp_2015());
        assert!(result.errors.is_empty(), "unexpected errors: {:?}", result.errors);
        let primary = result.primary.expect("expected a primary triplet with strict Hulse-Kemp defaults");
        assert!(primary.product_size >= 50 && primary.product_size <= 100);
    }

    #[test]
    fn design_marker_reports_flank_warning_for_short_flanks() {
        let variant = Variant {
            left_flank: "ACGT".to_string(),
            right_flank: "ACGT".to_string(),
            alleles: vec!["G".to_string(), "A".to_string()],
            source: VariantSource::Iupac('R'),
        };
        let result = design_marker("short", &variant, &relaxed_params());
        assert!(result.flank_warning.is_some());
    }

    #[test]
    fn design_marker_errors_on_single_allele_variant() {
        let variant = Variant {
            left_flank: "A".repeat(150),
            right_flank: "A".repeat(150),
            alleles: vec!["G".to_string()],
            source: VariantSource::Bracket,
        };
        let result = design_marker("bad", &variant, &relaxed_params());
        assert!(!result.errors.is_empty());
        assert!(result.primary.is_none());
    }
}
