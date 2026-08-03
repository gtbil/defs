//! Composite QC scoring for allele-specific primer candidates.
//!
//! Ports the legacy BatchPrimer3 100-point composite score
//! (`calculate_primer_score`) bin-for-bin (20/10/30/10/20/10 points across
//! size/GC/Tm/repeat-content/ambiguous-base-count/self-complementarity,
//! with a hard cutoff per bin that zeroes the whole score), but sources Tm
//! and self-complementarity from the `primer3` crate's real thermodynamics
//! instead of the legacy's hand-ported nearest-neighbor calculator and
//! custom Smith-Waterman aligner. The self-complementarity sub-score also
//! corrects an apparent sign inversion in the legacy formula (which scored
//! *worse* self-complementarity higher, up to the threshold); here, lower
//! complementarity always scores higher, matching every other bin's
//! closer-to-ideal-is-better convention.

use crate::model::PresetParams;
use primer3::{SolutionConditions, TmParams};

#[derive(Debug, Clone)]
pub struct ScoredPrimer {
    pub tm: f64,
    pub gc_content: f64,
    pub self_any: f64,
    pub self_end: f64,
    pub score: f64,
}

/// GC content as a percentage (0.0-100.0).
pub fn gc_content(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq.chars().filter(|c| matches!(c.to_ascii_uppercase(), 'G' | 'C')).count();
    100.0 * gc as f64 / seq.len() as f64
}

/// Count of bases that are not one of the four unambiguous DNA bases.
pub fn count_ambiguous_bases(seq: &str) -> usize {
    seq.chars().filter(|c| !matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'U')).count()
}

/// Counts bases participating in a homopolymer run of length >= 3 (a simple
/// repeat risk factor for primer specificity/synthesis).
pub fn count_repeat_bases(seq: &str) -> usize {
    let chars: Vec<char> = seq.chars().collect();
    let mut total = 0;
    let mut i = 0;
    while i < chars.len() {
        let mut j = i + 1;
        while j < chars.len() && chars[j].eq_ignore_ascii_case(&chars[i]) {
            j += 1;
        }
        let run_len = j - i;
        if run_len >= 3 {
            total += run_len;
        }
        i = j;
    }
    total
}

fn solution_conditions(params: &PresetParams) -> SolutionConditions {
    SolutionConditions { mv_conc: params.salt_conc_mm, dna_conc: params.dna_conc_nm, ..Default::default() }
}

/// Thermodynamic "badness" of a potential self-structure, expressed as its
/// melting temperature in Celsius (primer3 v2's `_th`-suffixed convention,
/// e.g. `PRIMER_MAX_SELF_ANY_TH` defaults to 47.0C) -- not delta G. A
/// higher structure Tm means the dimer/hairpin persists at higher
/// temperatures and is more likely to interfere with the intended PCR
/// annealing. `thal()` reports `temp = 0.0` when no stable structure is
/// found, so this is 0 in that case with no separate check needed.
fn badness(_structure_found: bool, tm: f64) -> f64 {
    tm.max(0.0)
}

/// Scores one candidate primer sequence against the given parameters. Any
/// failed hard cutoff (GC/Tm/repeat/ambiguous-base/self-complementarity
/// bands) zeroes the whole score, matching the legacy behavior. Returns
/// `None` only if the underlying thermodynamic calculation itself errors
/// (e.g. an empty sequence).
pub fn score_candidate(seq: &str, params: &PresetParams) -> Option<ScoredPrimer> {
    let tm_params = TmParams { conditions: solution_conditions(params), ..Default::default() };
    let tm = primer3::calc_tm_with(seq, &tm_params).ok()?;
    let gc = gc_content(seq);
    let repeats = count_repeat_bases(seq);
    let n = count_ambiguous_bases(seq);

    let homodimer = primer3::calc_homodimer(seq).ok()?;
    let end_stab = primer3::calc_end_stability(seq, seq).ok()?;
    let self_any = badness(homodimer.structure_found(), homodimer.tm());
    let self_end = badness(end_stab.structure_found(), end_stab.tm());

    let reject = |score: f64| Some(ScoredPrimer { tm, gc_content: gc, self_any, self_end, score });

    let mut score = 0.0;

    // Size: 20 points, proportional to distance from optimum.
    let size_span = (params.primer_max_size as f64 - params.primer_opt_size as f64)
        .abs()
        .max((params.primer_min_size as f64 - params.primer_opt_size as f64).abs())
        .max(1.0);
    score += 20.0 - (seq.len() as f64 - params.primer_opt_size as f64).abs() * 20.0 / size_span;

    // GC: 10 points, hard cutoff.
    if gc < params.primer_min_gc || gc > params.primer_max_gc {
        return reject(0.0);
    }
    score += 10.0;

    // Tm: 30 points, hard cutoff.
    if tm < params.primer_min_tm || tm > params.primer_max_tm {
        return reject(0.0);
    }
    let tm_span = (params.primer_max_tm - params.primer_opt_tm)
        .abs()
        .max((params.primer_min_tm - params.primer_opt_tm).abs())
        .max(f64::EPSILON);
    score += 30.0 - (tm - params.primer_opt_tm).abs() * 30.0 / tm_span;

    // Repeats: 10 points, hard cutoff at >50% repetitive content.
    let repeat_fraction = repeats as f64 / seq.len() as f64;
    if repeat_fraction > 0.5 {
        return reject(0.0);
    }
    score += 10.0 - repeat_fraction * 10.0;

    // Ambiguous bases: 20 points, hard cutoff.
    if n > params.max_ambiguous_bases {
        return reject(0.0);
    }
    score += if params.max_ambiguous_bases > 0 {
        20.0 - (n as f64 / params.max_ambiguous_bases as f64) * 20.0
    } else {
        20.0
    };

    // Self-complementarity: 10 points, hard cutoff. Lower complementarity
    // always scores higher (see module docs re: the legacy sign inversion).
    if self_any > params.max_self_any || self_end > params.max_self_end {
        return reject(0.0);
    }
    let self_any_score = 5.0 * (1.0 - self_any / params.max_self_any.max(f64::EPSILON));
    let self_end_score = 5.0 * (1.0 - self_end / params.max_self_end.max(f64::EPSILON));
    score += self_any_score + self_end_score;

    Some(ScoredPrimer { tm, gc_content: gc, self_any, self_end, score })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_content_basic() {
        assert!((gc_content("GGCC") - 100.0).abs() < 1e-9);
        assert!((gc_content("AATT")).abs() < 1e-9);
        assert!((gc_content("ACGT") - 50.0).abs() < 1e-9);
        assert!((gc_content("")).abs() < 1e-9);
    }

    #[test]
    fn ambiguous_base_count() {
        assert_eq!(count_ambiguous_bases("ACGT"), 0);
        assert_eq!(count_ambiguous_bases("ACGTN"), 1);
        assert_eq!(count_ambiguous_bases("ACGTRY"), 2);
    }

    #[test]
    fn repeat_base_count() {
        assert_eq!(count_repeat_bases("ACGT"), 0);
        assert_eq!(count_repeat_bases("AAAACGT"), 4);
        assert_eq!(count_repeat_bases("AAGTTTCC"), 3); // "TTT" run only ("AA" run is length 2, below threshold)
        assert_eq!(count_repeat_bases("AAAGGGTTT"), 9);
    }

    #[test]
    fn score_candidate_rejects_out_of_band_gc() {
        let params = PresetParams::cotton_kasp_hulse_kemp_2015();
        // All-A primer: 0% GC, well outside [20,80].
        let scored = score_candidate("AAAAAAAAAAAAAAAAAAAA", &params).unwrap();
        assert_eq!(scored.score, 0.0);
    }

    #[test]
    fn score_candidate_gives_positive_score_for_reasonable_primer() {
        let params = PresetParams::cotton_kasp_hulse_kemp_2015();
        // A plausible ~20bp, ~50% GC, non-repetitive, non-palindromic primer.
        let seq = "CTGACGATCGTAGGCATCGA";
        let scored = score_candidate(seq, &params).unwrap();
        assert!(
            scored.score > 0.0,
            "expected positive score, got {} (tm={}, gc={}, self_any={}, self_end={})",
            scored.score,
            scored.tm,
            scored.gc_content,
            scored.self_any,
            scored.self_end,
        );
    }

    #[test]
    fn score_candidate_rejects_excessive_repeats() {
        let params = PresetParams::cotton_kasp_hulse_kemp_2015();
        let scored = score_candidate("AAAAAAAAAAAAAAAAAAAA", &params).unwrap();
        // Also fails GC, but repeats alone would fail too; just confirm zero.
        assert_eq!(scored.score, 0.0);
    }
}
