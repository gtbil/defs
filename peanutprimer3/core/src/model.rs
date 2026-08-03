//! Shared domain types: presets, per-primer records, and per-marker results.

use serde::{Deserialize, Serialize};

/// Which strand/direction a primer reads in, relative to the input sequence
/// as given (not necessarily the genomic + strand).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Orientation {
    Forward,
    Reverse,
}

/// All tunable design parameters for a run. Two built-in presets are provided
/// (see [`PresetParams::cotton_kasp_hulse_kemp_2015`] and
/// [`PresetParams::batchprimer3_legacy_defaults`]); users may load/save their
/// own as JSON.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PresetParams {
    pub name: String,

    // Primer size (applies to both allele-specific and common primers)
    pub primer_min_size: usize,
    pub primer_opt_size: usize,
    pub primer_max_size: usize,

    // Melting temperature
    pub primer_min_tm: f64,
    pub primer_opt_tm: f64,
    pub primer_max_tm: f64,
    /// Maximum allowed Tm difference between the two allele-specific primers
    /// (this is our own QC check, not primer3's pair max_diff_tm, since the
    /// allele-specific primers are never picked as a primer3 "pair").
    pub max_tm_difference: f64,

    // GC content
    pub primer_min_gc: f64,
    pub primer_max_gc: f64,

    // Solution conditions
    pub salt_conc_mm: f64,
    pub dna_conc_nm: f64,

    // Complementarity / QC thresholds. These are thermodynamic structure
    // melting-temperature thresholds in Celsius (primer3 v2's `_th`
    // convention, e.g. `PRIMER_MAX_SELF_ANY_TH`, via the `primer3` crate's
    // `calc_homodimer`/`calc_end_stability`), not the old primer3 v1
    // alignment-score scale the legacy tool used -- the legacy numeric
    // defaults (8.0/3.0) don't carry over to this scale, so these default
    // to primer3 v2's own recommended value (47.0C) instead.
    pub max_self_any: f64,
    pub max_self_end: f64,
    pub max_ambiguous_bases: usize,

    // Product size (common/flanking primer pair)
    pub product_min_size: usize,
    pub product_opt_size: usize,
    pub product_max_size: usize,

    // ARMS second mismatch (SNP-only, opt-in)
    pub second_mismatch_enabled: bool,
    /// Position of the injected second mismatch, counted as a negative
    /// offset from the 3' terminal base (e.g. -3 = 3 bases in from the end).
    pub second_mismatch_pos: i32,

    // KASP fluorescent tails (opt-in)
    pub kasp_tails_enabled: bool,
    pub kasp_fam_tail: String,
    pub kasp_hex_tail: String,
}

impl PresetParams {
    /// Your working parameters (Hulse-Kemp et al. 2015, G3 5(6):1095-1105),
    /// set as the application default.
    pub fn cotton_kasp_hulse_kemp_2015() -> Self {
        Self {
            name: "Cotton SNP KASP (Hulse-Kemp 2015)".to_string(),
            primer_min_size: 15,
            primer_opt_size: 20,
            primer_max_size: 30,
            primer_min_tm: 55.0,
            primer_opt_tm: 57.0,
            primer_max_tm: 60.0,
            max_tm_difference: 2.0,
            primer_min_gc: 20.0,
            primer_max_gc: 80.0,
            salt_conc_mm: 50.0,
            dna_conc_nm: 50.0,
            max_self_any: 47.0,
            max_self_end: 47.0,
            max_ambiguous_bases: 0,
            product_min_size: 50,
            product_opt_size: 50,
            product_max_size: 100,
            second_mismatch_enabled: false,
            second_mismatch_pos: -3,
            kasp_tails_enabled: false,
            kasp_fam_tail: default_fam_tail(),
            kasp_hex_tail: default_hex_tail(),
        }
    }

    /// The original BatchPrimer3 tool's defaults, kept for reference/comparison.
    pub fn batchprimer3_legacy_defaults() -> Self {
        Self {
            name: "BatchPrimer3 legacy defaults".to_string(),
            primer_min_size: 15,
            primer_opt_size: 20,
            primer_max_size: 30,
            primer_min_tm: 50.0,
            primer_opt_tm: 60.0,
            primer_max_tm: 63.0,
            max_tm_difference: 10.0,
            primer_min_gc: 20.0,
            primer_max_gc: 80.0,
            salt_conc_mm: 50.0,
            dna_conc_nm: 50.0,
            max_self_any: 47.0,
            max_self_end: 47.0,
            max_ambiguous_bases: 0,
            // Confirmed from batchprimer3.cgi:208-210, 286-288 (generic
            // product-size form defaults, shown for mode 7 since it isn't in
            // the exclusion list) -- NOT the SNP_INNER_PRODUCT_* fields,
            // which are a red herring for mode 7 (only consumed by the
            // tetra-primer-ARMS-only design_ARMS_outer_primers path).
            product_min_size: 50,
            product_opt_size: 50,
            product_max_size: 100,
            second_mismatch_enabled: false,
            second_mismatch_pos: -3,
            kasp_tails_enabled: false,
            kasp_fam_tail: default_fam_tail(),
            kasp_hex_tail: default_hex_tail(),
        }
    }
}

impl Default for PresetParams {
    fn default() -> Self {
        Self::cotton_kasp_hulse_kemp_2015()
    }
}

/// Standard LGC KASP FAM-compatible tail.
pub fn default_fam_tail() -> String {
    "GAAGGTGACCAAGTTCATGCT".to_string()
}

/// Standard LGC KASP HEX-compatible tail.
pub fn default_hex_tail() -> String {
    "GAAGGTCGGAGTCAACGGATT".to_string()
}

/// One designed allele-specific primer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleSpecificPrimer {
    /// The allele this primer discriminates for (e.g. "A", "ATT", "" for a
    /// deletion allele).
    pub allele: String,
    pub orientation: Orientation,
    /// Untailed primer sequence (5'->3').
    pub sequence: String,
    /// `sequence` with the KASP FAM/HEX tail prepended, if enabled.
    pub tailed_sequence: Option<String>,
    /// Start offset relative to the variant's start (negative = upstream).
    pub pos_relative_to_variant: i64,
    pub length: usize,
    pub tm: f64,
    pub gc_content: f64,
    pub self_any: f64,
    pub self_end: f64,
    /// The 100-point composite QC score (0 = failed a hard cutoff).
    pub score: f64,
    pub has_second_mismatch: bool,
}

/// The shared/common flanking primer (designed by real primer3).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CommonPrimer {
    pub orientation: Orientation,
    pub sequence: String,
    pub tm: f64,
    pub gc_content: f64,
    pub self_any: f64,
    pub self_end: f64,
}

/// A single self-consistent, correctly-oriented triplet: one common primer
/// plus its matching allele-specific primers, grouped so no manual
/// re-pairing is ever needed.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MarkerTriplet {
    pub common: CommonPrimer,
    pub allele_specific: Vec<AlleleSpecificPrimer>,
    pub product_size: usize,
    pub tm_difference_between_alleles: f64,
    pub tm_balance_ok: bool,
    /// Human-readable notes on any heterodimer risk between the common
    /// primer and an allele-specific primer.
    pub heterodimer_warnings: Vec<String>,
    /// Combined score used to rank Option A vs Option B.
    pub combined_score: f64,
}

/// Full result for one input marker (one variant in one sequence record).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MarkerResult {
    pub sequence_id: String,
    pub variant_description: String,
    pub flank_warning: Option<String>,
    /// The auto-selected best triplet, if any valid one was found.
    pub primary: Option<MarkerTriplet>,
    /// The other orientation's triplet, kept for the detail view.
    pub alternative: Option<MarkerTriplet>,
    pub errors: Vec<String>,
}
