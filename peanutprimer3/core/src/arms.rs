//! ARMS second-mismatch injection (Little 1997), ported directly from the
//! legacy `replace_mismatch_base`/`add_second_mismatch` Perl functions
//! (`batchprimer3_results.cgi:3441-3520`). Opt-in, and scoped to true
//! biallelic SNPs only: the legacy function only has branches for the six
//! biallelic IUPAC codes (R/Y/S/W/K/M) and silently no-ops for anything
//! else (triallelic V/H/D/B, or by extension our indels) -- we make that
//! same restriction explicit and UI-visible instead of a silent no-op.
//!
//! The legacy logic, faithfully preserved: replace the base at the
//! mismatch position with `HASH[complement(original_base)]`, where `HASH`
//! is one of three fixed tables depending on the code group:
//! - R/Y: transition partner (A<->G, C<->T) -- `RY_hash` in the original.
//! - S/W: the complement itself, unchanged -- the original `SW_hash` is an
//!   identity map, so this collapses to a plain complement flip.
//! - K/M: transversion partner (A<->C, G<->T) -- `KM_hash` in the original.

use crate::variant::Variant;

fn complement(base: char) -> char {
    match base.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        other => other,
    }
}

/// R/Y-group second-mismatch table (transition partner).
fn ry_hash(b: char) -> char {
    match b {
        'G' => 'A',
        'T' => 'C',
        'A' => 'G',
        'C' => 'T',
        other => other,
    }
}

/// K/M-group second-mismatch table (transversion partner).
fn km_hash(b: char) -> char {
    match b {
        'G' => 'T',
        'T' => 'G',
        'A' => 'C',
        'C' => 'A',
        other => other,
    }
}

/// The IUPAC ambiguity code for an (unordered) pair of distinct bases, if
/// one of the six biallelic codes covers it. Used so this feature works
/// identically whether the variant came from an embedded IUPAC code or
/// bracket notation like `[A/G]` -- the destabilization rule only depends
/// on which two bases are being discriminated.
fn code_for_allele_pair(a: char, b: char) -> Option<char> {
    let mut pair = [a.to_ascii_uppercase(), b.to_ascii_uppercase()];
    pair.sort_unstable();
    match pair {
        ['A', 'G'] => Some('R'),
        ['C', 'T'] => Some('Y'),
        ['C', 'G'] => Some('S'),
        ['A', 'T'] => Some('W'),
        ['G', 'T'] => Some('K'),
        ['A', 'C'] => Some('M'),
        _ => None,
    }
}

/// Whether the ARMS second-mismatch feature applies to this variant at all:
/// exactly two single-base alleles forming one of the six biallelic codes.
pub fn supports_second_mismatch(variant: &Variant) -> bool {
    if !variant.is_snp_like() || variant.alleles.len() != 2 {
        return false;
    }
    let a = variant.alleles[0].chars().next();
    let b = variant.alleles[1].chars().next();
    matches!((a, b), (Some(a), Some(b)) if code_for_allele_pair(a, b).is_some())
}

/// Replaces the base at `mismatch_pos` (a negative offset from the 3' end,
/// e.g. -3 = 3 bases in from the end) with the ARMS-recommended
/// destabilizing base for the given IUPAC code. Returns the primer
/// unchanged if the position is out of range or the code isn't one of the
/// six biallelic codes this table covers.
pub fn inject_second_mismatch(primer_seq: &str, code: char, mismatch_pos: i32) -> String {
    let len = primer_seq.chars().count() as i64;
    let idx = if mismatch_pos < 0 { len + i64::from(mismatch_pos) } else { i64::from(mismatch_pos) };
    if idx < 0 || idx >= len {
        return primer_seq.to_string();
    }
    let idx = idx as usize;

    let mut chars: Vec<char> = primer_seq.chars().collect();
    let base = chars[idx].to_ascii_uppercase();
    let complemented = complement(base);

    let new_base = match code.to_ascii_uppercase() {
        'R' | 'Y' => ry_hash(complemented),
        'S' | 'W' => complemented,
        'K' | 'M' => km_hash(complemented),
        _ => return primer_seq.to_string(),
    };

    chars[idx] = new_base;
    chars.into_iter().collect()
}

/// Applies the second-mismatch injection to `primer_seq` if `variant`
/// supports it, given the two allele bases it discriminates. Returns the
/// (possibly unmodified) sequence and whether an injection actually
/// happened.
pub fn maybe_inject(variant: &Variant, primer_seq: &str, mismatch_pos: i32) -> (String, bool) {
    if !supports_second_mismatch(variant) {
        return (primer_seq.to_string(), false);
    }
    let a = variant.alleles[0].chars().next().unwrap();
    let b = variant.alleles[1].chars().next().unwrap();
    let code = code_for_allele_pair(a, b).unwrap();
    (inject_second_mismatch(primer_seq, code, mismatch_pos), true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::VariantSource;

    fn variant(alleles: &[&str], source: VariantSource) -> Variant {
        Variant {
            left_flank: "AAAA".to_string(),
            right_flank: "TTTT".to_string(),
            alleles: alleles.iter().map(|s| s.to_string()).collect(),
            source,
        }
    }

    #[test]
    fn code_for_allele_pair_matches_legacy_table() {
        assert_eq!(code_for_allele_pair('G', 'A'), Some('R'));
        assert_eq!(code_for_allele_pair('T', 'C'), Some('Y'));
        assert_eq!(code_for_allele_pair('G', 'C'), Some('S'));
        assert_eq!(code_for_allele_pair('A', 'T'), Some('W'));
        assert_eq!(code_for_allele_pair('G', 'T'), Some('K'));
        assert_eq!(code_for_allele_pair('A', 'C'), Some('M'));
    }

    #[test]
    fn inject_second_mismatch_ry_group() {
        // R/Y: new_base = ry_hash(complement(base)).
        // base='G' -> complement='C' -> ry_hash('C')='T'
        assert_eq!(inject_second_mismatch("AAAAG", 'R', -1), "AAAAT");
        // base='A' -> complement='T' -> ry_hash('T')='C'
        assert_eq!(inject_second_mismatch("AAAAA", 'Y', -1), "AAAAC");
    }

    #[test]
    fn inject_second_mismatch_sw_group_is_plain_complement() {
        // S/W: new_base = complement(base) directly (legacy SW_hash is identity).
        assert_eq!(inject_second_mismatch("AAAAG", 'S', -1), "AAAAC");
        assert_eq!(inject_second_mismatch("AAAAA", 'W', -1), "AAAAT");
    }

    #[test]
    fn inject_second_mismatch_km_group() {
        // K/M: new_base = km_hash(complement(base)).
        // base='G' -> complement='C' -> km_hash('C')='A'
        assert_eq!(inject_second_mismatch("AAAAG", 'K', -1), "AAAAA");
        // base='A' -> complement='T' -> km_hash('T')='G'
        assert_eq!(inject_second_mismatch("AAAAA", 'M', -1), "AAAAG");
    }

    #[test]
    fn inject_second_mismatch_respects_position_offset() {
        // 10-char primer, mismatch_pos=-3 -> index 7 (0-based).
        let primer = "AAAAAAAAAA";
        let result = inject_second_mismatch(primer, 'R', -3);
        assert_ne!(result.chars().nth(7), Some('A'));
        assert_eq!(result.chars().next(), Some('A'));
        assert_eq!(result.len(), primer.len());
    }

    #[test]
    fn inject_second_mismatch_out_of_range_is_noop() {
        assert_eq!(inject_second_mismatch("AAAA", 'R', -10), "AAAA");
    }

    #[test]
    fn supports_second_mismatch_true_for_iupac_snp() {
        let v = variant(&["G", "A"], VariantSource::Iupac('R'));
        assert!(supports_second_mismatch(&v));
    }

    #[test]
    fn supports_second_mismatch_true_for_bracket_snp() {
        let v = variant(&["A", "G"], VariantSource::Bracket);
        assert!(supports_second_mismatch(&v));
    }

    #[test]
    fn supports_second_mismatch_false_for_indel() {
        let v = variant(&["A", "ATT"], VariantSource::Bracket);
        assert!(!supports_second_mismatch(&v));
    }

    #[test]
    fn supports_second_mismatch_false_for_triallelic() {
        let v = variant(&["A", "C", "G"], VariantSource::Iupac('V'));
        assert!(!supports_second_mismatch(&v));
    }

    #[test]
    fn maybe_inject_noop_when_unsupported() {
        let v = variant(&["A", "ATT"], VariantSource::Bracket);
        let (seq, injected) = maybe_inject(&v, "AAAAG", -1);
        assert!(!injected);
        assert_eq!(seq, "AAAAG");
    }

    #[test]
    fn maybe_inject_applies_when_supported() {
        let v = variant(&["G", "A"], VariantSource::Iupac('R'));
        let (seq, injected) = maybe_inject(&v, "AAAAG", -1);
        assert!(injected);
        assert_eq!(seq, "AAAAT");
    }
}
