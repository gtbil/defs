//! Unified variant model: IUPAC ambiguity codes, bracket-notation SNPs, and
//! bracket-notation indels are all normalized into one `Variant` shape.

pub mod anchor;

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum VariantSource {
    /// A single IUPAC ambiguity code embedded in the sequence.
    Iupac(char),
    /// Bracket notation, e.g. `[A/G]`, `[A/ATT]`, `[AA/-]`.
    Bracket,
}

/// A variant site, represented as invariant flanking sequence plus a list of
/// allele strings (each possibly empty, for a deletion allele). There is no
/// shared coordinate space once allele lengths differ, so all downstream
/// logic works with `left_flank`/`right_flank` + per-allele substitution
/// rather than a single reference sequence with a point mutation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Variant {
    pub left_flank: String,
    pub right_flank: String,
    pub alleles: Vec<String>,
    pub source: VariantSource,
}

impl Variant {
    /// True only if every allele is exactly one base (a true SNP) — the
    /// only case the ARMS second-mismatch feature applies to.
    pub fn is_snp_like(&self) -> bool {
        self.alleles.iter().all(|a| a.len() == 1)
    }

    /// The full, allele-substituted sequence for allele `i`.
    pub fn full_seq(&self, allele_idx: usize) -> String {
        format!("{}{}{}", self.left_flank, self.alleles[allele_idx], self.right_flank)
    }

    /// The variant's start position (0-based, byte offset), identical for
    /// every allele since `left_flank` is invariant.
    pub fn variant_start(&self) -> usize {
        self.left_flank.len()
    }
}

/// Decodes a single IUPAC ambiguity code into its component bases. Returns
/// `None` for unambiguous bases (A/C/G/T/U/N) or non-DNA characters.
pub fn convert_code(c: char) -> Option<Vec<String>> {
    let alleles: &[&str] = match c.to_ascii_uppercase() {
        'S' => &["G", "C"],
        'W' => &["A", "T"],
        'R' => &["G", "A"],
        'Y' => &["T", "C"],
        'K' => &["G", "T"],
        'M' => &["A", "C"],
        'V' => &["A", "C", "G"],
        'H' => &["A", "C", "T"],
        'D' => &["A", "G", "T"],
        'B' => &["C", "G", "T"],
        _ => return None,
    };
    Some(alleles.iter().map(|s| (*s).to_string()).collect())
}

/// Scans a sequence for both bracket-notation variants (`[A/G]`,
/// `[A/ATT]`, `[AA/-]`) and embedded IUPAC ambiguity codes, returning all
/// variants found. Bracket regions are excluded from the IUPAC scan so a
/// base letter inside a bracket is never double-counted as its own variant.
pub fn find_variants(sequence: &str) -> Vec<Variant> {
    let mut variants = Vec::new();
    let mut bracket_ranges: Vec<(usize, usize)> = Vec::new();

    let mut i = 0;
    while i < sequence.len() {
        if sequence.as_bytes()[i] == b'['
            && let Some(close_rel) = sequence[i..].find(']') {
                let close = i + close_rel;
                let inner = &sequence[i + 1..close];
                let alleles: Vec<String> = inner
                    .split('/')
                    .map(|a| if a == "-" { String::new() } else { a.to_string() })
                    .collect();
                if alleles.len() >= 2 {
                    variants.push(Variant {
                        left_flank: sequence[..i].to_string(),
                        right_flank: sequence[close + 1..].to_string(),
                        alleles,
                        source: VariantSource::Bracket,
                    });
                }
                bracket_ranges.push((i, close));
                i = close + 1;
                continue;
            }
        i += 1;
    }

    for (idx, ch) in sequence.char_indices() {
        if bracket_ranges.iter().any(|&(s, e)| idx >= s && idx <= e) {
            continue;
        }
        if let Some(alleles) = convert_code(ch) {
            variants.push(Variant {
                left_flank: sequence[..idx].to_string(),
                right_flank: sequence[idx + ch.len_utf8()..].to_string(),
                alleles,
                source: VariantSource::Iupac(ch.to_ascii_uppercase()),
            });
        }
    }

    variants
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convert_code_biallelic() {
        assert_eq!(convert_code('R'), Some(vec!["G".to_string(), "A".to_string()]));
        assert_eq!(convert_code('Y'), Some(vec!["T".to_string(), "C".to_string()]));
        assert_eq!(convert_code('S'), Some(vec!["G".to_string(), "C".to_string()]));
        assert_eq!(convert_code('W'), Some(vec!["A".to_string(), "T".to_string()]));
        assert_eq!(convert_code('K'), Some(vec!["G".to_string(), "T".to_string()]));
        assert_eq!(convert_code('M'), Some(vec!["A".to_string(), "C".to_string()]));
    }

    #[test]
    fn convert_code_triallelic() {
        assert_eq!(convert_code('V'), Some(vec!["A".to_string(), "C".to_string(), "G".to_string()]));
        assert_eq!(convert_code('B'), Some(vec!["C".to_string(), "G".to_string(), "T".to_string()]));
    }

    #[test]
    fn convert_code_unambiguous_returns_none() {
        assert_eq!(convert_code('A'), None);
        assert_eq!(convert_code('N'), None);
    }

    #[test]
    fn find_variants_iupac_snp() {
        let vs = find_variants("AAAARCCCC");
        assert_eq!(vs.len(), 1);
        assert_eq!(vs[0].left_flank, "AAAA");
        assert_eq!(vs[0].right_flank, "CCCC");
        assert_eq!(vs[0].alleles, vec!["G".to_string(), "A".to_string()]);
        assert_eq!(vs[0].source, VariantSource::Iupac('R'));
        assert!(vs[0].is_snp_like());
    }

    #[test]
    fn find_variants_bracket_snp() {
        let vs = find_variants("AAAA[A/G]CCCC");
        assert_eq!(vs.len(), 1);
        assert_eq!(vs[0].left_flank, "AAAA");
        assert_eq!(vs[0].right_flank, "CCCC");
        assert_eq!(vs[0].alleles, vec!["A".to_string(), "G".to_string()]);
        assert_eq!(vs[0].source, VariantSource::Bracket);
        assert!(vs[0].is_snp_like());
    }

    #[test]
    fn find_variants_bracket_insertion() {
        let vs = find_variants("AAAA[A/ATT]CCCC");
        assert_eq!(vs.len(), 1);
        assert_eq!(vs[0].alleles, vec!["A".to_string(), "ATT".to_string()]);
        assert!(!vs[0].is_snp_like());
    }

    #[test]
    fn find_variants_bracket_deletion() {
        let vs = find_variants("AAAA[AA/-]CCCC");
        assert_eq!(vs.len(), 1);
        assert_eq!(vs[0].alleles, vec!["AA".to_string(), "".to_string()]);
        assert!(!vs[0].is_snp_like());
    }

    #[test]
    fn find_variants_bracket_does_not_double_count_as_iupac() {
        // The 'A' and 'G' inside the brackets must not also be scanned as
        // (non-ambiguous, so moot) IUPAC codes; more importantly, a
        // multi-allele IUPAC-looking letter placed *inside* brackets must
        // not produce a spurious second Variant.
        let vs = find_variants("AAAA[R/Y]CCCC");
        assert_eq!(vs.len(), 1);
        assert_eq!(vs[0].alleles, vec!["R".to_string(), "Y".to_string()]);
    }

    #[test]
    fn find_variants_none_found() {
        assert!(find_variants("ACGTACGTACGT").is_empty());
    }
}
