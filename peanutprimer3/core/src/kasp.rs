//! Standard KASP/PACE fluorescent tail sequences, applied to allele-specific
//! primers only when opted into. Tails are user-editable (see
//! [`crate::model::PresetParams`]) since exact tail sequences can vary by
//! vendor/assay.

/// Prepends `tail` to `primer_seq` (tails go on the 5' end).
pub fn apply_tail(primer_seq: &str, tail: &str) -> String {
    format!("{tail}{primer_seq}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{default_fam_tail, default_hex_tail};

    #[test]
    fn apply_tail_prepends() {
        assert_eq!(apply_tail("ACGT", "TAIL"), "TAILACGT");
    }

    #[test]
    fn default_tails_are_nonempty_and_distinct() {
        assert!(!default_fam_tail().is_empty());
        assert!(!default_hex_tail().is_empty());
        assert_ne!(default_fam_tail(), default_hex_tail());
    }
}
