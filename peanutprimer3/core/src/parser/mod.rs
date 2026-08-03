//! Flexible FASTA-like input parsing.
//!
//! Accepts classic FASTA (`>id description` header lines) but does not
//! require the `>` prefix: in its absence, a line is treated as a
//! header/description if it contains any character outside the allowed
//! "sequence alphabet" (IUPAC bases/ambiguity codes, digits, whitespace, and
//! the `[`, `]`, `/`, `-` bracket-variant symbols); otherwise it's treated as
//! more sequence for the current record. Blank lines always end the current
//! record.
//!
//! This heuristic is inherently fuzzy for the edge case of a short,
//! all-consonant description line that happens to use only IUPAC letters
//! (e.g. a bare accession like "RS12345" with no separators) — such a line
//! **will** be classified as sequence, not header, when no `>` is present.
//! This tradeoff is intentional and pinned down by the tests below rather
//! than left as an accident of implementation.

mod record;

pub use record::SequenceRecord;

/// Every character (case-insensitive) that's part of the accepted sequence
/// alphabet: the four bases, U, N, the ten IUPAC ambiguity codes, plus the
/// bracket-variant-notation symbols.
fn is_sequence_char(c: char) -> bool {
    if c.is_whitespace() || c.is_ascii_digit() {
        return true;
    }
    if matches!(c, '[' | ']' | '/' | '-') {
        return true;
    }
    matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'U' | 'N' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'V' | 'H' | 'D' | 'B')
}

fn looks_like_sequence(line: &str) -> bool {
    let trimmed = line.trim();
    !trimmed.is_empty() && trimmed.chars().all(is_sequence_char)
}

/// Strips whitespace (used to allow line-wrapped/blocked sequences like
/// `"Aaggagccca gcccatagaa"`), keeping only meaningful characters.
fn strip_sequence_whitespace(s: &str) -> String {
    s.chars().filter(|c| !c.is_whitespace()).collect()
}

/// Splits a header line into (id, description) on the first run of
/// whitespace, matching classic FASTA header conventions.
fn split_header(header: &str) -> (String, Option<String>) {
    match header.trim().split_once(char::is_whitespace) {
        Some((id, rest)) => {
            let rest = rest.trim();
            (id.to_string(), if rest.is_empty() { None } else { Some(rest.to_string()) })
        }
        None => (header.trim().to_string(), None),
    }
}

/// A pending record's header, tracking whether it came from an explicit
/// `>` line (in which case it gets id/description splitting, matching
/// classic FASTA) or was inferred from a header-like line with no `>`
/// (in which case the whole line becomes the id, unsplit, since there's no
/// reliable way to tell an id from a description without the marker).
enum Header {
    Explicit(String),
    Implicit(String),
}

fn finish_record(
    header: Option<Header>,
    seq_accum: String,
    auto_count: &mut usize,
) -> Option<SequenceRecord> {
    if seq_accum.is_empty() {
        return None;
    }
    let (id, description) = match header {
        Some(Header::Explicit(h)) => split_header(&h),
        Some(Header::Implicit(h)) => (h, None),
        None => {
            *auto_count += 1;
            (format!("Sequence_{auto_count}"), None)
        }
    };
    Some(SequenceRecord { id, description, sequence: seq_accum })
}

/// Parses input text into records. See module docs for the header-optional
/// heuristic.
pub fn parse_records(input: &str) -> Vec<SequenceRecord> {
    let mut records = Vec::new();
    let mut auto_count = 0usize;
    let mut cur_header: Option<Header> = None;
    let mut cur_seq = String::new();
    let mut have_current = false;

    for raw_line in input.lines() {
        let line = raw_line.trim_end_matches('\r');
        let trimmed = line.trim();

        if trimmed.is_empty() {
            if have_current {
                if let Some(rec) =
                    finish_record(cur_header.take(), std::mem::take(&mut cur_seq), &mut auto_count)
                {
                    records.push(rec);
                }
                have_current = false;
            }
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix('>') {
            if have_current
                && let Some(rec) =
                    finish_record(cur_header.take(), std::mem::take(&mut cur_seq), &mut auto_count)
                {
                    records.push(rec);
                }
            cur_header = Some(Header::Explicit(rest.trim().to_string()));
            cur_seq = String::new();
            have_current = true;
        } else if looks_like_sequence(trimmed) {
            if !have_current {
                have_current = true;
                cur_header = None;
            }
            cur_seq.push_str(&strip_sequence_whitespace(trimmed));
        } else {
            // Header-like line without a leading '>'.
            if have_current
                && let Some(rec) =
                    finish_record(cur_header.take(), std::mem::take(&mut cur_seq), &mut auto_count)
                {
                    records.push(rec);
                }
            cur_header = Some(Header::Implicit(trimmed.to_string()));
            cur_seq = String::new();
            have_current = true;
        }
    }

    if have_current
        && let Some(rec) = finish_record(cur_header.take(), std::mem::take(&mut cur_seq), &mut auto_count) {
            records.push(rec);
        }

    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classic_fasta_single_record() {
        let input = ">Chr1_12345 some description\nACGTACGTACGT\nACGTACGT\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id, "Chr1_12345");
        assert_eq!(recs[0].description.as_deref(), Some("some description"));
        assert_eq!(recs[0].sequence, "ACGTACGTACGTACGTACGT");
    }

    #[test]
    fn classic_fasta_multiple_records() {
        let input = ">a\nACGT\n>b desc\nTTTT\nGGGG\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "a");
        assert_eq!(recs[0].sequence, "ACGT");
        assert_eq!(recs[1].id, "b");
        assert_eq!(recs[1].sequence, "TTTTGGGG");
    }

    #[test]
    fn header_optional_blank_line_separated() {
        let input = "Marker one\nACGTRACGT\n\nMarker two\nTTTTYTTTT\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "Marker one");
        assert_eq!(recs[0].sequence, "ACGTRACGT");
        assert_eq!(recs[1].id, "Marker two");
        assert_eq!(recs[1].sequence, "TTTTYTTTT");
    }

    #[test]
    fn header_optional_no_header_at_all() {
        let input = "ACGTRACGT\n\nTTTTYTTTT\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "Sequence_1");
        assert_eq!(recs[1].id, "Sequence_2");
    }

    #[test]
    fn header_optional_description_line_with_underscore_is_header() {
        // "Chr_Position" style ID without '>' - underscore breaks the
        // sequence-alphabet check, so it's correctly classified as a header.
        let input = "Chr1_98765\nACGTRACGT\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id, "Chr1_98765");
        assert_eq!(recs[0].sequence, "ACGTRACGT");
    }

    #[test]
    fn mixed_header_presence_in_one_file() {
        let input = ">explicit_header\nACGT\nimplicit header line\nTTTT\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "explicit_header");
        assert_eq!(recs[1].id, "implicit header line");
    }

    #[test]
    fn adversarial_all_iupac_description_line_is_classified_as_sequence() {
        // Documented tie-break: a header-less line made up entirely of
        // IUPAC letters (no separators) is treated as sequence, even though
        // a human might have meant it as an id.
        let input = "RSKM\nACGT\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id, "Sequence_1");
        assert_eq!(recs[0].sequence, "RSKMACGT");
    }

    #[test]
    fn spaced_and_lowercase_blocks_like_legacy_dbsnp_export() {
        let input = ">gnl|dbSNP|rs1\nAaggagccca gcccatagaa tgatgtcttc\nR\nAGTATGATTG TGGGGGTGGG\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id, "gnl|dbSNP|rs1");
        assert_eq!(recs[0].sequence, "AaggagcccagcccatagaatgatgtcttcRAGTATGATTGTGGGGGTGGG");
    }

    #[test]
    fn bracket_notation_sequence_is_accepted() {
        let input = ">snp1\nACGTACGT[A/G]ACGTACGT\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].sequence, "ACGTACGT[A/G]ACGTACGT");
    }

    #[test]
    fn indel_bracket_notation_is_accepted() {
        let input = ">indel1\nACGTACGT[A/ATT]ACGTACGT\n>indel2\nACGTACGT[AA/-]ACGTACGT\n";
        let recs = parse_records(input);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].sequence, "ACGTACGT[A/ATT]ACGTACGT");
        assert_eq!(recs[1].sequence, "ACGTACGT[AA/-]ACGTACGT");
    }

    #[test]
    fn empty_input_yields_no_records() {
        assert!(parse_records("").is_empty());
        assert!(parse_records("   \n\n  \n").is_empty());
    }
}
