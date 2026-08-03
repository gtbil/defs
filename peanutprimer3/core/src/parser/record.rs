//! A single parsed input record: an identifier plus its sequence (which may
//! contain embedded IUPAC ambiguity codes and/or bracket-notation variants).

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    pub id: String,
    pub description: Option<String>,
    /// Raw sequence text, case preserved, whitespace/digits stripped.
    pub sequence: String,
}
