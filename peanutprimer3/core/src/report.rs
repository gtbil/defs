//! CSV/TSV and Excel export. Both formats use the same grouped-per-marker
//! layout (common primer + its matching allele-specific primers together)
//! rather than the legacy tool's separate, unordered listing -- so there's
//! never a reason to manually re-pair rows after export. The alternative
//! (non-selected) orientation, when one was found, can optionally be
//! included as a second, `"Alt "`-prefixed block per marker.

use csv::WriterBuilder;
use rust_xlsxwriter::Workbook;

use crate::model::{MarkerResult, MarkerTriplet};

/// The IUPAC codes this tool supports cover at most 3 alleles (V/H/D/B),
/// so the report reserves a fixed 3-allele-wide block of columns.
const MAX_ALLELES_IN_REPORT: usize = 3;

fn fmt_f64(v: f64) -> String {
    format!("{v:.2}")
}

/// Header columns for one triplet block (common primer + product/Tm-diff +
/// up to `MAX_ALLELES_IN_REPORT` allele-specific primers), `prefix`ed so the
/// same shape can be reused for the primary and (optionally) alternative
/// blocks without duplicating column names.
fn triplet_header(prefix: &str) -> Vec<String> {
    let mut headers: Vec<String> = [
        "Common Orientation",
        "Common Sequence",
        "Common Tm",
        "Common GC%",
        "Product Size",
        "Tm Difference",
        "Tm Balance OK",
        "Heterodimer Warnings",
    ]
    .into_iter()
    .map(|s| format!("{prefix}{s}"))
    .collect();

    for i in 1..=MAX_ALLELES_IN_REPORT {
        headers.extend([
            format!("{prefix}Allele{i} Base"),
            format!("{prefix}Allele{i} Orientation"),
            format!("{prefix}Allele{i} Primer"),
            format!("{prefix}Allele{i} Tailed Primer"),
            format!("{prefix}Allele{i} Tm"),
            format!("{prefix}Allele{i} GC%"),
            format!("{prefix}Allele{i} Score"),
            format!("{prefix}Allele{i} 2nd Mismatch"),
        ]);
    }
    headers
}

fn header_row(include_alternative: bool) -> Vec<String> {
    let mut headers: Vec<String> =
        ["Sequence ID", "Variant", "Flank Warning", "Errors", "Alternative Available"]
            .into_iter()
            .map(String::from)
            .collect();
    headers.extend(triplet_header(""));
    if include_alternative {
        headers.extend(triplet_header("Alt "));
    }
    headers
}

/// Values for one triplet block, in the same order as [`triplet_header`].
/// Blank for every column when `triplet` is `None` (no valid design in that
/// orientation).
fn triplet_values(triplet: Option<&MarkerTriplet>) -> Vec<String> {
    let mut row = vec![
        triplet.map(|t| format!("{:?}", t.common.orientation)).unwrap_or_default(),
        triplet.map(|t| t.common.sequence.clone()).unwrap_or_default(),
        triplet.map(|t| fmt_f64(t.common.tm)).unwrap_or_default(),
        triplet.map(|t| fmt_f64(t.common.gc_content)).unwrap_or_default(),
        triplet.map(|t| t.product_size.to_string()).unwrap_or_default(),
        triplet.map(|t| fmt_f64(t.tm_difference_between_alleles)).unwrap_or_default(),
        triplet.map(|t| t.tm_balance_ok.to_string()).unwrap_or_default(),
        triplet.map(|t| t.heterodimer_warnings.join("; ")).unwrap_or_default(),
    ];

    for i in 0..MAX_ALLELES_IN_REPORT {
        match triplet.and_then(|t| t.allele_specific.get(i)) {
            Some(p) => {
                row.push(p.allele.clone());
                row.push(format!("{:?}", p.orientation));
                row.push(p.sequence.clone());
                row.push(p.tailed_sequence.clone().unwrap_or_default());
                row.push(fmt_f64(p.tm));
                row.push(fmt_f64(p.gc_content));
                row.push(fmt_f64(p.score));
                row.push(p.has_second_mismatch.to_string());
            }
            None => row.extend(std::iter::repeat_n(String::new(), 8)),
        }
    }

    row
}

fn build_row(result: &MarkerResult, include_alternative: bool) -> Vec<String> {
    let mut row = vec![
        result.sequence_id.clone(),
        result.variant_description.clone(),
        result.flank_warning.clone().unwrap_or_default(),
        result.errors.join("; "),
        result.alternative.is_some().to_string(),
    ];
    row.extend(triplet_values(result.primary.as_ref()));
    if include_alternative {
        row.extend(triplet_values(result.alternative.as_ref()));
    }
    row
}

fn write_delimited(results: &[MarkerResult], delimiter: u8, include_alternative: bool) -> Result<String, String> {
    let mut wtr = WriterBuilder::new().delimiter(delimiter).from_writer(vec![]);
    wtr.write_record(header_row(include_alternative)).map_err(|e| e.to_string())?;
    for r in results {
        wtr.write_record(build_row(r, include_alternative)).map_err(|e| e.to_string())?;
    }
    let bytes = wtr.into_inner().map_err(|e| e.to_string())?;
    String::from_utf8(bytes).map_err(|e| e.to_string())
}

pub fn write_csv_string(results: &[MarkerResult], include_alternative: bool) -> Result<String, String> {
    write_delimited(results, b',', include_alternative)
}

pub fn write_tsv_string(results: &[MarkerResult], include_alternative: bool) -> Result<String, String> {
    write_delimited(results, b'\t', include_alternative)
}

pub fn write_xlsx_bytes(results: &[MarkerResult], include_alternative: bool) -> Result<Vec<u8>, String> {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();

    for (col, header) in header_row(include_alternative).iter().enumerate() {
        worksheet.write(0, col as u16, header).map_err(|e| e.to_string())?;
    }
    for (row_idx, result) in results.iter().enumerate() {
        for (col, value) in build_row(result, include_alternative).iter().enumerate() {
            worksheet.write((row_idx + 1) as u32, col as u16, value).map_err(|e| e.to_string())?;
        }
    }

    workbook.save_to_buffer().map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{AlleleSpecificPrimer, CommonPrimer, MarkerTriplet, Orientation};

    fn sample_triplet(common_tm: f64) -> MarkerTriplet {
        MarkerTriplet {
            common: CommonPrimer {
                orientation: Orientation::Reverse,
                sequence: "ACGTACGTACGTACGTACGT".to_string(),
                tm: common_tm,
                gc_content: 50.0,
                self_any: 0.0,
                self_end: 0.0,
            },
            allele_specific: vec![
                AlleleSpecificPrimer {
                    allele: "G".to_string(),
                    orientation: Orientation::Forward,
                    sequence: "GGGGCATCGATCG".to_string(),
                    tailed_sequence: None,
                    pos_relative_to_variant: -10,
                    length: 13,
                    tm: 58.0,
                    gc_content: 55.0,
                    self_any: 10.0,
                    self_end: 5.0,
                    score: 90.0,
                    has_second_mismatch: false,
                },
                AlleleSpecificPrimer {
                    allele: "A".to_string(),
                    orientation: Orientation::Forward,
                    sequence: "AGGGCATCGATCG".to_string(),
                    tailed_sequence: None,
                    pos_relative_to_variant: -10,
                    length: 13,
                    tm: 57.5,
                    gc_content: 50.0,
                    self_any: 8.0,
                    self_end: 4.0,
                    score: 88.0,
                    has_second_mismatch: false,
                },
            ],
            product_size: 75,
            tm_difference_between_alleles: 0.5,
            tm_balance_ok: true,
            heterodimer_warnings: vec![],
            combined_score: 178.0,
        }
    }

    fn sample_result() -> MarkerResult {
        MarkerResult {
            sequence_id: "marker_1".to_string(),
            variant_description: "G/A (variant start 100)".to_string(),
            flank_warning: None,
            primary: Some(sample_triplet(59.0)),
            alternative: None,
            errors: vec![],
        }
    }

    #[test]
    fn csv_export_has_header_and_one_data_row() {
        let csv = write_csv_string(&[sample_result()], false).unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("Sequence ID,Variant,"));
        assert!(lines[1].starts_with("marker_1,"));
        assert!(lines[1].contains("GGGGCATCGATCG"));
        assert!(lines[1].contains("AGGGCATCGATCG"));
    }

    #[test]
    fn tsv_export_uses_tab_delimiter() {
        let tsv = write_tsv_string(&[sample_result()], false).unwrap();
        let first_line = tsv.lines().next().unwrap();
        assert!(first_line.contains('\t'));
        assert!(!first_line.contains(','));
    }

    #[test]
    fn xlsx_export_produces_nonempty_bytes() {
        let bytes = write_xlsx_bytes(&[sample_result()], false).unwrap();
        assert!(!bytes.is_empty());
        // .xlsx files are zip archives; check the zip magic bytes.
        assert_eq!(&bytes[0..2], b"PK");
    }

    #[test]
    fn missing_primary_produces_blank_primer_columns() {
        let mut result = sample_result();
        result.primary = None;
        result.errors = vec!["no valid triplet".to_string()];
        let csv = write_csv_string(&[result], false).unwrap();
        let data_line = csv.lines().nth(1).unwrap();
        assert!(data_line.contains("no valid triplet"));
    }

    #[test]
    fn excluding_alternative_omits_alt_columns_even_when_present() {
        let mut result = sample_result();
        result.alternative = Some(sample_triplet(61.0));
        let csv = write_csv_string(&[result], false).unwrap();
        let header = csv.lines().next().unwrap();
        assert!(!header.contains("Alt "));
        // The "Alternative Available" marker-level column still reflects it.
        let data_line = csv.lines().nth(1).unwrap();
        let alt_available_idx = header.split(',').position(|c| c == "Alternative Available").unwrap();
        assert_eq!(data_line.split(',').nth(alt_available_idx).unwrap(), "true");
    }

    #[test]
    fn including_alternative_adds_alt_prefixed_columns_with_its_data() {
        let mut result = sample_result();
        result.alternative = Some(sample_triplet(61.0));
        let csv = write_csv_string(&[result], true).unwrap();
        let header = csv.lines().next().unwrap();
        assert!(header.contains("Alt Common Orientation"));
        assert!(header.contains("Alt Allele1 Base"));

        let data_line = csv.lines().nth(1).unwrap();
        let alt_tm_idx = header.split(',').position(|c| c == "Alt Common Tm").unwrap();
        assert_eq!(data_line.split(',').nth(alt_tm_idx).unwrap(), "61.00");
    }

    #[test]
    fn including_alternative_blanks_alt_columns_when_none_present() {
        let result = sample_result(); // alternative: None
        let csv = write_csv_string(&[result], true).unwrap();
        let header = csv.lines().next().unwrap();
        let data_line = csv.lines().nth(1).unwrap();
        let alt_tm_idx = header.split(',').position(|c| c == "Alt Common Tm").unwrap();
        assert_eq!(data_line.split(',').nth(alt_tm_idx).unwrap(), "");
    }
}
