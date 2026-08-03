//! End-to-end integration tests: raw file input -> batch design -> export,
//! exercising the whole public pipeline the way the app actually drives it
//! (as opposed to the per-module unit tests that check individual pieces
//! in isolation).

use peanutprimer3_core::batch::run_batch;
use peanutprimer3_core::model::PresetParams;
use peanutprimer3_core::report;

fn fixture(name: &str) -> String {
    let path = format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read fixture {path}: {e}"))
}

#[test]
fn legacy_style_iupac_input_designs_a_full_triplet() {
    let input = fixture("legacy_style_iupac.fasta");
    let results = run_batch(&input, &PresetParams::cotton_kasp_hulse_kemp_2015(), None);

    assert_eq!(results.len(), 1);
    let marker = &results[0];
    assert!(marker.errors.is_empty(), "unexpected errors: {:?}", marker.errors);
    assert!(marker.flank_warning.is_none(), "150bp flanks should not trigger the short-flank warning");

    let triplet = marker.primary.as_ref().expect("expected a primary triplet");
    assert_eq!(triplet.allele_specific.len(), 2);
    let alleles: Vec<&str> = triplet.allele_specific.iter().map(|p| p.allele.as_str()).collect();
    assert!(alleles.contains(&"G"));
    assert!(alleles.contains(&"A"));
    // Both allele-specific primers share one orientation; the common
    // primer takes the opposite one -- the core "no manual re-pairing"
    // guarantee.
    assert_eq!(triplet.allele_specific[0].orientation, triplet.allele_specific[1].orientation);
    assert_ne!(triplet.allele_specific[0].orientation, triplet.common.orientation);
}

#[test]
fn cotton_bracket_notation_input_with_chr_position_header_designs_a_triplet() {
    let input = fixture("cotton_bracket_notation.fasta");
    let results = run_batch(&input, &PresetParams::cotton_kasp_hulse_kemp_2015(), None);

    assert_eq!(results.len(), 1);
    assert_eq!(results[0].sequence_id, "Chr1_45678");
    assert!(results[0].errors.is_empty(), "unexpected errors: {:?}", results[0].errors);
    let triplet = results[0].primary.as_ref().expect("expected a primary triplet");
    assert!(triplet.product_size >= 50 && triplet.product_size <= 100);
}

#[test]
fn short_flank_input_still_designs_but_warns() {
    let input = fixture("short_flank_warning.fasta");
    // Relax product size range since there isn't 50-100bp of room on a
    // 20bp flank; this test is specifically about the warning, not about
    // whether the strict Hulse-Kemp product-size window is achievable here.
    let mut params = PresetParams::cotton_kasp_hulse_kemp_2015();
    params.product_min_size = 10;
    params.product_max_size = 60;

    let results = run_batch(&input, &params, None);
    assert_eq!(results.len(), 1);
    assert!(results[0].flank_warning.is_some(), "expected a short-flank warning");
}

#[test]
fn batch_results_export_to_csv_and_xlsx_without_error() {
    let input = fixture("legacy_style_iupac.fasta");
    let results = run_batch(&input, &PresetParams::cotton_kasp_hulse_kemp_2015(), None);

    let csv = report::write_csv_string(&results, false).expect("csv export should succeed");
    assert!(csv.contains("Sequence ID"));
    assert!(csv.contains("gnl|dbSNP|rs1234567"));

    let xlsx = report::write_xlsx_bytes(&results, false).expect("xlsx export should succeed");
    assert_eq!(&xlsx[0..2], b"PK");
}
