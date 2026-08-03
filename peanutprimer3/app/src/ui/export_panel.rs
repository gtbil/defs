use peanutprimer3_core::report;

use crate::state::AppState;

pub fn show(ui: &mut egui::Ui, state: &mut AppState) {
    if state.results.is_empty() {
        ui.label("Run a batch first to enable export.");
        return;
    }

    let include_alternative = state.show_alternative_orientation;

    ui.horizontal(|ui| {
        if ui.button("Export CSV...").clicked() {
            export_with(state, "csv", |results| {
                report::write_csv_string(results, include_alternative).map(String::into_bytes)
            });
        }
        if ui.button("Export TSV...").clicked() {
            export_with(state, "tsv", |results| {
                report::write_tsv_string(results, include_alternative).map(String::into_bytes)
            });
        }
        if ui.button("Export Excel (.xlsx)...").clicked() {
            export_with(state, "xlsx", move |results| report::write_xlsx_bytes(results, include_alternative));
        }
    });
}

fn export_with(
    state: &mut AppState,
    extension: &str,
    build: impl FnOnce(&[peanutprimer3_core::model::MarkerResult]) -> Result<Vec<u8>, String>,
) {
    let Some(path) = rfd::FileDialog::new()
        .set_file_name(format!("peanutprimer3_results.{extension}"))
        .save_file()
    else {
        return;
    };

    match build(&state.results).and_then(|bytes| std::fs::write(&path, bytes).map_err(|e| e.to_string())) {
        Ok(()) => state.status = format!("Exported to {}", path.display()),
        Err(e) => state.status = format!("Export failed: {e}"),
    }
}
