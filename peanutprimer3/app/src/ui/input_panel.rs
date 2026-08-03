use crate::state::AppState;

pub fn show(ui: &mut egui::Ui, state: &mut AppState) {
    ui.horizontal(|ui| {
        if ui.button("Load FASTA file...").clicked()
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("FASTA/text", &["fasta", "fa", "txt"])
                .pick_file()
            {
                match std::fs::read_to_string(&path) {
                    Ok(contents) => {
                        state.input_text = contents;
                        state.refresh_snp_check();
                        state.status = format!("Loaded {}", path.display());
                    }
                    Err(e) => state.status = format!("Failed to read file: {e}"),
                }
            }
        if ui.button("Clear").clicked() {
            state.input_text.clear();
            state.refresh_snp_check();
        }
    });

    ui.label(
        "Paste one or more marker sequences below. Each SNP/indel can be marked either as an \
         embedded IUPAC ambiguity code (e.g. R for G/A) or bracket notation (e.g. [A/G], \
         [A/ATT], [AA/-]). A leading '>' header is optional -- blank lines and header-like \
         lines are used to split records automatically.",
    );

    let response = ui.add(
        egui::TextEdit::multiline(&mut state.input_text)
            .desired_rows(20)
            .desired_width(f32::INFINITY)
            .font(egui::TextStyle::Monospace),
    );
    if response.changed() {
        state.refresh_snp_check();
    }
}
