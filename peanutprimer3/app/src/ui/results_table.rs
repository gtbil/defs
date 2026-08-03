use egui_extras::{Column, TableBuilder, TableRow};
use peanutprimer3_core::model::MarkerTriplet;

use crate::state::AppState;

const ALT_ROW_COLOR: egui::Color32 = egui::Color32::from_rgb(150, 150, 150);

fn muted_or_plain(ui: &mut egui::Ui, text: impl Into<String>, muted: bool) {
    let text = text.into();
    if muted {
        ui.colored_label(ALT_ROW_COLOR, text);
    } else {
        ui.label(text);
    }
}

/// Renders the six data columns (Common/Product/Allele primers/Common
/// sequence/Tm diff OK/Notes) for one triplet, shared between the primary
/// and alternative rows. `extra_notes` carries row-specific context (the
/// flank warning and "alternative available" marker) that only makes sense
/// attached to the primary row.
fn render_triplet_columns(row: &mut TableRow, triplet: Option<&MarkerTriplet>, extra_notes: &[String], muted: bool) {
    match triplet {
        Some(t) => {
            row.col(|ui| muted_or_plain(ui, format!("{:?}", t.common.orientation), muted));
            row.col(|ui| muted_or_plain(ui, t.product_size.to_string(), muted));
            row.col(|ui| {
                let text = t
                    .allele_specific
                    .iter()
                    .map(|p| format!("{}: {}", p.allele, p.tailed_sequence.as_deref().unwrap_or(&p.sequence)))
                    .collect::<Vec<_>>()
                    .join("\n");
                muted_or_plain(ui, text, muted);
            });
            row.col(|ui| muted_or_plain(ui, t.common.sequence.clone(), muted));
            row.col(|ui| muted_or_plain(ui, if t.tm_balance_ok { "yes" } else { "no" }, muted));
            row.col(|ui| {
                let mut notes: Vec<String> = extra_notes.to_vec();
                notes.extend(t.heterodimer_warnings.iter().cloned());
                muted_or_plain(ui, notes.join(" | "), muted);
            });
        }
        None => {
            for _ in 0..4 {
                row.col(|ui| muted_or_plain(ui, "-", muted));
            }
            // A missing triplet only reaches this branch for the primary
            // row (the alternative is only ever rendered when `Some`), so
            // this is always the "no valid design found" error case --
            // always red, regardless of `muted`.
            row.col(|ui| {
                ui.colored_label(egui::Color32::RED, extra_notes.join(" | "));
            });
        }
    }
}

pub fn show(ui: &mut egui::Ui, state: &AppState) {
    if state.results.is_empty() {
        ui.label("No results yet -- load an input and click Run.");
        return;
    }

    TableBuilder::new(ui)
        .striped(true)
        .resizable(true)
        .column(Column::auto().at_least(100.0))
        .column(Column::auto().at_least(80.0))
        .column(Column::auto().at_least(70.0))
        .column(Column::remainder().at_least(180.0))
        .column(Column::remainder().at_least(180.0))
        .column(Column::auto().at_least(70.0))
        .column(Column::remainder().at_least(200.0))
        .header(22.0, |mut header| {
            header.col(|ui| { ui.strong("Sequence ID"); });
            header.col(|ui| { ui.strong("Common"); });
            header.col(|ui| { ui.strong("Product"); });
            header.col(|ui| { ui.strong("Allele-specific primers"); });
            header.col(|ui| { ui.strong("Common primer sequence"); });
            header.col(|ui| { ui.strong("Tm diff OK"); });
            header.col(|ui| { ui.strong("Notes"); });
        })
        .body(|mut body| {
            for result in &state.results {
                body.row(20.0, |mut row| {
                    row.col(|ui| { ui.label(&result.sequence_id); });

                    if result.primary.is_some() {
                        let mut notes = Vec::new();
                        if let Some(w) = &result.flank_warning {
                            notes.push(w.clone());
                        }
                        if result.alternative.is_some() {
                            notes.push("alternative orientation also available below".to_string());
                        }
                        render_triplet_columns(&mut row, result.primary.as_ref(), &notes, false);
                    } else {
                        render_triplet_columns(&mut row, None, &[result.errors.join("; ")], false);
                    }
                });

                if state.show_alternative_orientation
                    && let Some(alt) = &result.alternative {
                        body.row(20.0, |mut row| {
                            row.col(|ui| {
                                ui.colored_label(ALT_ROW_COLOR, format!("\u{21b3} {} (alternative)", result.sequence_id));
                            });
                            render_triplet_columns(&mut row, Some(alt), &[], true);
                        });
                    }
            }
        });
}
