// Suppresses the console window on Windows release builds; harmless on
// other platforms (the attribute is a no-op there).
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod state;
mod ui;
mod worker;

use state::{AppState, Tab};

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1000.0, 700.0]),
        // Glow (OpenGL) instead of the wgpu default: this is a simple
        // forms/table UI with no need for wgpu's modern-GPU features, and
        // OpenGL is far more universally supported -- notably including
        // Wine, where wgpu's DX12 backend (translated via vkd3d) is one of
        // the less mature paths and can hang during shader compilation.
        renderer: eframe::Renderer::Glow,
        ..Default::default()
    };
    eframe::run_native(
        "PeanutPrimer3 - KASP/PACE allele-specific primer designer",
        options,
        Box::new(|_cc| Ok(Box::new(PeanutPrimer3App::default()))),
    )
}

#[derive(Default)]
struct PeanutPrimer3App {
    state: AppState,
}

impl eframe::App for PeanutPrimer3App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        if self.state.poll_worker() {
            ui.ctx().request_repaint();
        }

        egui::Panel::top("top_panel").show(ui, |ui| {
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.state.tab, Tab::Input, "Input");
                ui.selectable_value(&mut self.state.tab, Tab::Parameters, "Parameters");
                ui.selectable_value(&mut self.state.tab, Tab::Results, "Results");

                ui.separator();

                let run_enabled = !self.state.running && !self.state.input_text.trim().is_empty();
                if ui.add_enabled(run_enabled, egui::Button::new("Run")).clicked() {
                    self.state.running = true;
                    self.state.results.clear();
                    self.state.status = "Running...".to_string();
                    self.state.worker_rx =
                        Some(worker::spawn_batch_job(self.state.input_text.clone(), self.state.params.clone()));
                }

                if self.state.running {
                    ui.add(egui::widgets::Spinner::new());
                    if self.state.total_markers > 0 {
                        ui.add(egui::ProgressBar::new(
                            self.state.completed_markers as f32 / self.state.total_markers as f32,
                        ).text(format!("{}/{}", self.state.completed_markers, self.state.total_markers)));
                    }
                }

                ui.separator();
                ui.label(&self.state.status);
            });
            ui.add_space(4.0);
        });

        egui::CentralPanel::default().show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| match self.state.tab {
                Tab::Input => ui::input_panel::show(ui, &mut self.state),
                Tab::Parameters => ui::params_panel::show(ui, &mut self.state),
                Tab::Results => {
                    ui::export_panel::show(ui, &mut self.state);
                    ui.separator();
                    ui::results_table::show(ui, &self.state);
                }
            });
        });
    }
}
