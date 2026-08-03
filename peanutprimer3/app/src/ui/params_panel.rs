use crate::state::{AppState, PresetChoice};

fn mark_custom_if_changed(state: &mut AppState, changed: bool) {
    if changed {
        state.preset_choice = PresetChoice::Custom;
    }
}

pub fn show(ui: &mut egui::Ui, state: &mut AppState) {
    ui.horizontal(|ui| {
        ui.label("Preset:");
        egui::ComboBox::from_id_salt("preset_choice")
            .selected_text(match state.preset_choice {
                PresetChoice::CottonKasp => "Cotton SNP KASP (Hulse-Kemp 2015)",
                PresetChoice::LegacyDefaults => "BatchPrimer3 legacy defaults",
                PresetChoice::Custom => "Custom",
            })
            .show_ui(ui, |ui| {
                if ui.selectable_label(state.preset_choice == PresetChoice::CottonKasp, "Cotton SNP KASP (Hulse-Kemp 2015)").clicked() {
                    state.apply_preset(PresetChoice::CottonKasp);
                }
                if ui.selectable_label(state.preset_choice == PresetChoice::LegacyDefaults, "BatchPrimer3 legacy defaults").clicked() {
                    state.apply_preset(PresetChoice::LegacyDefaults);
                }
            });
    });

    ui.separator();
    ui.heading("Primer size (bp)");
    ui.horizontal(|ui| {
        let mut changed = false;
        ui.label("Min");
        changed |= ui.add(egui::DragValue::new(&mut state.params.primer_min_size)).changed();
        ui.label("Opt");
        changed |= ui.add(egui::DragValue::new(&mut state.params.primer_opt_size)).changed();
        ui.label("Max");
        changed |= ui.add(egui::DragValue::new(&mut state.params.primer_max_size)).changed();
        mark_custom_if_changed(state, changed);
    });

    ui.heading("Melting temperature (C)");
    ui.horizontal(|ui| {
        let mut changed = false;
        ui.label("Min");
        changed |= ui.add(egui::DragValue::new(&mut state.params.primer_min_tm).speed(0.1)).changed();
        ui.label("Opt");
        changed |= ui.add(egui::DragValue::new(&mut state.params.primer_opt_tm).speed(0.1)).changed();
        ui.label("Max");
        changed |= ui.add(egui::DragValue::new(&mut state.params.primer_max_tm).speed(0.1)).changed();
        ui.label("Max diff between alleles");
        changed |= ui.add(egui::DragValue::new(&mut state.params.max_tm_difference).speed(0.1)).changed();
        mark_custom_if_changed(state, changed);
    });

    ui.heading("GC content (%)");
    ui.horizontal(|ui| {
        let mut changed = false;
        ui.label("Min");
        changed |= ui.add(egui::DragValue::new(&mut state.params.primer_min_gc).speed(0.5)).changed();
        ui.label("Max");
        changed |= ui.add(egui::DragValue::new(&mut state.params.primer_max_gc).speed(0.5)).changed();
        mark_custom_if_changed(state, changed);
    });

    ui.heading("Product size (common/flanking primer pair, bp)");
    ui.horizontal(|ui| {
        let mut changed = false;
        ui.label("Min");
        changed |= ui.add(egui::DragValue::new(&mut state.params.product_min_size)).changed();
        ui.label("Opt");
        changed |= ui.add(egui::DragValue::new(&mut state.params.product_opt_size)).changed();
        ui.label("Max");
        changed |= ui.add(egui::DragValue::new(&mut state.params.product_max_size)).changed();
        mark_custom_if_changed(state, changed);
    });

    ui.heading("Solution conditions");
    ui.horizontal(|ui| {
        let mut changed = false;
        ui.label("Salt (mM)");
        changed |= ui.add(egui::DragValue::new(&mut state.params.salt_conc_mm).speed(0.5)).changed();
        ui.label("DNA (nM)");
        changed |= ui.add(egui::DragValue::new(&mut state.params.dna_conc_nm).speed(0.5)).changed();
        mark_custom_if_changed(state, changed);
    });

    ui.heading("QC thresholds");
    ui.horizontal(|ui| {
        let mut changed = false;
        ui.label("Max self-complementarity Tm (C)");
        changed |= ui.add(egui::DragValue::new(&mut state.params.max_self_any).speed(0.5)).changed();
        ui.label("Max 3' self-complementarity Tm (C)");
        changed |= ui.add(egui::DragValue::new(&mut state.params.max_self_end).speed(0.5)).changed();
        ui.label("Max ambiguous bases");
        changed |= ui.add(egui::DragValue::new(&mut state.params.max_ambiguous_bases)).changed();
        mark_custom_if_changed(state, changed);
    });

    ui.separator();
    ui.heading("ARMS second mismatch (opt-in, biallelic SNPs only)");
    ui.add_enabled_ui(!state.has_non_snp_variant, |ui| {
        let mut changed = false;
        changed |= ui.checkbox(&mut state.params.second_mismatch_enabled, "Inject a deliberate second mismatch").changed();
        ui.horizontal(|ui| {
            ui.label("Position from 3' end");
            changed |= ui.add(egui::DragValue::new(&mut state.params.second_mismatch_pos)).changed();
        });
        mark_custom_if_changed(state, changed);
    });
    if state.has_non_snp_variant {
        ui.colored_label(
            egui::Color32::from_rgb(200, 120, 0),
            "Disabled: the loaded input contains an indel or triallelic variant, which this feature doesn't apply to.",
        );
    }

    ui.separator();
    ui.heading("KASP fluorescent tails (opt-in)");
    let mut changed = false;
    changed |= ui.checkbox(&mut state.params.kasp_tails_enabled, "Add standard KASP FAM/HEX tails to allele-specific primers").changed();
    ui.horizontal(|ui| {
        ui.label("FAM tail");
        changed |= ui.text_edit_singleline(&mut state.params.kasp_fam_tail).changed();
    });
    ui.horizontal(|ui| {
        ui.label("HEX tail");
        changed |= ui.text_edit_singleline(&mut state.params.kasp_hex_tail).changed();
    });
    mark_custom_if_changed(state, changed);

    ui.separator();
    ui.heading("Display options");
    ui.checkbox(
        &mut state.show_alternative_orientation,
        "Show the alternative orientation too, when one was found",
    );
}
