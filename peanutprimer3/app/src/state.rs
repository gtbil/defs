use std::sync::mpsc::Receiver;

use peanutprimer3_core::model::{MarkerResult, PresetParams};

use crate::worker::WorkerEvent;

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum Tab {
    Input,
    Parameters,
    Results,
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum PresetChoice {
    CottonKasp,
    LegacyDefaults,
    Custom,
}

pub struct AppState {
    pub input_text: String,
    pub params: PresetParams,
    pub preset_choice: PresetChoice,
    pub tab: Tab,

    pub results: Vec<MarkerResult>,
    pub status: String,
    pub running: bool,
    pub total_markers: usize,
    pub completed_markers: usize,
    pub worker_rx: Option<Receiver<WorkerEvent>>,

    /// True when any variant in the currently-loaded input isn't a
    /// biallelic SNP -- used to disable the ARMS second-mismatch toggle
    /// (rather than silently ignoring it) per the domain's constraint that
    /// the feature only applies to true biallelic SNPs.
    pub has_non_snp_variant: bool,

    /// Whether to show/export the alternative (non-selected) orientation's
    /// triplet alongside the primary one, when one was found. A display/
    /// export preference, not a `PresetParams` field -- it doesn't affect
    /// primer design at all.
    pub show_alternative_orientation: bool,
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            input_text: String::new(),
            params: PresetParams::default(),
            preset_choice: PresetChoice::CottonKasp,
            tab: Tab::Input,
            results: Vec::new(),
            status: String::new(),
            running: false,
            total_markers: 0,
            completed_markers: 0,
            worker_rx: None,
            has_non_snp_variant: false,
            show_alternative_orientation: false,
        }
    }
}

impl AppState {
    pub fn apply_preset(&mut self, choice: PresetChoice) {
        self.preset_choice = choice;
        match choice {
            PresetChoice::CottonKasp => self.params = PresetParams::cotton_kasp_hulse_kemp_2015(),
            PresetChoice::LegacyDefaults => self.params = PresetParams::batchprimer3_legacy_defaults(),
            PresetChoice::Custom => {}
        }
    }

    /// Re-scans the current input for any non-biallelic-SNP variant, used
    /// to gate the ARMS toggle in the parameters panel.
    pub fn refresh_snp_check(&mut self) {
        use peanutprimer3_core::arms::supports_second_mismatch;
        use peanutprimer3_core::parser::parse_records;
        use peanutprimer3_core::variant::find_variants;

        self.has_non_snp_variant = parse_records(&self.input_text).iter().any(|record| {
            let variants = find_variants(&record.sequence);
            variants.iter().any(|v| !supports_second_mismatch(v))
        });
    }

    /// Poll the worker channel for new events. Returns true if a repaint
    /// should be requested (a job is still running).
    pub fn poll_worker(&mut self) -> bool {
        let Some(rx) = self.worker_rx.take() else { return false };
        let mut still_running = self.running;
        let mut done = false;
        while let Ok(event) = rx.try_recv() {
            match event {
                WorkerEvent::Progress(peanutprimer3_core::batch::ProgressEvent::Started { total_markers }) => {
                    self.total_markers = total_markers;
                    self.completed_markers = 0;
                }
                WorkerEvent::Progress(peanutprimer3_core::batch::ProgressEvent::MarkerDone { .. }) => {
                    self.completed_markers += 1;
                }
                WorkerEvent::Progress(peanutprimer3_core::batch::ProgressEvent::Complete) => {}
                WorkerEvent::Finished(results) => {
                    self.results = results;
                    self.status = format!("Done: {} marker(s) processed", self.results.len());
                    still_running = false;
                    done = true;
                }
            }
        }
        if !done {
            self.worker_rx = Some(rx);
        }
        self.running = still_running;
        self.running
    }
}
