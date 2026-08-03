//! Runs a batch design job on a background thread so the egui UI thread
//! never blocks, forwarding progress events and the final result back over
//! a channel the UI polls non-blockingly each frame.

use std::sync::mpsc::{channel, Receiver};
use std::thread;

use peanutprimer3_core::batch::{run_batch, ProgressEvent};
use peanutprimer3_core::model::{MarkerResult, PresetParams};

pub enum WorkerEvent {
    Progress(ProgressEvent),
    Finished(Vec<MarkerResult>),
}

pub fn spawn_batch_job(input: String, params: PresetParams) -> Receiver<WorkerEvent> {
    let (out_tx, out_rx) = channel::<WorkerEvent>();
    let (progress_tx, progress_rx) = channel::<ProgressEvent>();

    // Forward progress events into the same outward channel as they arrive.
    {
        let out_tx = out_tx.clone();
        thread::spawn(move || {
            for ev in progress_rx {
                if out_tx.send(WorkerEvent::Progress(ev)).is_err() {
                    break;
                }
            }
        });
    }

    thread::spawn(move || {
        let results = run_batch(&input, &params, Some(progress_tx));
        let _ = out_tx.send(WorkerEvent::Finished(results));
    });

    out_rx
}
