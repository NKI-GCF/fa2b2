// Roel Kluin, 2023, GPL v3

use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::marker::KmerIter;
use crate::new_types::{extended_position::ExtPosEtc, position::Position};
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use crate::xmerhasher::XmerHasher;
use anyhow::Result;
use crossbeam_channel::{unbounded, Receiver, Sender};
use itertools::izip;
use noodles_fasta as fasta;
use std::io::BufRead;
use std::iter::repeat;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;
use std::thread::spawn;
use std::thread::JoinHandle;

pub(crate) struct MarkContigThreads {
    pub(crate) threads: Vec<JoinHandle<Result<(), anyhow::Error>>>,
    pub(crate) rx_in_main: Receiver<(Vec<ExtPosEtc>, Vec<XmerLoc>)>,
    tx_to_thread: Vec<Sender<XmerLoc>>,
}

fn gen_xloc_threads(ct: usize) -> impl Iterator<Item = (Sender<XmerLoc>, Receiver<XmerLoc>)> {
    repeat(()).take(ct).map(|_| unbounded::<XmerLoc>())
}

impl MarkContigThreads {
    pub(crate) fn new(ct: usize, kmerlen: usize, rep_max_dist: Position) -> Self {
        // channels to send XmerLocs to threads
        let (tx_to_thread, rx_from_main): (Vec<_>, Vec<_>) = gen_xloc_threads(ct).unzip();

        // channels for threads to send XmerLocs to one another when extended past their responsibility
        let (mut tx_inter_thread, rx_inter_thread): (Vec<_>, Vec<_>) = gen_xloc_threads(ct).unzip();

        // rotate transmitter back so that its corresponding receiver is passed to the next thread.
        let first = tx_inter_thread.remove(0);
        tx_inter_thread.push(first);

        // to collect the results
        let (tx_to_main, rx_in_main) = unbounded::<(Vec<ExtPosEtc>, Vec<XmerLoc>)>();
        let shutdown_poll = Arc::new(AtomicUsize::new(0));

        let threads: Vec<_> = izip!(
            0..ct,
            rx_from_main.into_iter(),
            tx_inter_thread.into_iter(),
            rx_inter_thread.into_iter(),
        )
        .map(|xmer_channels| {
            let tx_to_main = tx_to_main.clone();
            let shutdown_poll = shutdown_poll.clone();
            let mut xh = XmerHasher::new(ct, kmerlen, rep_max_dist, xmer_channels, tx_to_main);
            spawn(move || xh.work(shutdown_poll))
        })
        .collect();
        MarkContigThreads {
            threads,
            rx_in_main,
            tx_to_thread,
        }
    }
    pub(crate) fn mark_contig<T>(
        &self,
        ks: &mut KmerStore,
        kc: KmerConst,
        mut fa: fasta::Reader<T>,
    ) -> Result<()>
    where
        T: BufRead,
    {
        // The main thread to read from fasta and prefilter
        let mut kmi = KmerIter::new(ks, &kc)?;

        for res in fa.records() {
            let record = res?;
            dbg_print!("Starting with record {}.", record.name());
            kmi.markcontig(&self.tx_to_thread, &kc, record)?;
            dbg_print!("Finished with record.");
        }
        dbg_print!("Ending transmission to threads.");
        Ok(())
    }
}
