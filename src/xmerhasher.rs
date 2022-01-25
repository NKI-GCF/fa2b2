use crate::kmerconst::XmerHash;
use crate::new_types::extended_position::ExtPosEtc;
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use anyhow::{ensure, Result};
use crossbeam_channel::{Receiver, Sender, TryRecvError};
use std::collections::VecDeque;
use std::sync::Arc;

pub(crate) struct XmerHasher {
    thread_nr: usize,
    thread_max: usize,
    kmerlen: usize,
    overbit: usize,
    rx_from_main: Receiver<XmerLoc>,
    tx_inter_thread: Sender<XmerLoc>,
    rx_inter_thread: Receiver<XmerLoc>,
    tx_to_main: Sender<(Vec<ExtPosEtc>, Vec<XmerLoc>)>,
    to_next: VecDeque<XmerLoc>,
    shutdown_poll: Arc<u64>,
}

impl XmerHasher {
    pub(crate) fn new(
        thread_nr: usize,
        thread_max: usize,
        kmerlen: usize,
        rx_from_main: Receiver<XmerLoc>,
        tx_inter_thread: Sender<XmerLoc>,
        rx_inter_thread: Receiver<XmerLoc>,
        tx_to_main: Sender<(Vec<ExtPosEtc>, Vec<XmerLoc>)>,
        shutdown_poll: Arc<u64>,
    ) -> Result<XmerHasher> {
        ensure!(kmerlen < 16);
        ensure!(thread_max.is_power_of_two());

        Ok(XmerHasher {
            thread_nr,
            thread_max,
            kmerlen,
            overbit: 1 << (kmerlen * 2 + 1),
            rx_from_main,
            tx_inter_thread,
            rx_inter_thread,
            tx_to_main,
            to_next: VecDeque::<XmerLoc>::with_capacity(1 << 16),
            shutdown_poll,
        })
    }
    pub(crate) fn work(&mut self) -> Result<()> {
        let bitlen = self.kmerlen * 2;
        let shl = bitlen - 2 - usize::try_from(self.thread_max.trailing_zeros()).unwrap();
        let mut kmp = vec![ExtPosEtc::default(); 1 << shl]; // position + strand per xmer handled in this thread.
        let mut max_extended = Vec::with_capacity(1 << 8);
        loop {
            // non blocking first
            match self.rx_from_main.try_recv().or_else(|e| {
                if e == TryRecvError::Empty {
                    self.rx_inter_thread.try_recv()
                } else {
                    Err(e)
                }
            }) {
                Ok(mark) => self.store_mark_and_extend(mark, &mut max_extended, &mut kmp)?,
                Err(TryRecvError::Disconnected) => break,
                Err(TryRecvError::Empty) => {
                    if let Some(err) = self
                        .to_next
                        .pop_front()
                        .and_then(|m| self.tx_inter_thread.try_send(m).err())
                    {
                        self.to_next.push_front(err.into_inner());
                    } else if let Ok(mark) = self.rx_from_main.recv() {
                        // We seem starved, block until more work.
                        self.store_mark_and_extend(mark, &mut max_extended, &mut kmp)?;
                    }
                }
            }
        }
        // The main channel has shutdown.
        loop {
            if let Ok(mark) = self.rx_inter_thread.try_recv() {
                if let Some(sp) = Arc::get_mut(&mut self.shutdown_poll) {
                    *sp &= !(1 << self.thread_nr);
                }
                self.store_mark_and_extend(mark, &mut max_extended, &mut kmp)?;
            } else if self.to_next.is_empty() {
                if let Some(sp) = Arc::get_mut(&mut self.shutdown_poll) {
                    *sp |= 1 << self.thread_nr;
                    if *sp == (1 << self.thread_max) - 1 {
                        break;
                    }
                }
            } else if let Some(err) = self
                .to_next
                .pop_front()
                .and_then(|m| self.tx_inter_thread.try_send(m).err())
            {
                self.to_next.push_front(err.into_inner());
            }
        }
        // send the packages back to main in order of threads.
        if self.thread_nr == 0 {
            self.tx_to_main.send((kmp, max_extended))?;
        } else if let Ok(_) = self.rx_inter_thread.recv() {
            self.tx_to_main.send((kmp, max_extended))?;
        }
        if self.thread_nr + 1 != self.thread_max {
            // message the next thread it can start sending.
            self.tx_inter_thread.send(XmerLoc::default())?;
        }
        Ok(())
    }
    /// The main work for the thread, hashing and extending index and updating its p, on collision.
    fn store_mark_and_extend(
        &mut self,
        mut mark: XmerLoc,
        max_extended: &mut Vec<XmerLoc>,
        kmp: &mut Vec<ExtPosEtc>,
    ) -> Result<()> {
        while self.try_store_mark(&mut mark, kmp)? {
            // From the XmerHash trait, extension fails on last, when all ext bits are set.
            if self.extend_xmer(&mut mark).is_err() {
                max_extended.push(mark);
                break;
            }
            if self.is_for_next_thread(&mark) {
                self.to_next.push_back(mark);
                break;
            }
        }
        Ok(())
    }
    /// the top bits of the xmerhash indicate whether the extension can be stored in this kmp
    /// or the kmp of the next thread.
    pub(crate) fn is_for_next_thread(&self, mark: &XmerLoc) -> bool {
        let kmer = self.unhash_and_uncompress_to_kmer(mark.idx, mark.p.x());
        kmer & (self.thread_max - 1) != self.thread_nr
    }
    /// Try to store the mark in kmp (key map or kmer position) returns true while
    /// we need to extend to keep trying. The actual mark that is extended may change
    /// in the process.
    fn try_store_mark(&mut self, mark: &mut XmerLoc, kmp: &mut Vec<ExtPosEtc>) -> Result<bool> {
        let old_stored_p = kmp[mark.idx];

        if old_stored_p.is_zero() {
            kmp[mark.idx].set(mark.p);
            return Ok(false);
        }
        if old_stored_p.is_replaceable_by(mark.p) {
            if old_stored_p.pos() == mark.p.pos() {
                // set and already mark.p. Leave the bit states.
                return Ok(false);
            }
            kmp[mark.idx].set(mark.p);
            if old_stored_p.x() == mark.p.x() {
                // same extension means same base k-mer origin. this is a duplicate.
                kmp[mark.idx].mark_more_recurs_upseq();
            }
            dbg_print!("{} -> ?", mark);
            mark.p = old_stored_p; // to be extended next
        } else if old_stored_p.extension() == mark.p.extension() {
            // Note: same extension and hash means same k-mer origin: identical k-mer sequence.
            kmp[mark.idx].mark_more_recurs_upseq();
        }
        // collision between hashes, the one in mark.p will be extended and tried again.
        Ok(true)
    }
}

impl XmerHash for XmerHasher {
    fn get_kmerlen(&self) -> usize {
        self.kmerlen
    }
    fn get_overbit(&self) -> usize {
        self.overbit
    }
}
