use crate::kmerconst::XmerHash;
use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::position::Position;
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use anyhow::{ensure, Result};
use crossbeam_channel::{Receiver, Sender, TryRecvError};
use std::collections::VecDeque;
use std::sync;
use std::sync::atomic::AtomicUsize;

pub(crate) struct XmerHasher {
    thread_nr: usize,
    no_threads: usize,
    kmerlen: usize,
    overbit: usize,
    shift: usize,
    rep_max_dist: Position,
    rx_from_main: Receiver<XmerLoc>,
    tx_inter_thread: Sender<XmerLoc>,
    rx_inter_thread: Receiver<XmerLoc>,
    tx_to_main: Sender<(Vec<ExtPosEtc>, Vec<XmerLoc>)>,
    to_next: VecDeque<XmerLoc>,
}

type XmerChannels = (usize, Receiver<XmerLoc>, Sender<XmerLoc>, Receiver<XmerLoc>);

impl XmerHasher {
    pub(crate) fn new(
        no_threads: usize,
        kmerlen: usize,
        rep_max_dist: Position,
        xmer_channels: XmerChannels,
        tx_to_main: Sender<(Vec<ExtPosEtc>, Vec<XmerLoc>)>,
    ) -> XmerHasher {
        let bitlen = kmerlen * 2;
        let shift = bitlen - 1 - usize::try_from(no_threads.trailing_zeros()).unwrap();
        dbg_print!("thread {}", xmer_channels.0);

        XmerHasher {
            thread_nr: xmer_channels.0,
            no_threads,
            kmerlen,
            overbit: 1 << (kmerlen * 2 + 1),
            shift,
            rep_max_dist,
            rx_from_main: xmer_channels.1,
            tx_inter_thread: xmer_channels.2,
            rx_inter_thread: xmer_channels.3,
            tx_to_main,
            to_next: VecDeque::<XmerLoc>::with_capacity(1 << 16),
        }
    }
    pub(crate) fn work(&mut self, shutdown_poll: sync::Arc<AtomicUsize>) -> Result<()> {
        let mut kmp = vec![ExtPosEtc::default(); 1 << self.shift]; // position + strand per xmer handled in this thread.
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
                Ok(mark) => {
                    dbg_assert!(
                        self.is_for_this_thread(&mark),
                        "Thread {} received {} ({:#x}) from main or inter thread, but it seems for thread {}",
                        self.thread_nr,
                        mark,
                        self.unhash_and_uncompress_to_kmer(mark.idx, mark.p.x()),
                        self.unhash_and_uncompress_to_kmer(mark.idx, mark.p.x()) >> self.shift
                    );
                    self.store_mark_and_extend(mark, &mut max_extended, &mut kmp);
                }
                Err(TryRecvError::Disconnected) => break,
                Err(TryRecvError::Empty) => {
                    if let Some(err) = self
                        .to_next
                        .pop_front()
                        .and_then(|m| self.tx_inter_thread.try_send(m).err())
                    {
                        self.to_next.push_front(err.into_inner());
                    } else if let Ok(mark) = self.rx_from_main.recv() {
                        dbg_assert!(
                            self.is_for_this_thread(&mark),
                            "Thread {} received {} ({:#x}) from main thread (after block), but it seems for thread {}",
                            self.thread_nr,
                            mark,
                            self.unhash_and_uncompress_to_kmer(mark.idx, mark.p.x()),
                            self.unhash_and_uncompress_to_kmer(mark.idx, mark.p.x()) >> self.shift
                        );
                        // We seem starved, block receive until more work.
                        // XXX if we have something else to do in the future this should change.
                        self.store_mark_and_extend(mark, &mut max_extended, &mut kmp);
                    }
                }
            }
        }
        dbg_print!("The main channel has shutdown.");
        let end_state = (1 << self.no_threads) - 1;
        loop {
            if let Ok(mark) = self.rx_inter_thread.try_recv() {
                // indicate we had something to process still by unsetting a bit (if it was set).
                let x_ptr = sync::Arc::as_ptr(&shutdown_poll);
                unsafe { &*x_ptr }
                    .fetch_and(!(1 << self.thread_nr), sync::atomic::Ordering::Relaxed);

                self.store_mark_and_extend(mark, &mut max_extended, &mut kmp);
            } else if self.to_next.is_empty() {
                //This may be undone upon new received..
                dbg_print!("Thread {} has nothing left to process.", self.thread_nr);
                let x_ptr = sync::Arc::as_ptr(&shutdown_poll);

                if unsafe { &*x_ptr }.fetch_or(1 << self.thread_nr, sync::atomic::Ordering::Relaxed)
                    == end_state
                {
                    break;
                }
            } else if let Some(err) = self
                .to_next
                .pop_front()
                .and_then(|m| self.tx_inter_thread.try_send(m).err())
            {
                self.to_next.push_front(err.into_inner());
            }
        }
        dbg_print!(
            "blocking send the packages back to the main thread in order of processing threads."
        );
        if self.thread_nr == 0 || self.rx_inter_thread.recv().is_ok() {
            dbg_print!("Now sending for thread {}", self.thread_nr);
            let x = self.tx_to_main.send((kmp, max_extended));
        }
        if self.thread_nr + 1 != self.no_threads {
            // message the next thread it can start sending.
            let x = self.tx_inter_thread.send(XmerLoc::default());
        }
        Ok(())
    }
    /// The actual work for the thread, hashing and extending index and updating its p, on collision.
    fn store_mark_and_extend(
        &mut self,
        mut mark: XmerLoc,
        max_extended: &mut Vec<XmerLoc>,
        kmp: &mut Vec<ExtPosEtc>,
    ) {
        while self.try_store_mark(&mut mark, kmp) {
            dbg_print!("mark {} needs extension", mark);
            // From the XmerHash trait, extension fails on last, when all ext bits are set.
            if self.extend_xmer(&mut mark).is_err() {
                max_extended.push(mark);
                break;
            }
            if !self.is_for_this_thread(&mark) {
                self.to_next.push_back(mark);
                break;
            }
        }
    }
    /// the top bits of the xmerhash indicate whether the extension can be stored in this kmp
    /// or the kmp of the next thread.
    pub(crate) fn is_for_this_thread(&self, mark: &XmerLoc) -> bool {
        (mark.idx >> self.shift) == self.thread_nr || (!mark.idx >> self.shift) == self.thread_nr
    }
    /// Try to store the mark in kmp (key map or kmer position) returns true while
    /// we need to extend to keep trying. The actual mark that is extended may change
    /// in the process.
    fn try_store_mark(&mut self, mark: &mut XmerLoc, kmp: &mut Vec<ExtPosEtc>) -> bool {
        let kmp_mask = kmp.len() - 1;
        let old_stored_p = kmp[mark.idx & kmp_mask];

        if old_stored_p.is_zero() {
            kmp[mark.idx & kmp_mask].set(mark.p);
            return false;
        }
        if old_stored_p.gt_ext_or_eq_and_le_pos(mark.p) {
            if old_stored_p.pos() == mark.p.pos() {
                // set and already mark.p. Leave the bit states.
                return false;
            }
            if old_stored_p.x() == mark.p.x() {
                // same extension means same base k-mer origin. this is a duplicate.
                if mark.p.pos() > self.rep_max_dist + old_stored_p.pos() {
                    mark.p.set_dup();
                } else {
                    kmp[mark.idx & kmp_mask].set_repetitive();
                    return false;
                }
            }
            kmp[mark.idx & kmp_mask].set(mark.p);
            dbg_print!("{} -> ?", mark);
            mark.p = old_stored_p; // to be extended next
        } else if old_stored_p.extension() == mark.p.extension() {
            // Note: same extension and hash means same k-mer origin: identical k-mer sequence.
            if mark.p.pos() > self.rep_max_dist + old_stored_p.pos() {
                kmp[mark.idx & kmp_mask].set_dup();
            } else {
                kmp[mark.idx & kmp_mask].set_repetitive();
                return false;
            }
        }
        // collision between hashes, the one in mark.p will be extended and tried again.
        true
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
