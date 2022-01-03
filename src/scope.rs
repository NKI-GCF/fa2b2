use crate::kmer::Kmer;
use crate::kmerconst::KmerConst;
use crate::kmerloc::{KmerLoc, PriExtPosOri};
use crate::kmerstore::KmerStore;
use crate::past_scope::PastScope;
use crate::rdbg::STAT_DB;
use anyhow::{ensure, Result};
use std::{cmp, fmt};

pub trait Scope {
    fn get_mark(&self) -> Option<&KmerLoc<u64>>;
    fn get_kc(&self) -> &KmerConst;
    fn get_p(&self) -> u64;
    fn get_i(&self) -> usize;
    fn get_d(&self, i: usize) -> &Kmer<u64>;
    fn is_repetitive(&self) -> bool;
    fn get_plim(&self) -> (u64, u64);
    fn clear_p_extension(&mut self);
    fn increment(&mut self, b2: u8) -> bool;
    fn extend_p(&mut self);
    fn mark_reset(&mut self);
    fn set_mark(&mut self, idx: usize, p: u64, x: usize);
    fn set_period(&mut self, period: u64);

    /// is occ complete? call na complete_kmer() - self.i increment.
    fn all_kmers(&self) -> bool {
        self.get_i() >= self.get_kc().venster
    }

    fn is_mark_out_of_scope(&self, mark: &KmerLoc<u64>) -> bool {
        let p = self.get_p();
        dbg_assert!(p.pos() >= mark.p.pos(), "{:#x}, {:#x}", p, mark.p);
        p.pos() >= mark.p.pos() + (self.get_kc().no_xmers(p.x()) << 1) as u64
    }

    /// sufficient kmers to have mimimum and do we have minimum?
    fn remark<T: PriExtPosOri>(
        &mut self,
        ks: &KmerStore<T>,
        reset_extension_if_leaving: bool,
    ) -> Result<bool> {
        if self.all_kmers() {
            if let Some(mark) = self.get_mark() {
                if self.is_mark_out_of_scope(mark) {
                    if reset_extension_if_leaving {
                        self.clear_p_extension();
                    }
                    return self.set_next_mark(ks);
                }
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn get_if_repetitive(&self, stored_pos: u64, mark_pos: u64) -> Option<u64> {
        assert!(mark_pos > stored_pos);
        let dist = mark_pos - stored_pos;
        let plim = self.get_plim();
        if dbgx!(stored_pos >= plim.0.pos() && stored_pos < plim.1.pos())
            && dbgx!(dist < self.get_kc().repetition_max_dist)
        {
            Some(dist)
        } else {
            None
        }
    }
    fn can_extend(&self) -> bool {
        self.get_p().x() + 1 < self.get_kc().extent.len()
    }

    fn is_before_end_of_contig(&self) -> bool {
        self.get_p().pos() < self.get_plim().1
    }

    fn increment_for_extension<T>(&mut self, ks: &KmerStore<T>) -> Result<()>
    where
        T: PriExtPosOri + fmt::LowerHex + Copy,
    {
        let b2 = ks.b2_for_p(self.get_p(), Some("(extension)"))?;
        ensure!(self.is_before_end_of_contig(), "running into end of contig");
        dbg_assert!(self.increment(b2));
        Ok(())
    }

    fn handle_mark<T>(&mut self, ks: &mut KmerStore<T>) -> Result<()>
    where
        T: PriExtPosOri + fmt::LowerHex + Copy,
    {
        if let Some(mark) = self.get_mark() {
            let (min_idx, min_p) = mark.get();
            let stored_p = ks.kmp[min_idx];

            if stored_p.is_replaceable_by(min_p) {
                if dbgx!(stored_p.is_set_and_not(min_p)) {
                    let mut new_scp = PastScope::new(ks, self.get_kc(), &stored_p, min_idx)?;

                    // unset when kmer is not observed before bound.1:
                    if new_scp.get_p().is_set() {
                        dbg_print!(
                            "resolving past for [{:x}], {:#x} <= {:#x}",
                            min_idx,
                            stored_p,
                            min_p
                        );
                        ks.kmp[min_idx].set(min_p);
                        new_scp.handle_mark(ks)?;
                    }
                } else if ks.kmp[min_idx].is_no_pos() && !self.is_repetitive() {
                    ks.kmp[min_idx].set(min_p);
                }
                // .. else set and already min_p. Then leave dupbit state.
                return Ok(());
            }
            if stored_p.extension() == min_p.extension() {
                // If a kmer occurs multiple times within an extending readlength (repetition),
                // only the first gets a position. During mapping this should be kept in mind.
                let mark_pos = mark.p.pos();
                if let Some(dist) = self.get_if_repetitive(stored_p.pos(), mark_pos) {
                    self.set_period(dist);
                    return Ok(());
                }
                if !self.is_repetitive() {
                    ks.kmp[min_idx].set_dup();
                }
            } else {
                dbg_print!("\t\t<!>");
            }
        }
        while self.can_extend() {
            self.extend_p();
            if self
                .get_mark()
                .filter(|m| self.is_mark_out_of_scope(m))
                .is_some()
            {
                if let Err(e) = self.increment_for_extension(ks) {
                    dbg_print!("{}", e);
                    self.clear_p_extension();
                    break;
                }
            }
            // every time we extend, we reiterate all for mark. May fail in repetitive DNA.
            if self.set_next_mark(ks)? {
                break;
            }
        }
        Ok(())
    }

    /// bepaal van alle kmers het nieuwe minimum/optimum (na leaving mark of extensie)
    /// bij repetitive DNA kan geen enkel optimum de uitkomst zijn.
    fn set_next_mark<T: PriExtPosOri>(&mut self, ks: &KmerStore<T>) -> Result<bool> {
        self.mark_reset();
        let x = self.get_p().x();
        let kc = self.get_kc();
        let bin = kc.get_kmers(x);
        // reverse is logischer omdat we bij gelijke extensie voor een lagere positie kiezen.
        for i in (0..kc.no_xmers(x)).rev() {
            let _ = self.set_if_optimum(i, x, bin, Some(ks));
        }
        Ok(self.get_mark().is_some())
    }

    fn is_xmer_complete(&self, x: usize) -> bool {
        let kc = self.get_kc();
        self.get_i() >= kc.kmerlen + kc.afstand(x)
    }
    /// voor een offset i en extensie x, maak de kmer/hash en zet mark + return true als optimum.
    fn set_if_optimum<T: PriExtPosOri>(
        &mut self,
        i: usize,
        x: usize,
        bin: (usize, usize),
        oks: Option<&KmerStore<T>>,
    ) -> bool {
        // XXX function is hot
        if self.is_xmer_complete(x) {
            let kc = self.get_kc();
            let base = self.get_i() - kc.kmerlen;
            let d = i + bin.0;
            let nk = kc.no_kmers;
            let kmer1 = self.get_d(base.wrapping_sub(d) % nk);
            let mut hash = kmer1.get_idx(bin.0 <= bin.1);
            let mut p = if let Some(mark) = self.get_mark() {
                match hash.cmp(&mark.get_idx()) {
                    cmp::Ordering::Less => self.get_p().with_ext(x) - (d << 1) as u64,
                    cmp::Ordering::Greater => return false,
                    cmp::Ordering::Equal => {
                        let p = self.get_p().with_ext(x) - (d << 1) as u64;
                        if p >= mark.p {
                            return false;
                        }
                        p
                    }
                }
            } else {
                self.get_p().with_ext(x) - (d << 1) as u64
            };
            p ^= match bin.0.cmp(&bin.1) {
                cmp::Ordering::Less => {
                    hash ^= self.get_d(base.wrapping_sub(i + bin.1) % nk).get_idx(true);
                    p & 1
                }
                cmp::Ordering::Greater => {
                    hash ^= self.get_d(base.wrapping_sub(i + bin.1) % nk).get_idx(false);
                    !p & 1
                }
                cmp::Ordering::Equal => (p ^ kmer1.dna) & 1,
            };
            if let Some(rep) = oks.and_then(|ks| ks.repeat.get(&hash)) {
                let ks = oks.unwrap();
                let stored_p = &ks.kmp[hash];
                if stored_p.is_set()
                    && p.pos() != stored_p.pos()
                    && dbgf!(
                        p.pos() < stored_p.pos() + rep.1 as u64,
                        "{} {:#x} Rep:{:?}",
                        hash,
                        rep
                    )
                {
                    return false;
                }
            }
            self.set_mark(hash, p, x);
            // XXX maybe return self.all_kmers() ??
            true
        } else {
            false
        }
    }
}
