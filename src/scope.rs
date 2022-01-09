use crate::kmer::Kmer;
use crate::kmer::TwoBit;
use crate::kmerconst::KmerConst;
use crate::kmerloc::{ExtPosEtc, KmerLoc};
use crate::kmerstore::KmerStore;
use crate::new_types::position::Position;
use crate::past_scope::PastScope;
use crate::rdbg::STAT_DB;
use anyhow::Result;
use std::cmp;

pub trait Scope {
    fn get_mark(&self) -> Option<&KmerLoc>;
    fn get_kc(&self) -> &KmerConst;
    fn get_p(&self) -> u64;
    fn get_i(&self) -> usize;
    fn get_d(&self, i: usize) -> &Kmer<u64>;
    fn is_repetitive(&self) -> bool;
    fn clear_p_extension(&mut self);
    fn increment(&mut self, b2: TwoBit) -> bool;
    fn extend_p(&mut self);
    fn mark_reset(&mut self);
    fn set_mark(&mut self, idx: usize, p: ExtPosEtc, x: usize);
    fn set_period(&mut self, period: Position);
    fn increment_for_extension(&mut self, ks: &KmerStore) -> Result<()>;
    fn dist_if_repetitive(
        &self,
        stored_p: ExtPosEtc,
        mark_p: ExtPosEtc,
        max_dist: Position,
    ) -> Option<Position>;
    fn unset_period(&mut self);

    /// is occ complete? call na complete_kmer() - self.i increment.
    fn all_kmers(&self) -> bool {
        self.get_i() >= self.get_kc().venster
    }

    fn is_mark_out_of_scope(&self, mark: &KmerLoc) -> bool {
        let p = self.get_p();
        dbg_assert!(p.pos() >= mark.p.pos(), "{:#x}, {:#x}", p, mark.p);
        p.pos() >= mark.p.pos() + (self.get_kc().no_xmers(p.x()) as u64).basepos_to_pos()
    }

    /// Manage mark, do we have any minimum?
    fn remark(&mut self, reset_extension_if_leaving: bool) -> Result<bool> {
        if self.all_kmers() {
            if let Some(mark) = self.get_mark() {
                if self.is_mark_out_of_scope(mark) {
                    if reset_extension_if_leaving {
                        self.clear_p_extension();
                    }
                    return self.set_next_mark();
                }
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn can_extend(&self) -> bool {
        self.get_p().x() + 1 < self.get_kc().extent.len()
    }

    fn handle_mark(&mut self, ks: &mut KmerStore) -> Result<()> {
        if let Some(mark) = self.get_mark() {
            let (min_idx, min_p) = mark.get();
            if self.is_repetitive() && ks.kmp[min_idx].is_set() {
                ks.kmp[min_idx].set_repetitive();
                return Ok(());
            }
            let stored_p = ks.kmp[min_idx];

            if ks.kmp[min_idx].is_zero() {
                ks.set_kmp(min_idx, min_p);
                return Ok(());
            } else if stored_p.is_replaceable_by(min_p) {
                if dbgx!(stored_p.is_set_and_not(min_p)) {
                    ks.pos_max = cmp::max(ks.pos_max, self.get_p().pos());
                    let mut new_scp = PastScope::new(ks, self.get_kc(), &stored_p, min_idx)?;

                    // unset when kmer is not observed before bound.1:
                    if new_scp.get_p().is_set() {
                        dbg_print!(
                            "resolving past for [{:x}], {:#x} <= {:#x}",
                            min_idx,
                            stored_p,
                            min_p
                        );
                        ks.set_kmp(min_idx, min_p);
                        new_scp.handle_mark(ks)?;
                    }
                }
                // .. else set and already min_p. Then leave dupbit state.
                return Ok(());
            }
            if stored_p.extension() == min_p.extension() {
                // If a kmer occurs multiple times within an extending readlength (repetition),
                // only the first gets a position. During mapping this should be kept in mind.
                if let Some(dist) = self.dist_if_repetitive(stored_p, mark.p, ks.rep_max_dist) {
                    self.set_period(dist);
                    ks.kmp[min_idx].set_repetitive();
                    return Ok(());
                }
                ks.kmp[min_idx].set_dup();
            } else {
                dbg_print!("not replacable, extend..");
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
            if self.set_next_mark()? {
                self.handle_mark(ks)?;
                break;
            }
        }
        Ok(())
    }

    /// bepaal van alle kmers het nieuwe minimum/optimum (na leaving mark of extensie)
    /// bij repetitive DNA kan geen enkel optimum de uitkomst zijn.
    fn set_next_mark(&mut self) -> Result<bool> {
        self.mark_reset();
        let x = self.get_p().x();
        let kc = self.get_kc();
        let base = self.get_i() - kc.kmerlen;
        let bin = kc.get_kmers(x);
        // reverse is logischer omdat we bij gelijke extensie voor een lagere positie kiezen.
        for i in (0..kc.no_xmers(x)).rev() {
            let _ = self.set_if_optimum(x, base, (bin.0 + i, bin.1 + i));
        }
        Ok(self.get_mark().is_some())
    }

    fn is_xmer_complete(&self, x: usize) -> bool {
        let kc = self.get_kc();
        self.get_i() >= kc.kmerlen + kc.afstand(x)
    }
    /// voor een offset i en extensie x, maak de kmer/hash en zet mark + return true als optimum.
    fn set_if_optimum(&mut self, x: usize, base: usize, bin: (usize, usize)) -> bool {
        if self.is_xmer_complete(x) {
            let kc = self.get_kc();
            let base = self.get_i() - kc.kmerlen;
            let nk = kc.no_kmers;
            let kmer1 = self.get_d(base.wrapping_sub(bin.0) % nk);
            let mut hash = kmer1.get_idx(bin.0 <= bin.1);
            let mut p = if let Some(mark) = self.get_mark() {
                match hash.cmp(&mark.get_idx()) {
                    cmp::Ordering::Less => {
                        self.get_p().pos_with_ext(x) - (bin.0 as u64).basepos_to_pos().as_u64()
                    }
                    cmp::Ordering::Greater => return false,
                    cmp::Ordering::Equal => {
                        let p =
                            self.get_p().pos_with_ext(x) - (bin.0 as u64).basepos_to_pos().as_u64();
                        if p >= mark.p {
                            return false;
                        }
                        p
                    }
                }
            } else {
                self.get_p().pos_with_ext(x) - (bin.0 as u64).basepos_to_pos().as_u64()
            };
            p ^= match bin.0.cmp(&bin.1) {
                cmp::Ordering::Less => {
                    hash ^= self.get_d(base.wrapping_sub(bin.1) % nk).get_idx(true);
                    p & 1
                }
                cmp::Ordering::Greater => {
                    hash ^= self.get_d(base.wrapping_sub(bin.1) % nk).get_idx(false);
                    !p & 1
                }
                cmp::Ordering::Equal => (p ^ kmer1.dna) & 1,
            };
            if self.is_repetitive() {
                p.set_repetitive();
            }
            self.set_mark(hash, p, x);
            // XXX maybe return self.all_kmers() ??
            true
        } else {
            false
        }
    }
}
