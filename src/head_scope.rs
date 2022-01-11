use crate::kmer::{Kmer, TwoBit};
use crate::kmerconst::KmerConst;
use crate::kmerloc::{ExtPosEtc, KmerLoc};
use crate::kmerstore::KmerStore;
use crate::new_types::position::Position;
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::Result;
use std::fmt;

pub struct HeadScope<'a> {
    kc: &'a KmerConst,
    pub p: ExtPosEtc,
    d: Vec<Kmer>,
    z: Vec<usize>,
    pub mark: KmerLoc,
    pub i: usize,
    pub mod_i: usize,
    pub period: Position,
}

impl<'a> HeadScope<'a> {
    pub fn new(kc: &'a KmerConst) -> Self {
        HeadScope {
            kc,
            p: ExtPosEtc::zero(),
            d: vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers],
            z: (0..kc.no_kmers).into_iter().collect(),
            mark: KmerLoc::new(usize::max_value(), ExtPosEtc::zero()),
            i: 0,
            mod_i: 0,
            period: Position::zero(),
        }
    }

    // hier krijgen we nieuwe sequence, zijn geen past scope aan het behandelen, of zo.
    // .i & .p increments en kmer .d[] updates vinden plaats.
    pub fn complete_and_update_mark(&mut self, b2: TwoBit, ks: &mut KmerStore) -> Result<()> {
        if self.increment(b2) {
            let base = self.get_i() - self.kc.kmerlen;
            for x in 0..=self.p.x() {
                if self.set_if_optimum(x, base, self.kc.get_kmers(x)) {
                    break;
                }
            }
            if self.remark(true)? {
                self.handle_mark(ks)?;
            }
        }
        Ok(())
    }
    fn is_on_last_contig(&self, ks: &KmerStore, pos: Position) -> bool {
        pos >= ks.contig.last().unwrap().twobit
    }
    fn pick_mark(&mut self, x: usize) -> (usize, ExtPosEtc) {
        let med = self.kc.no_kmers >> 1;
        let i = self
            .z
            .select_nth_unstable_by(med, |&a, &b| self.d[a].cmp(&self.d[b]))
            .1;
        self.d[*i].get_hash_and_p(x)
    }
    // TODO: binary search to highest set for baseidx. with dupbit 0.
}

impl<'a> Scope for HeadScope<'a> {
    fn get_mark(&self) -> Option<(usize, ExtPosEtc)> {
        self.mark.get()
    }
    fn get_kc(&self) -> &KmerConst {
        self.kc
    }
    fn get_p(&self) -> ExtPosEtc {
        self.p
    }
    fn get_i(&self) -> usize {
        self.i
    }
    fn get_d(&self, i: usize) -> &Kmer {
        &self.d[i]
    }
    fn is_repetitive(&self) -> bool {
        self.period.is_set()
    }
    fn set_period(&mut self, period: Position) {
        dbg_assert!(self.p.is_zero() || period < self.p.pos());
        self.period = period;
    }
    fn unset_period(&mut self) {
        self.set_period(Position::zero());
    }
    fn clear_p_extension(&mut self) {
        self.p.clear_extension();
    }
    fn increment_for_extension(&mut self, ks: &KmerStore) -> Result<()> {
        let b2 = ks.b2_for_p(self.get_p().pos(), false)?;
        dbg_assert!(self.increment(b2));
        Ok(())
    }

    fn dist_if_repetitive(&self, ks: &KmerStore, stored_p: ExtPosEtc) -> Option<Position> {
        // FIXME: the contig check was removed. for the upper bound that makes sense for head
        // scope, as upperbound is pos_max rather than previous
        let stored_pos = stored_p.pos();
        let mark_pos = self.mark.p.pos();
        dbg_assert!(mark_pos > stored_pos);
        if self.is_on_last_contig(ks, stored_pos) {
            let dist = mark_pos - stored_pos;
            if dist < ks.rep_max_dist {
                return Some(dist);
            }
        }
        None
    }

    /// add twobit to k-mers, update k-mer vec, increment pos and update orientation
    /// true if we have at least one kmer.
    fn increment(&mut self, b2: TwoBit) -> bool {
        // XXX: function is hot
        if self.i >= self.kc.kmerlen {
            let old_d = self.d[self.mod_i];
            self.mod_i += 1;
            if self.mod_i == self.kc.no_kmers {
                self.mod_i = 0;
            }
            self.d[self.mod_i] = old_d;
        }
        // first bit is strand bit, set according to kmer orientation bit.
        self.p.set_ori(self.d[self.mod_i].update(b2));
        self.p.incr_pos();
        self.i += 1;
        self.i >= self.kc.kmerlen
    }
    fn extend_p(&mut self) {
        self.p.extend();
    }
    fn mark_reset(&mut self) {
        self.mark.reset();
    }
    fn set_mark(&mut self, idx: usize, p: ExtPosEtc, x: usize) {
        format!("[{:x}] = {:?} | x({})", idx, p, x);
        self.mark.set(idx, p, x);
    }
}

impl<'a> fmt::Display for HeadScope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.unshift_pos() as usize;
        let mp = self.mark.p.unshift_pos() as usize;
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat((p - self.kc.venster) * 5);
        let r = p - mp;
        if r == 0 {
            let x = self.kc.venster - n;
            let s = if x != 0 {
                " ".repeat(x << 2) + "|"
            } else {
                String::from("")
            };
            write!(f, "{2}<{3}{: ^1$x}>", self.mark.get_idx(), n << 2, o, s)
        } else if r + self.kc.kmerlen == self.kc.venster {
            let x = self.kc.venster - n;
            let s = if x != 0 {
                String::from("|") + &" ".repeat(x << 2)
            } else {
                String::from("")
            };
            write!(f, "{2}<{: ^1$x}{3}>", self.mark.get_idx(), n << 2, o, s)
        } else {
            //let l = self.kc.venster - r - n;
            //let ls = if o {" ".repeat(o) + "|"} else {String::from("")};
            //let rs = if l {String::from("|") + &" ".repeat(l << 2)} else {String::from("")};
            //write!(f, "{2}<{3}|{: ^1$x}|{4}>", self.mark.idx, n << 2, o, ls, rs)
            write!(f, "")
        }
    }
}
