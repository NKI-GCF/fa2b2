use std::collections::VecDeque;

use crate::kmer::Kmer;
pub use crate::kmerconst::KmerConst;
use crate::kmerloc::{KmerLoc, PriExtPosOri};
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use anyhow::{ensure, Result};
use std::fmt;

pub struct Scope<'a> {
    pub kc: &'a KmerConst,
    pub p: u64,
    d: VecDeque<Kmer<u64>>, // misschien is deze on the fly uit ks te bepalen?
    pub mark: KmerLoc<u64>,
    pub i: usize,
    pub plim: (u64, u64),
}

impl<'a> Scope<'a> {
    pub fn new(plim: (u64, u64), kc: &'a KmerConst, p: u64) -> Self {
        Scope {
            kc,
            p,
            d: VecDeque::from(vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers]),
            mark: KmerLoc::new(usize::max_value(), p.extension()),
            i: 0,
            plim,
        }
    }

    pub fn set(&mut self, plim: (u64, u64), ext: u64) {
        self.p = ext | plim.0;
        self.plim = plim;
        self.i = 0;
        self.mark.idx = usize::max_value();
        self.mark.p = ext;
    }

    fn is_xmer_complete(&self, x: usize) -> bool {
        /*dbgf!(
            self.i < self.kc.kmerlen + self.kc.afstand(x),
            "{}: {} < {} + {}?",
            self.i,
            self.kc.kmerlen,
            self.kc.afstand(x)
        )*/
        self.i < self.kc.kmerlen + self.kc.afstand(x)
    }

    pub fn increment(&mut self, b2: u8) {
        let mut kmx = 0;
        if let Some(i) = self.i.checked_sub(self.kc.kmerlen) {
            let old_idx = i % self.kc.no_kmers;
            if old_idx + 1 != self.kc.no_kmers {
                kmx = old_idx + 1;
            }
            self.d[kmx] = self.d[old_idx];
        }
        // first bit is strand bit, set according to kmer deviant bit.
        self.p += if self.d[kmx].update(b2) { 3 } else { 2 } - (1 & self.p);
        self.i += 1;
    }

    /// add b2 to kmer, move .p according to occ orientation
    /// return whether required nr of kmers for struct scope were seen.
    // hierin vinden .i & .p increments and kmer .d[] update plaats. geen kmer ? false.
    fn complete(&mut self, b2: u8, x_start: usize) -> bool {
        self.increment(b2);
        if self.i >= self.kc.kmerlen {
            for x in x_start..=self.p.x() {
                if self.is_xmer_complete(x) || self.set_if_optimum(0, x) {
                    break;
                }
            }
            self.mark.is_set() && self.all_kmers()
        } else {
            false
        }
    }

    pub fn complete_and_update_mark(&mut self, b2: u8, x_start: usize) -> Result<bool> {
        if self.complete(b2, x_start) {
            if !self.mark_remains() {
                self.set_next_mark()?;
            }
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// geef de current kmer.
    pub fn kmer(&self) -> Kmer<u64> {
        // Todo: on the fly uit ks.
        self.d[(self.i - self.kc.kmerlen) % self.kc.no_kmers]
    }

    /// extend positie (als kmer/hash niet replaceble was); true indien mogelijk.
    // waar wordt dit ongedaan gemaakt??
    pub fn extend(&mut self) -> bool {
        let ret = self.p.x() + 1 < self.kc.ext_max;
        if ret {
            self.p.extend();
        }
        ret
    }

    /// is occ complete? call na complete_kmer() - self.i increment.
    pub fn all_kmers(&self) -> bool {
        self.i >= self.kc.readlen
    }

    /// true als voor de occurance de hash een nieuw minimum/optimum is.
    fn hash_is_optimum(&mut self, hash: usize, p: u64) -> bool {
        let xh = hash ^ self.kc.ext_domain(p.x());
        /*dbg_print!(
            "{:x} ^ {:x} = {:x} < {:x}? || (== && {:x} < {:x})",
            hash,
            self.kc.ext_domain(p.x()),
            xh,
            self.mark.idx,
            p,
            self.mark.p
        );*/
        xh < self.mark.idx || (xh == self.mark.idx && p < self.mark.p)
    }

    /// is de minimum/optimum leaving?
    pub fn mark_is_leaving(&self) -> bool {
        //let afs = self.kc.afstand(self.mark.p.x());
        let p_pos = self.p.pos();
        let dist = (self.kc.no_kmers << 1) as u64;
        p_pos == self.mark.p.pos() + dist
    }

    pub fn mark_remains(&self) -> bool {
        !self.mark_is_leaving()
    }

    /// bepaal van alle kmers het nieuwe minimum/optimum (na leaving mark of extensie)
    pub fn set_next_mark(&mut self) -> Result<()> {
        self.mark.reset();
        let x = self.p.x();
        let afs = self.kc.afstand(x);
        ensure!(self.kc.no_kmers > afs, "Couldn't obtain new mark.");
        for i in 0..(self.kc.no_kmers - afs) {
            let _ = self.set_if_optimum(i, x);
        }
        dbg_assert!(self.mark.is_set());
        Ok(())
    }

    /// continue rebuild
    pub fn extend_kmer_stack<T: PriExtPosOri>(&mut self, ks: &KmerStore<T>) -> Result<()> {
        let x = self.p.x();

        // set mark after extension, if possible
        self.mark.reset();
        if !self.is_xmer_complete(x) {
            let _ = self.set_if_optimum(0, x);
        }

        while self.p.pos() < self.plim.1 {
            let b2 = ks.b2_for_p(self.p).unwrap();
            if self.complete(b2, x) && self.mark_is_leaving() {
                self.set_next_mark()?;
                if ks.kmp[self.mark.idx].is_same(self.mark.p) {
                    // retracked => finished
                    break;
                }
            }
        }
        Ok(())
    }
    fn get_xmer_and_strand(&self, i: usize, afs: usize) -> (usize, bool) {
        let base = self.i - self.kc.kmerlen;
        let d_i = base.wrapping_sub(i) % self.kc.no_kmers;
        let mut kmer = self.d[d_i];
        if afs > 0 {
            let d_i2 = base.wrapping_sub(afs + i) % self.kc.no_kmers;
            dbg_print!("dna:{:#x}, dna2:{:#x}", self.d[d_i2].dna, kmer.dna);
            kmer.hash(self.d[d_i2]);
        }
        (kmer.get_idx(), kmer.is_template())
    }

    /// voor een offset i en extensie x, maak de kmer/hash en zet mark + return true als optimum.
    fn set_if_optimum(&mut self, i: usize, x: usize) -> bool {
        let afs = self.kc.afstand(x);
        let (hash, is_template) = self.get_xmer_and_strand(i, afs);
        let mut p = self.p.with_ext(x) - ((afs + i) << 1) as u64;
        if is_template {
            p |= 1;
        }
        let ret = self.hash_is_optimum(hash, p);
        if ret {
            dbgf!(
                self.mark.set(hash, p, x),
                "{:?} [{:x}] = {:x} | x({})",
                hash,
                p,
                x
            );
        }
        ret
    }
}

impl<'a> fmt::Display for Scope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.pos() as usize;
        let mp = self.mark.p.pos() as usize;
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat(((p >> 1) - self.kc.readlen) * 5);
        let r = (p - mp) >> 1;
        if r == 0 {
            let x = self.kc.readlen - n;
            let s = if x != 0 {
                " ".repeat(x << 2) + "|"
            } else {
                String::from("")
            };
            write!(f, "{2}<{3}{: ^1$x}>", self.mark.idx, n << 2, o, s)
        } else if r + self.kc.kmerlen == self.kc.readlen {
            let x = self.kc.readlen - n;
            let s = if x != 0 {
                String::from("|") + &" ".repeat(x << 2)
            } else {
                String::from("")
            };
            write!(f, "{2}<{: ^1$x}{3}>", self.mark.idx, n << 2, o, s)
        } else {
            //let l = self.kc.readlen - r - n;
            //let ls = if o {" ".repeat(o) + "|"} else {String::from("")};
            //let rs = if l {String::from("|") + &" ".repeat(l << 2)} else {String::from("")};
            //write!(f, "{2}<{3}|{: ^1$x}|{4}>", self.mark.idx, n << 2, o, ls, rs)
            write!(f, "")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const READLEN: usize = 16;
    const SEQLEN: usize = 250;

    #[test]
    fn test_push_b2() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut occ = Scope::new((0, 100), &kc, 0);
        for _ in 0..6 {
            occ.complete(1, 0);
        }
        let mut kmer = occ.kmer();
        dbg_assert_eq!(kmer.dna, 0x55);
        dbg_assert_eq!(kmer.rc, 0xff);
        dbg_assert_eq!(occ.p, 0xd);

        occ.complete(2, 0);
        kmer = occ.kmer();
        dbg_assert_eq!(kmer.dna, 0x95);
        dbg_assert_eq!(kmer.rc, 0xfc);
        dbg_assert_eq!(occ.p, 0xf);
        dbg_assert_eq!(kmer.get_idx(), 0x6a);
    }
    #[test]
    fn occurrence() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut occ = Scope::new((0, 100), &kc, 0);
        for _ in 0..8 {
            occ.complete(0, 0);
        }
        for _ in 0..8 {
            occ.complete(1, 0);
        }
        for _ in 0..8 {
            occ.complete(2, 0);
        }
        for _ in 0..8 {
            occ.complete(3, 0);
        }
    }
}
