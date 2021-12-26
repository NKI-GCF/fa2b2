use crate::kmer::Kmer;
use crate::kmerconst::KmerConst;
use crate::kmerloc::{KmerLoc, PriExtPosOri};
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use anyhow::{ensure, Result};
use std::{cmp, fmt};

pub struct Scope<'a> {
    pub kc: &'a KmerConst,
    pub p: u64,
    d: Vec<Kmer<u64>>, // misschien is deze on the fly uit ks te bepalen?
    pub mark: KmerLoc<u64>,
    pub i: usize,
    pub mod_i: usize,
    pub plim: (u64, u64),
}

impl<'a> Scope<'a> {
    pub fn new(plim: (u64, u64), kc: &'a KmerConst, p: u64) -> Self {
        Scope {
            kc,
            p,
            d: vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers],
            mark: KmerLoc::new(usize::max_value(), p.extension()),
            i: 0,
            mod_i: 0,
            plim,
        }
    }

    pub fn rebuild<T: PriExtPosOri>(
        &mut self,
        ks: &KmerStore<T>,
        plim: (u64, u64),
        p: u64,
    ) -> Result<bool> {
        /* Er gaat iets mis, blijkt uit repetitie telling die afwijkt afhankelijk
        van de scope heropbouw. Indien met plim.0, geeft dit andere repetition
        aantallen dan wanneer de eerste 2bit in scope wordt gebruikt. Voor de
        meeste contigs meer repetitie, enkele minder. Duizendtallen op de primary
        contigs.
        De positie waarvoor herbouwt wordt zou een repetition of duplicaat kunnen
        zijn, maar repetition wordt alleen geteld bij nieuwe sequentie. Dus afh
        van scope rebuilding moeten er verschillende xmers 'bezet' zijn, wat dan
        invloed heeft op de telling bij nieuwe sequence, lijkt mij de logischte
        verklaring. */

        // FIXME: dit zou hetzelfde resultaat moeten geven!
        // TODO: verbeter rebuild code en test impact hierop.
        // Minder of geen verschil bij deze twee alternatieven is beter:
        // self.p = p.extension() | plim.0;
        self.p = self.kc.leftmost_of_scope(p, plim.0);

        self.plim = plim;
        self.i = 0;
        self.mod_i = 0;
        self.mark.reset();

        let x = p.x();
        loop {
            let b2 = ks.b2_for_p(self.p).unwrap();
            dbg_print!("=> b2 {:x}, p: {:#x} <=", b2, self.p);
            if self.complete_and_update_mark(b2, x)? {
                if self.mark.p == p {
                    return Ok(true);
                }
                if self.p.pos() >= self.plim.1 {
                    self.p.clear();
                    return Ok(false);
                }
            }
        }
    }

    fn is_xmer_complete(&self, x: usize) -> bool {
        /*dbgf!(
            self.i >= self.kc.kmerlen + self.kc.afstand(x),
            "{}: {} >= {} + {}?",
            self.i,
            self.kc.kmerlen,
            self.kc.afstand(x)
        )*/
        self.i >= self.kc.kmerlen + self.kc.afstand(x)
    }

    pub fn increment(&mut self, b2: u8) {
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
        self.p += 2;
        self.p ^= 1 & (self.p ^ self.d[self.mod_i].update(b2));
        self.i += 1;
    }

    /// add b2 to kmer, move .p according to occ orientation
    /// return whether required nr of kmers for struct scope were seen.
    fn update_mark(&mut self, x_start: usize) -> bool {
        // XXX: function is very hot
        if self.i >= self.kc.kmerlen {
            for x in x_start..=self.p.x() {
                if self.set_if_optimum(0, x) {
                    return dbgf!(self.all_kmers(), "[update_mark()]: {}");
                }
            }
            self.mark.is_set() && self.all_kmers()
        } else {
            false
        }
    }

    // hierin vinden .i & .p increments and kmer .d[] update plaats.
    pub fn complete_and_update_mark(&mut self, b2: u8, x_start: usize) -> Result<bool> {
        // XXX: function is very hot
        self.increment(b2);
        if self.i >= self.kc.kmerlen {
            for x in x_start..=self.p.x() {
                if self.set_if_optimum(0, x) {
                    return Ok(dbgf!(self.all_kmers(), "[complete_and_update_mark()]: {}"));
                }
            }
            if self.mark_is_first_or_leaving() {
                self.set_next_mark()?;
                return Ok(true);
            }
        }
        return Ok(false);
    }

    /// extend positie (als kmer/hash niet replaceble was); true indien mogelijk.
    pub fn extend(&mut self) -> bool {
        let ret = self.p.x() + 1 < self.kc.ext_max;
        if ret {
            self.p.extend();
        }
        ret
    }

    /// is occ complete? call na complete_kmer() - self.i increment.
    fn all_kmers(&self) -> bool {
        self.i >= self.kc.readlen
    }

    /// is de minimum/optimum leaving? eerste is speciaal.
    fn mark_is_first_or_leaving(&self) -> bool {
        match self.i.cmp(&self.kc.readlen) {
            cmp::Ordering::Greater => {
                let p_pos = self.p.pos();
                let dist = (self.kc.no_kmers << 1) as u64;
                p_pos >= self.mark.p.pos() + dist
            }
            cmp::Ordering::Less => false,
            cmp::Ordering::Equal => true,
        }
    }

    /// bepaal van alle kmers het nieuwe minimum/optimum (na leaving mark of extensie)
    pub fn set_next_mark(&mut self) -> Result<()> {
        self.mark.reset();
        let x = self.p.x();
        let afs = self.kc.afstand(x);
        ensure!(self.kc.no_kmers > afs, "Couldn't obtain new mark.");
        for i in (0..(self.kc.no_kmers - afs)).rev() {
            let _ = self.set_if_optimum(i, x);
        }
        dbg_assert!(self.mark.is_set());
        Ok(())
    }

    /// continue rebuild //FIXME: rewrite, if this function is necessary.
    pub fn extend_kmer_stack<T: PriExtPosOri>(&mut self, ks: &KmerStore<T>) -> Result<()> {
        let x = self.p.x();

        // set mark after extension, if possible
        self.mark.reset();
        if self.is_xmer_complete(x) {
            let _ = self.set_if_optimum(0, x);
        }

        while self.p.pos() < self.plim.1 {
            let b2 = ks.b2_for_p(self.p).unwrap();
            self.increment(b2);
            if self.update_mark(x) && self.mark_is_first_or_leaving() {
                self.set_next_mark()?;
                if ks.kmp[self.mark.idx].is_same(self.mark.p) {
                    // retracked => finished
                    break;
                }
            }
        }
        Ok(())
    }

    /// voor een offset i en extensie x, maak de kmer/hash en zet mark + return true als optimum.
    fn set_if_optimum(&mut self, i: usize, x: usize) -> bool {
        // XXX function is hot
        if self.i >= self.kc.kmerlen + self.kc.afstand(x) {
            let e = self.kc.extension[x];
            let base = self.i - self.kc.kmerlen;
            let d = i + e.0 as usize;
            let nk = self.kc.no_kmers;
            let kmer1 = self.d[base.wrapping_sub(d) % nk];
            let mut hash = kmer1.get_idx(e.0 <= e.1);
            let mut p = match hash.cmp(&self.mark.idx) {
                cmp::Ordering::Less => self.p.with_ext(x) - (d << 1) as u64,
                cmp::Ordering::Greater => return false,
                cmp::Ordering::Equal => {
                    let p = self.p.with_ext(x) - (d << 1) as u64;
                    if p >= self.mark.p {
                        return false;
                    }
                    p
                }
            };
            p ^= match e.0.cmp(&e.1) {
                cmp::Ordering::Less => {
                    hash ^= self.d[base.wrapping_sub(i + e.1 as usize) % nk].get_idx(true);
                    p & 1
                }
                cmp::Ordering::Greater => {
                    hash ^= self.d[base.wrapping_sub(i + e.1 as usize) % nk].get_idx(false);
                    !p & 1
                }
                cmp::Ordering::Equal => (p ^ kmer1.dna) & 1,
            };
            dbgf!(
                self.mark.set(hash, p, x),
                "{:?} [{:x}] = {:x} | x({})",
                hash,
                p,
                x
            );
            true
        } else {
            false
        }
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

    impl<'a> Scope<'a> {
        /// de current kmer.
        fn kmer(&self) -> Kmer<u64> {
            self.d[self.mod_i]
        }
    }

    const READLEN: usize = 16;
    const SEQLEN: usize = 250;

    #[test]
    fn test_push_b2() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut occ = Scope::new((0, 100), &kc, 0);
        for _ in 0..6 {
            occ.increment(1);
            occ.update_mark(0);
        }
        let mut kmer = occ.kmer();
        dbg_assert_eq!(kmer.dna, 0x55);
        dbg_assert_eq!(kmer.rc, 0xff);
        dbg_assert_eq!(occ.p, 0xd);

        occ.increment(2);
        occ.update_mark(0);
        kmer = occ.kmer();
        dbg_assert_eq!(kmer.dna, 0x95);
        dbg_assert_eq!(kmer.rc, 0xfc);
        dbg_assert_eq!(occ.p, 0xf);
        dbg_assert_eq!(kmer.get_idx(true), 0x6a);
    }
    #[test]
    fn occurrence() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut occ = Scope::new((0, 100), &kc, 0);
        for _ in 0..8 {
            occ.increment(0);
            occ.update_mark(0);
        }
        for _ in 0..8 {
            occ.increment(1);
            occ.update_mark(0);
        }
        for _ in 0..8 {
            occ.increment(2);
            occ.update_mark(0);
        }
        for _ in 0..8 {
            occ.increment(3);
            occ.update_mark(0);
        }
    }
}
