use crate::kmer::Kmer;
use crate::kmerconst::KmerConst;
use crate::kmerloc::{KmerLoc, PriExtPosOri};
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::Result;
use std::fmt;

pub struct HeadScope<'a> {
    kc: &'a KmerConst,
    pub p: u64,
    d: Vec<Kmer<u64>>, // misschien is deze on the fly uit ks te bepalen?
    pub mark: KmerLoc<u64>,
    pub i: usize,
    pub mod_i: usize,
    pub plim: (u64, u64),
    pub period: u64,
}

impl<'a> HeadScope<'a> {
    pub fn new(plim: (u64, u64), kc: &'a KmerConst, p: u64) -> Self {
        HeadScope {
            kc,
            p,
            d: vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers],
            mark: KmerLoc::new(usize::max_value(), p.extension()),
            i: 0,
            mod_i: 0,
            plim,
            period: 0,
        }
    }

    // hier krijgen we nieuwe sequence, zijn geen past scope aan het behandelen, of zo.
    // .i & .p increments en kmer .d[] updates vinden plaats.
    pub fn complete_and_update_mark<T>(&mut self, b2: u8, ks: &mut KmerStore<T>) -> Result<()>
    where
        T: PriExtPosOri + fmt::LowerHex + Copy,
    {
        // XXX: function is very hot
        if self.increment(b2) {
            // mark.p hoeft niet gezet te zijn.
            for x in 0..=self.p.x() {
                //
                if self.set_if_optimum(0, x, Some(ks)) {
                    if self.all_kmers() && !self.handle_mark(ks)? {
                        dbg_print!("Unable to find mark");
                    }
                    return Ok(());
                }
            }
            if self.mark_set_considering_leaving(ks, true)? && !self.handle_mark(ks)? {
                dbg_print!("Unable to find mark");
            }
        }
        Ok(())
    }
}

impl<'a> Scope for HeadScope<'a> {
    fn get_mark(&self) -> Option<&KmerLoc<u64>> {
        if self.mark.is_set() {
            Some(&self.mark)
        } else {
            None
        }
    }
    fn get_kc(&self) -> &KmerConst {
        self.kc
    }
    fn get_p(&self) -> u64 {
        self.p
    }
    fn get_plim(&self) -> (u64, u64) {
        self.plim
    }
    fn get_i(&self) -> usize {
        self.i
    }
    fn get_d(&self, i: usize) -> &Kmer<u64> {
        &self.d[i]
    }
    fn is_repetitive(&self) -> bool {
        self.period != 0
    }
    fn set_period(&mut self, period: u64) {
        self.period = period;
    }
    fn clear_p_extension(&mut self) {
        self.p.clear_extension();
    }

    /// add twobit to k-mers, update k-mer vec, increment pos and update orientation
    /// true if we have at least one kmer.
    fn increment(&mut self, b2: u8) -> bool {
        // XXX: function is hot
        if self.i >= self.kc.kmerlen {
            let old_d = self.d[self.mod_i];
            self.mod_i += 1;
            if self.mod_i == self.kc.no_kmers {
                self.mod_i = 0;
            }
            self.d[self.mod_i] = old_d;
        }
        dbg_assert!(self.p.is_set());
        // first bit is strand bit, set according to kmer orientation bit.
        self.p &= !1;
        self.p += 2 + self.d[self.mod_i].update(b2);
        self.i += 1;
        self.i >= self.kc.kmerlen
    }
    fn extend_p(&mut self) {
        self.p.extend();
    }
    fn mark_reset(&mut self) {
        self.mark.reset();
    }
    fn set_mark(&mut self, idx: usize, p: u64, x: usize) {
        self.mark.set(idx, p, x);
    }
}

impl<'a> fmt::Display for HeadScope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.pos() as usize;
        let mp = self.mark.p.pos() as usize;
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat(((p >> 1) - self.kc.venster) * 5);
        let r = (p - mp) >> 1;
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
