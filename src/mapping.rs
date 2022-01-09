use crate::kmer::{Kmer, ThreeBit, TwoBit};
use crate::kmerconst::KmerConst;
use crate::kmerloc::{ExtPosEtc, KmerLoc};
use crate::kmerstore::KmerStore;
use crate::new_types::position::Position;
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::Result;
use noodles_fastq as fastq;
use std::fmt;

pub struct Mapping<'a> {
    kc: &'a KmerConst,
    p: u64,
    i: usize,
    mod_i: usize,
    mark: KmerLoc,
    d: Vec<Kmer<u64>>, // misschien is deze on the fly uit ks te bepalen?
    z: Vec<usize>,
}

impl<'a> Mapping<'a> {
    pub fn new<T: ExtPosEtc + fmt::LowerHex>(
        ks: &KmerStore,
        kc: &'a KmerConst,
        record: fastq::Record,
    ) -> Result<Self> {
        let mut scp = Mapping {
            kc,
            p: 0,
            i: 0,
            mod_i: 0,
            mark: KmerLoc::new(usize::max_value(), 0),
            d: vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers],
            z: (0..kc.no_kmers).collect(),
        };
        let mut x = 0;
        let mut bin = kc.get_kmers(x);
        for b3 in record.sequence().iter().map(ThreeBit::from) {
            //let x = scp.p.x();
            //let bin = kc.get_kmers(x);
            if let Some(b2) = b3.as_twobit_if_not_n() {
                if scp.increment(b2) {
                    // we weten extension op voorhand.
                    let base = scp.i - kc.kmerlen;
                    if scp.set_if_optimum(x, base, bin) && scp.all_kmers() {}
                }
            }
        }
        Ok(scp)
    }
    fn pick_mark(&mut self, x: usize) -> (usize, Position) {
        let med = self.kc.no_kmers >> 1;
        let i = self
            .z
            .select_nth_unstable_by(med, |&a, &b| self.d[a].cmp(&self.d[b]))
            .1;
        self.d[*i].get_hash(x)
    }
}

impl<'a> Scope for Mapping<'a> {
    fn is_repetitive(&self) -> bool {
        panic!();
    }
    fn set_period(&mut self, _period: Position) {
        panic!();
    }
    fn unset_period(&mut self) {
        panic!();
    }
    fn get_mark(&self) -> Option<&KmerLoc> {
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
    fn clear_p_extension(&mut self) {
        self.p.clear_extension();
    }
    fn get_i(&self) -> usize {
        self.i
    }
    fn get_d(&self, i: usize) -> &Kmer<u64> {
        &self.d[i]
    }

    fn increment_for_extension(&mut self, ks: &KmerStore) -> Result<()> {
        let b2 = ks.b2_for_p(self.get_p(), false)?;
        dbg_assert!(self.increment(b2));
        Ok(())
    }
    fn dist_if_repetitive(
        &self,
        stored_p: u64,
        mark_p: u64,
        max_dist: Position,
    ) -> Option<Position> {
        None
    }

    /// add twobit to k-mers, update k-mer vec, incr. pos and update ori
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
        // first bit is strand orientation (ori), set according to k-mer.
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
    fn set_mark(&mut self, idx: usize, p: u64, x: usize) {
        dbg_print!("{:<30}<P>", format!("[{:x}] = {:x} | x({})", idx, p, x));
        self.mark.set(idx, p, x);
    }
}
