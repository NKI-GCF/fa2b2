use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::{ExtPosEtc, KmerLoc, EXT_MAX};
use crate::new_types::{
    position::Position,
    twobit::{ThreeBit, TwoBit},
    xmer::Xmer,
};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::{anyhow, Result};
use noodles_fastq as fastq;
use std::fmt;

pub struct Mapper<'a> {
    ks: &'a KmerStore,
    kc: &'a KmerConst,
    p: ExtPosEtc,
    i: usize,
    mod_i: usize,
    mark: KmerLoc,
    d: Vec<Xmer>, // misschien is deze on the fly uit ks te bepalen?
    z: Vec<usize>,
}

impl<'a> Mapper<'a> {
    pub fn new(ks: &'a KmerStore, kc: &'a KmerConst) -> Self {
        Mapper {
            ks,
            kc,
            p: ExtPosEtc::zero(),
            i: 0,
            mod_i: 0,
            mark: KmerLoc::new(usize::max_value(), ExtPosEtc::zero()),
            d: vec![Xmer::new(kc.kmerlen as u32); kc.no_kmers],
            z: (0..kc.no_kmers).collect(),
        }
    }

    // This function evaluates the number of multimappers.
    // An error indicates an xmer inconsistency: should be a mismatch to ref
    // in the window for this xmer. Read on until next xmer? else, try adjacent to picked median?
    fn se_mapping(&mut self, i: usize) -> Result<()> {
        let xmer = &self.d[i];
        // binary search for dupbit 0 status
        let mut last_pos = Position::zero();
        for x in 0..=EXT_MAX {
            let hash = xmer.get_hash(x);
            let test_p = self.ks.kmp[hash];
            if test_p.x() == x {
                let pos = test_p.pos();
                // Because we haven't seen a non-DUPLICATE, we know there should follow more.
                // x and hash lead to one baseindex. That should be position ordered.
                dbg_assert!(last_pos < pos, "Xmer invalid by ref coordinate order.");

                last_pos = pos;
                if test_p.is_repetitive() {}
                if test_p.is_last_on_ref() {
                    // no DUPLICATE: => last sequence that corresponded with this xmer.
                    break;
                }
            } else if test_p.x() < x {
                return Err(anyhow!("Xmer invalid by non-pos."));
            }
            // test_p.x() > x can happen: collisions
            // since the other extension was greater, this xmer & p got extended.
        }
        Ok(())
    }
    fn new_xmer_median(&mut self) -> Option<usize> {
        let i = self.pick_mark();
        if self.d[i].pos == self.mark.p.pos() {
            None
        } else {
            Some(i)
        }
    }

    pub fn read_record(&mut self, record: fastq::Record) -> Result<()> {
        for b3 in record.sequence().iter().map(ThreeBit::from) {
            if let Some(b2) = b3.as_twobit_if_not_n() {
                if self.update() {
                    // we weten extension op voorhand.
                    if let Some(i) = self.new_xmer_median() {
                        self.se_mapping(i)?;
                    }
                }
                self.increment(b2)
            }
        }
        Ok(())
    }
}

impl<'a> Scope for Mapper<'a> {
    fn get_kc(&self) -> &KmerConst {
        self.kc
    }
    fn get_d(&self, i: usize) -> &Xmer {
        &self.d[i]
    }
    fn pick_mark(&mut self) -> usize {
        let med = self.kc.no_kmers >> 1;
        let i = self
            .z
            .select_nth_unstable_by(med, |&a, &b| self.d[a].cmp(&self.d[b]))
            .1;
        *i
    }
    fn dist_if_repetitive(
        &self,
        _ks: &KmerStore,
        _sp: ExtPosEtc,
        _mp: ExtPosEtc,
    ) -> Option<Position> {
        None
    }

    /// add twobit to k-mers, update k-mer vec, incr. pos and update ori
    /// true if we have at least one kmer.
    fn update(&mut self) -> bool {
        // XXX: function is hot
        if self.i >= self.kc.kmerlen {
            let old_d = self.d[self.mod_i];
            self.mod_i += 1;
            if self.mod_i == self.kc.no_kmers {
                self.mod_i = 0;
            }
            self.d[self.mod_i] = old_d;
            self.d[self.mod_i].pos = self.p.pos();
            true
        } else {
            false
        }
    }
    fn increment(&mut self, b2: TwoBit) {
        // first bit is strand orientation (ori), set according to k-mer.
        self.p.set_ori(self.d[self.mod_i].update(b2));
        self.p.incr_pos();
        self.i += 1;
    }
    fn set_mark(&mut self, idx: usize, p: ExtPosEtc) {
        dbg_print!("[{:x}] = {:?}", idx, p);
        self.mark.set(idx, p);
    }
}
