use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::{ExtPosEtc, EXT_MAX};
use crate::new_types::{
    position::Position,
    twobit::{ThreeBit, TwoBit},
    xmer::Xmer,
};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use crate::xmer_location::XmerLoc;
use anyhow::{anyhow, Result};
use noodles_fastq as fastq;

pub struct Mapper<'a> {
    ks: &'a KmerStore,
    kc: &'a KmerConst,
    p: ExtPosEtc,
    i: usize,
    mod_i: usize,
    mark: XmerLoc,
    d: Vec<Xmer>, // misschien is deze on the fly uit ks te bepalen?
    z: Vec<usize>,
}

impl<'a> Mapper<'a> {
    pub(crate) fn new(ks: &'a KmerStore, kc: &'a KmerConst) -> Self {
        Mapper {
            ks,
            kc,
            p: ExtPosEtc::default(),
            i: 0,
            mod_i: 0,
            mark: XmerLoc::new(usize::max_value(), ExtPosEtc::default()),
            d: vec![Xmer::new(); kc.no_kmers],
            z: (0..kc.no_kmers).collect(),
        }
    }

    // This function evaluates the number of multimappers.
    // An error indicates an xmer inconsistency: should be a mismatch to ref
    // in the window for this xmer. Read on until next xmer? else, try adjacent to picked median?
    fn se_mapping(&mut self, i: usize) -> Result<()> {
        let xmer = &self.d[i];
        // binary search for dupbit 0 status
        let mut last_pos = Position::default();
        for x in 0..=EXT_MAX {
            let hash = xmer.get_hash(self.kc, x);
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

    pub(crate) fn read_record(&mut self, record: fastq::Record) -> Result<()> {
        for b in record.sequence().iter() {
            if let Ok(b2) = TwoBit::try_from(ThreeBit::from(b)) {
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
        self.p.set_ori(self.d[self.mod_i].update(self.kc, b2));
        self.p.incr_pos();
        self.i += 1;
    }
    fn set_mark(&mut self, mark: &XmerLoc) {
        dbg_print!("{} (mapping mark)", mark);
        self.mark = *mark;
    }
}
