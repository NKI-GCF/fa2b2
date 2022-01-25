use crate::kmerconst::KmerConst;
use crate::kmerconst::XmerHash;
use crate::kmerstore::KmerStore;
use crate::new_types::{
    extended_position::{ExtPosEtc, EXT_MAX},
    extension::Extension,
    position::Position,
    twobit::{TwoBit, TwoBitDna, TwoBitRcDna},
};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use crate::scope::NO_XMERS;
use crate::xmer_location::XmerLoc;
use anyhow::{anyhow, Result};
use moveslice::Moveslice;
use noodles_fastq as fastq;
use std::cmp::Ordering;

pub struct Mapper<'a> {
    ks: &'a KmerStore,
    kc: &'a KmerConst,
    xmer_loc: [XmerLoc; NO_XMERS], // misschien is deze on the fly uit ks te bepalen?
    xmer_loc_ord: [usize; NO_XMERS],
    pos_lookup: [ExtPosEtc; NO_XMERS],
    pos: Position,
    dna: TwoBitDna,
    rc: TwoBitRcDna,
    i: usize,
    rotation: usize,
    mark_i: usize,
}

impl<'a> Mapper<'a> {
    pub(crate) fn new(ks: &'a KmerStore, kc: &'a KmerConst) -> Self {
        Mapper {
            ks,
            kc,
            xmer_loc: Default::default(),
            xmer_loc_ord: (0..NO_XMERS).collect::<Vec<_>>().try_into().unwrap(),
            pos_lookup: Default::default(),
            pos: Position::default(),
            dna: TwoBitDna::new(0),
            rc: TwoBitRcDna::new(0),
            i: 0,
            rotation: 0,
            mark_i: usize::MAX,
        }
    }
    /// get hashed k-mer, with no compression of the orientation
    fn get_hashed_kmer(&self, is_template: bool) -> usize {
        let seq = if is_template {
            self.dna.to_usize()
        } else {
            self.rc.to_usize()
        };
        self.kc.xmer_hash(seq, self.kc.seed)
    }
    fn get_xmer_loc_ord_leaving_and_inserted_indices(
        &self,
        inserted_value: usize,
    ) -> (usize, usize) {
        let rotation = self.rotation;
        let xmer_loc = &self.xmer_loc;
        let xmer_loc_ord = &self.xmer_loc_ord;
        // The leaving xmer_loc to search is implied by the rotation.
        // another array.
        let leaving_index = xmer_loc_ord
            .binary_search_by(|&i| {
                xmer_loc[i]
                    .idx
                    .cmp(&xmer_loc[rotation].idx)
                    .then(i.cmp(&rotation))
            })
            .unwrap();
        // provide insertion value.
        match xmer_loc_ord
            .binary_search_by(|&i| xmer_loc[i].idx.cmp(&inserted_value).then(i.cmp(&rotation)))
        {
            Ok(v) => (leaving_index, v + 1),
            Err(v) => (leaving_index, v),
        }
    }

    // This function evaluates the number of multimappers.
    // An error indicates an xmer inconsistency: should be a mismatch to ref
    // in the window for this xmer. Read on until next xmer? else, try adjacent to picked median?
    fn se_mapping(&mut self, i: usize) -> Result<()> {
        let xmer_loc = &self.xmer_loc[i];
        // binary search for dupbit 0 status
        let mut last_pos = Position::default();
        for x in 0..=EXT_MAX {
            let hash = xmer_loc.get_hash(self.kc, x);
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

    /// If the median xmer changes, its position may already have been queued for this index. This
    /// happens because a xmer can first be median, then not and later again be median within NO_KMERS.
    ///
    /// Secondly, only pass the first idx & pos for regularly repetitive sequences. Store the
    /// number of repetitions is stored in a hashmap. Repetitions are also deemed 'unworthy' for storage.
    ///
    fn get_worthy_median_xmer_idx(&mut self) -> Option<usize> {
        let mut ret = None;
        let i = self.xmer_loc_ord[self.kc.no_kmers >> 1];
        if i != self.mark_i {
            // use first bits of basepos | ori in array for lookup. Sufficient for within scope of NO_KMERS.
            let scope_idx = self.xmer_loc[i].get_scope_idx();
            if self.pos_lookup[scope_idx] != self.xmer_loc[i].p {
                self.pos_lookup[scope_idx] = self.xmer_loc[i].p;
                ret = Some(i)
            }
        }
        ret
    }
    pub fn updated_median_xmer(&mut self, b2: TwoBit) -> Option<XmerLoc> {
        let mut ret = None;
        if self.i >= self.kc.no_kmers {
            self.update();
            // two marks are added for both strands, and two leave. either can become the new mark.
            if let Some(i) = self.get_worthy_median_xmer_idx() {
                let mut median_xmer = self.xmer_loc[i];

                // update mark so we can skip until next median xmer
                self.set_mark(i);

                // pass through the xmer_loc with hashing undone.
                // FIXME: put the other strand in the idx top bits. We need more alternative XmerLoc types
                // for distinction along with traits !!
                median_xmer.idx = self.kc.xmer_hash(median_xmer.idx, self.kc.seed);
                ret = Some(median_xmer);
            }
        } else {
            self.get_ready();
        }
        self.increment(b2);
        ret
    }
    fn new_xmer_median(&mut self) -> Option<usize> {
        let i = self.pick_mark();
        if self.xmer_loc[i].pos == self.mark.p.pos() {
            None
        } else {
            Some(i)
        }
    }
    fn get_ready(&mut self) {
        // Use the seed to hash. TODO: This hash should be undone before colision detection.
        let ext = self.kc.seed;
        let mut p = ExtPosEtc::from((Extension::from(ext), self.pos));

        let template_hash = self.kc.xmer_hash(self.dna.to_usize(), ext);
        self.xmer_loc[self.rotation].set(template_hash, p);
        self.rotation += 1;

        let reverse_complement_hash = self.kc.xmer_hash(self.dna.to_usize(), ext);
        self.xmer_loc[self.rotation].set(reverse_complement_hash, p);
        self.rotation += 1;
    }
    pub(crate) fn read_record(&mut self, record: fastq::Record) -> Result<()> {
        for b in record.sequence().into_iter() {
            if let Some(b2) = TwoBit::from_u8(*b) {
                if let Some(mark) = self.updated_median_xmer() {
                    self.se_mapping(i)?;
                }
            } else {
                //FIXME: N's ??
            }
        }
        Ok(())
    }
}

impl<'a> Scope for Mapper<'a> {
    /// Two marks are added for both strands, and two leave. Order is maintained in kmer_idx for
    /// median selection.
    fn update(&mut self) {
        // The seed is used to hash the xmer before median selection. A good seed should help select
        // against repetitive sequences as median. On the other hand, the flipping of bits
        // could double the effect of mismatches.
        let ext = self.kc.seed;
        let template_p = ExtPosEtc::from((Extension::from(ext), self.pos));

        // binary search of the 2 insert sites for the pending in xmer_loc_ord.
        for p in [template_p, template_p.get_rc()] {
            // here the
            let inserted_value = self.get_hashed_kmer(p.is_template());
            let (leaving_index, inserted_index) =
                self.get_xmer_loc_ord_leaving_and_inserted_indices(inserted_value);

            match inserted_index.cmp(&leaving_index) {
                Ordering::Less => {
                    self.xmer_loc_ord
                        .moveslice(inserted_index..leaving_index, inserted_index + 1);
                    self.xmer_loc_ord[inserted_index] = self.rotation;
                }
                Ordering::Greater => {
                    self.xmer_loc_ord
                        .moveslice((leaving_index + 1)..inserted_index, leaving_index);
                    self.xmer_loc_ord[inserted_index - 1] = self.rotation;
                }
                Ordering::Equal => {}
            };
            if self.mark_i == leaving_index {
                self.mark_i = usize::MAX
            }
            self.xmer_loc[self.rotation].set(inserted_value, p);
            self.rotation += 1;
        }
        if self.rotation == self.kc.no_kmers {
            self.rotation = 0;
        }
    }
    fn dist_if_repetitive(
        &self,
        _ks: &KmerStore,
        _sp: ExtPosEtc,
        _mp: ExtPosEtc,
    ) -> Option<Position> {
        None
    }

    fn increment(&mut self, b2: TwoBit) {
        // first bit is strand orientation (ori), set according to k-mer.
        self.dna.add(b2, self.kc.dna_topb2_shift);
        self.rc.add(b2, self.kc.rc_mask);
        self.pos.incr()
    }

    fn set_mark(&mut self, i: usize) {
        dbg_print!("{}", self.xmer_loc[i]);
        self.mark_i = i;
    }
}
