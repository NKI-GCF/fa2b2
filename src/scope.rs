use crate::kmerconst::{KmerConst, XmerHash};
use crate::kmerstore::KmerStore;
use crate::new_types::{
    extended_position::ExtPosEtc,
    extension::Extension,
    position::Position,
    twobit::{TwoBit, TwoBitDna, TwoBitRcDna},
};
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use anyhow::{ensure, Result};
use bitvec::prelude::{BitSlice, Lsb0};
use moveslice::Moveslice;
use std::cmp::Ordering;
use std::iter::repeat;

// for this nr of xmers a median will be selected. Two strands, so 16 basepositions + kmer length.
pub struct Scope<'a> {
    kc: &'a KmerConst,
    xmer_loc: Vec<XmerLoc>,
    xmer_loc_ord: Vec<usize>,
    pos_lookup: Vec<ExtPosEtc>,
    mark_i: usize,
    pos: Position,
    dna: TwoBitDna,
    rc: TwoBitRcDna,
    i: usize,
    rotation: usize,
}

/// Processes the sequence head. Filters ambiguous and passes through only median xmers within the
/// scope NO_XMERS.
impl<'a> Scope<'a> {
    pub(crate) fn new(kc: &'a KmerConst) -> Result<Self> {
        ensure!(kc.no_kmers & 1 == 0, "Even: req to store template & revcmp");
        Ok(Scope {
            kc,
            xmer_loc: repeat(XmerLoc::default())
                .take(kc.no_kmers)
                .collect::<Vec<_>>(),
            xmer_loc_ord: (0..kc.no_kmers).collect::<Vec<_>>(),
            pos_lookup: repeat(ExtPosEtc::default())
                .take(kc.no_kmers)
                .collect::<Vec<_>>(),
            mark_i: usize::MAX,
            pos: Default::default(),
            dna: TwoBitDna::new(0),
            rc: TwoBitRcDna::new(0),
            i: 0,
            rotation: 0,
        })
    }
    /// clear for the processing of more sequence.
    pub(crate) fn reset(&mut self) {
        self.i = 0;
        self.rotation = 0;
    }
    pub(crate) fn get_xmer_loc(&self) -> Option<&XmerLoc> {
        self.xmer_loc.get(self.mark_i)
    }
    /// If the median xmer (mark) changes, its position may already have been queued for this index. This
    /// happens because a xmer can first be median, then not and later again be median within NO_KMERS.
    ///
    /// Secondly, only pass the first idx & pos for regularly repetitive sequences. Store the
    /// number of repetitions is stored in a hashmap. Repetitions are also deemed 'unworthy' for storage.

    pub fn updated_median_xmer(&mut self, b2: &BitSlice<u8, Lsb0>) -> Option<XmerLoc> {
        let mut ret = None;
        if self.is_xmer_ready_to_estimate_optima() {
            self.update();
            // two marks are added for both strands, and two leave. either can become the new mark.
            let i = self.xmer_loc_ord[self.kc.no_kmers >> 1];
            if self.rotation == 0 {
                // TODO: set bit in ExtPosEtc to prioritize one in no_kmers against extension.
                // This should allow for smaller kmers and reduce memory.
            }
            if i != self.mark_i || self.rotation == 0 {
                // skip until next median xmer
                self.mark_i = i;
                let test = &self.xmer_loc[i];
                // use first bits of basepos | ori in array for lookup. Sufficient for within scope of NO_KMERS.
                let scope_idx = test.get_scope_idx(self.kc.no_kmers);
                if self.pos_lookup[scope_idx] != test.p {
                    // TODO: use self.mini_kmp to filter out duplicates
                    self.pos_lookup[scope_idx] = test.p;
                    let mut median_xmer = XmerLoc::new(test.idx, test.p);
                    // FIXME: put the other strand in the idx top bits (saves reverse complementing).
                    // Also, we need more alternative XmerLoc types for distinction along with traits !!

                    // pass through the xmer_loc with hashing undone.
                    median_xmer.idx = self.kc.xmer_hash(median_xmer.idx, self.kc.seed);
                    median_xmer.p.clear_extension();
                    dbg_print!("Median xmer: {}", median_xmer);
                    ret = Some(median_xmer);
                }
            }
        } else {
            self.get_ready();
        }
        self.increment(b2);
        ret
    }
    pub(crate) fn is_past_contig(&self) -> bool {
        self.i != 0
    }
    pub(crate) fn get_pos(&self) -> Position {
        self.pos
    }
    fn is_xmer_ready_to_estimate_optima(&self) -> bool {
        self.i >= self.kc.no_kmers
    }

    fn get_ready(&mut self) {
        // Use the seed to hash. TODO: This hash should be undone before colision detection.
        let ext = self.kc.seed;
        let p = ExtPosEtc::from((Extension::from(ext), self.pos));

        let template_hash = self.kc.xmer_hash(self.dna.to_usize(), ext);
        self.xmer_loc[self.rotation].set(template_hash, p);

        let reverse_complement_hash = self.kc.xmer_hash(self.dna.to_usize(), ext);
        self.xmer_loc[self.rotation + 1].set(reverse_complement_hash, p);
        self.rotation += 2;
        self.i += 2;
        if self.rotation == self.kc.no_kmers {
            self.rotation = 0;
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
        ks: &KmerStore,
        stored_p: ExtPosEtc,
        min_p: ExtPosEtc,
    ) -> Option<Position> {
        // FIXME: the contig check was removed. for the upper bound that makes sense for head
        // scope, as upperbound is pos_max rather than previous
        let stored_pos = stored_p.pos();
        let mark_pos = min_p.pos();
        dbg_assert!(mark_pos > stored_pos);
        if ks.is_on_last_contig(stored_pos) {
            let dist = mark_pos - stored_pos;
            if dist < ks.rep_max_dist {
                return Some(dist);
            }
        }
        None
    }

    /// add twobit to k-mers, increment pos for next median selection
    fn increment(&mut self, b2: &BitSlice<u8, Lsb0>) {
        self.dna.add(TwoBit::from(b2), self.kc.dna_topb2_shift);
        self.rc.add(TwoBit::from(b2), self.kc.rc_mask);
        self.pos.incr()
    }
}
