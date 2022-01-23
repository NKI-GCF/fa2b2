use crate::kmerconst::{KmerConst, XmerHash};
use crate::kmerstore::KmerStore;
use crate::new_types::{
    extended_position::ExtPosEtc,
    extension::Extension,
    position::Position,
    twobit::{TwoBit, TwoBitDna, TwoBitRcDna},
};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use crate::xmer_location::XmerLoc;
use moveslice::Moveslice;
use std::cmp::Ordering;
//std:sync::{mpsc, Arc},
//std:thread,

const NTHREADS: usize = 16;

// for this nr of xmers a median will be selected. Two strands, so 16 basepositions + kmer length.
pub(crate) const NO_XMERS: usize = 32;

pub struct HeadScope<'a> {
    kc: &'a KmerConst,
    xmer_loc: [XmerLoc; NO_XMERS],
    xmer_loc_ord: [usize; NO_XMERS],
    pos_lookup: [ExtPosEtc; NO_XMERS],
    mark_i: usize,
    pos: Position,
    dna: TwoBitDna,
    rc: TwoBitRcDna,
    i: usize,
    rotation: usize,
}

/// Processes the sequence head. Filters ambiguous and passes through only median xmers within the
/// scope NO_XMERS.
impl<'a> HeadScope<'a> {
    pub(crate) fn new(kc: &'a KmerConst) -> Self {
        HeadScope {
            kc,
            xmer_loc: Default::default(),
            xmer_loc_ord: (0..NO_XMERS).collect::<Vec<_>>().try_into().unwrap(),
            pos_lookup: Default::default(),
            mark_i: usize::MAX,
            pos: Default::default(),
            dna: TwoBitDna::new(0),
            rc: TwoBitRcDna::new(0),
            i: 0,
            rotation: 0,
        }
    }
    /// clear sequence accounting for the processing of more sequence.
    pub(crate) fn reset_get_pos(&mut self) -> Position {
        self.i = 0;
        self.rotation = 0;
        self.pos
    }
    pub(crate) fn get_xmer_loc(&self) -> Option<&XmerLoc> {
        self.xmer_loc.get(self.mark_i)
    }
    pub(crate) fn is_coding_sequence_pending(&self, ks: &KmerStore) -> bool {
        ks.get_pos() != self.pos
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

    pub fn set_mark(&mut self) {}

    /// If the median xmer changes, its position may already have been queued for this index. This
    /// happens because a xmer can first be median, then not and later again be median within NO_KMERS.
    ///
    /// Secondly, only pass the first idx & pos for regularly repetitive sequences. Store the
    /// number of repetitions is stored in a hashmap. Repetitions are also deemed 'unworthy' for storage.
    ///

    pub fn updated_median_xmer(&mut self, b2: TwoBit) -> Option<usize> {
        let mut ret = None;
        if self.i >= self.kc.no_kmers {
            self.update();
            // two marks are added for both strands, and two leave. either can become the new mark.
            let i = self.xmer_loc_ord[self.kc.no_kmers >> 1];
            if i != self.mark_i {
                let test = &self.xmer_loc[i];
                // use first bits of basepos | ori in array for lookup. Sufficient for within scope of NO_KMERS.
                let scope_idx = test.get_scope_idx();
                if self.pos_lookup[scope_idx] != test.p {
                    // TODO: use self.mini_kmp to filter out duplicates
                    self.pos_lookup[scope_idx] = test.p;
                    ret = Some(i);
                }
            }
        } else {
            self.get_ready();
        }
        self.increment(b2);
        ret
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
}
impl<'a> Scope for HeadScope<'a> {
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

    fn set_mark(&mut self, i: usize) {
        dbg_print!("{}", self.xmer_loc[i]);
        self.mark_i = i;
    }

    /// add twobit to k-mers, increment pos for next median selection
    fn increment(&mut self, b2: TwoBit) {
        self.dna.add(b2, self.kc.dna_topb2_shift);
        self.rc.add(b2, self.kc.rc_mask);
        self.pos.incr()
    }
}

/*pub struct OldHeadScope<'a> {
    kc: &'a KmerConst,
    p: ExtPosEtc,
    mark: XmerLoc,
}

impl<'a> OldHeadScope<'a> {
    pub(crate) fn new(kc: &'a KmerConst) -> Self {
        OldHeadScope {
            kc,
            p: ExtPosEtc::default(),
            mark: XmerLoc::default(),
        }
    }
    pub(crate) fn get_pos(&self) -> Position {
        self.p.pos()
    }

    pub(crate) fn is_coding_sequence_pending(&self, ks: &KmerStore) -> bool {
        ks.get_pos() != self.get_pos()
    }

    // .i & .p increments en kmer .d[] updates vinden plaats.
    pub(crate) fn complete_and_update_mark(
        &mut self,
        ks: &mut KmerStore,
        b2: TwoBit,
    ) -> Result<()> {
        if self.period.is_set() {
            let pd = self.period;
            let pos = self.p.pos();
            dbg_assert!(pd <= pos, "{} {}", pd, self.p);
            if ks.b2_for_pos(pos - pd, true) == b2 {
                self.update_repetitive(ks, pd);
            } else {
                self.period.to_default();
            }
        }
        if self.update() {
            // one mark is added, and one leaving. both influence mark (and order).
            let i = self.pick_mark();
            if self.d[i].pos == self.mark.p.pos() {
                return Ok(());
            }
            let mut mark = self.d[i].get_hash_and_p(self.kc, 0);
            self.set_mark(&mark);
            let orig_pos = mark.p.pos();
            // this seems to be hotlooping
            while self.try_store_mark(ks, &mut mark)? {
                if mark.p.pos() != orig_pos {
                    if self.kc.extend_xmer(&mut mark).is_ok() {
                        // extending some pase baseidx. TODO: if frequently the same recurs,
                        // it might be worthwhile to store the reverse complement in a temp
                        // we should not store a past index in self.mark.p !!
                    } else {
                        dbg_print!("couldn't extend: {} ..?", mark);
                        break;
                    }
                } else {
                    if mark.p.extend().is_err() {
                        break;
                    }
                    mark = self.d[i].get_hash_and_p(self.kc, mark.p.x());
                    self.set_mark(&mark);
                }
            }
        }
        self.increment(b2);
        Ok(())
    }

    fn try_store_mark(&mut self, ks: &mut KmerStore, mark: &mut XmerLoc) -> Result<bool> {
        if self.period.is_set() && ks.kmp[mark.idx].is_set() {
            ks.kmp[mark.idx].set_repetitive();
            return Ok(false);
        }
        let old_stored_p = ks.kmp[mark.idx];

        if old_stored_p.is_zero() {
            ks.set_kmp(&mark);
            return Ok(false);
        }
        if old_stored_p.is_replaceable_by(mark.p) {
            if old_stored_p.pos() == mark.p.pos() {
                // set and already mark.p. Leave the bit states.
                return Ok(false);
            }
            ks.set_kmp(&mark);
            if old_stored_p.x() == mark.p.x() {
                // same extension means same base k-mer origin. this is a duplicate.
                ks.kmp[mark.idx].mark_more_recurs_upseq();
            }
            dbg_print!("{} -> ?", mark);
            mark.p = old_stored_p; // to be extended next
        } else if old_stored_p.extension() == mark.p.extension() {
            // Note: same extension and hash means same k-mer origin: identical k-mer sequence.

            // If a kmer occurs multiple times within an extending readlength (repetition),
            // only the first gets a position. During mapping this should be kept in mind.
            if let Some(dist) = self.dist_if_repetitive(ks, old_stored_p, mark.p) {
                dbg_assert!(dist < self.p.pos() || self.p.is_zero());
                self.period = dist;
                ks.kmp[mark.idx].set_repetitive();
                return Ok(false);
            }
            ks.kmp[mark.idx].mark_more_recurs_upseq();
        }
        // collision between hashes, the one in mark.p will be extended and tried again.
        Ok(true)
    }
}

impl<'a> Scope for OldHeadScope<'a> {
    /// add twobit to k-mers, update k-mer vec, increment pos and update orientation
    /// true if we have at least one kmer.
    fn update(&mut self) -> bool {
        if self.i >= self.kc.kmerlen {
            let old_d = self.d[self.rotation];
            self.rotation += 1;
            if self.rotation == self.kc.no_kmers {
                self.rotation = 0;
            }
            self.d[self.rotation] = old_d;
            self.d[self.rotation].pos = self.p.pos();
            true
        } else {
            self.i += 1;
            false
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

    fn set_mark(&mut self, mark: &XmerLoc) {
        dbg_print!("{}", mark);
        self.mark = *mark;
    }

    fn increment(&mut self, b2: TwoBit) {
        // first bit is strand bit, set according to kmer orientation bit.
        self.p.set_ori(self.d[self.rotation].update(self.kc, b2));
        self.p.incr_pos();
    }
}*/
