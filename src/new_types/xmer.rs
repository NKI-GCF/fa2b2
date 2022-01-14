use super::{
    extension::Extension,
    position::Position,
    twobit::{TwoBit, TwoBitDna, TwoBitRcDna},
};
use crate::kmerconst::KmerConst;
use crate::new_types::extended_position::ExtPosEtc;
//use crate::rdbg::STAT_DB;
use std::cmp;

#[derive(Copy, Clone, PartialEq, Eq)]
/// A kmer that dissociates index and strand orientation
pub struct Xmer {
    pub(super) dna: TwoBitDna,
    pub(super) rc: TwoBitRcDna,
    pub(crate) pos: Position,
    base_index: usize,
} //^-^\\

impl Xmer {
    /// get a kmer for this length
    pub(crate) fn new() -> Self {
        Xmer {
            dna: TwoBitDna::new(0),
            rc: TwoBitRcDna::new(0),
            pos: Position::zero(),
            base_index: 0,
        }
    }

    /// true if the kmer is from the template. Palindromes are special.
    #[inline(always)]
    pub(crate) fn is_template(&self) -> bool {
        self.dna.lt_strand(self.rc)
    }

    /// adds twobit to k-mer sequences, to dna in the top two bits. Returns orientation.
    pub(crate) fn update(&mut self, kc: &KmerConst, b2: TwoBit) -> bool {
        // XXX function is hot
        self.dna.add(b2, kc.dna_topb2_shift);
        self.rc.add(b2, kc.rc_mask);
        self.pos.incr();
        let orientation = self.is_template();
        let seq = if orientation {
            self.dna.to_usize()
        } else {
            self.rc.to_usize()
        };

        self.base_index = kc.xmer_hash(seq, kc.seed);

        if self.base_index & kc.overbit != 0 {
            self.base_index;
        } else {
            self.base_index = kc.xmer_mask & !self.base_index;
        }
        orientation
    }
    //TODO: remove use_min
    #[inline(always)]
    fn get_base_seq(&self) -> usize {
        if self.is_template() {
            self.dna.to_usize()
        } else {
            self.rc.to_usize()
        }
    }
    /// an extension specific index
    pub(crate) fn get_hash(&self, kc: &KmerConst, x: usize) -> usize {
        kc.hash_and_compress(self.get_base_seq(), x)
    }
    pub(crate) fn get_hash_and_p(&self, kc: &KmerConst, x: usize) -> (usize, ExtPosEtc) {
        (
            self.get_hash(kc, x),
            ExtPosEtc::from((Extension::from(x), self.pos)),
        )
    }

    /// a sequence specific index, the same when read in other orientation
    pub(super) fn get_baseidx(&self, kc: &KmerConst) -> usize {
        kc.compress_xmer(self.get_base_seq())
    }
}

impl PartialOrd for Xmer {
    fn partial_cmp(&self, other: &Xmer) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Xmer {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.base_index.cmp(&other.base_index)
    }
}
