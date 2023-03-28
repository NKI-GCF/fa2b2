// Roel Kluin, 2023, GPL v3

use crate::new_types::position::{BasePos, PosRange, Position};
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use anyhow::Result;
use num::{FromPrimitive, PrimInt};

pub struct KmerConst {
    pub(crate) no_kmers: usize,
    pub(crate) kmerlen: usize,
    pub(crate) bitlen: usize,
    pub(crate) read_len: usize,
    pub(crate) dna_topb2_shift: u32,
    pub(crate) rc_mask: u64,
    // put these in kc:
    pub(crate) half_mask: usize,
    pub(crate) overbit: usize,
    pub(crate) xmer_mask: usize,
    pub(crate) over_mask: usize,
    // XXX!seed influences order of selection, a good seed may improve indexing / mapping.
    // XXX!indexing and mapping should use the same seed !!
    pub(crate) seed: usize,
}

pub trait RevCmp<T: PrimInt + FromPrimitive> {
    fn revcmp(self, kmerlen: usize) -> T;
}

/// create bitmask. e.g. dvm::<u32>(0xf0  0xff) => 0xf0_f0_f0_f0
#[inline]
fn dvm<T: PrimInt + FromPrimitive>(numerator: u32, divisor: u32) -> T {
    let base = T::max_value() / T::from_u32(divisor).unwrap();
    T::from_u32(numerator).unwrap() * base
}
//TODO: lookup table per byte?
macro_rules! implement_revcmp { ($($ty:ty),*) => ($(
    /// give twobit reverse complent for given kmerlen
    impl RevCmp<$ty> for $ty {
        #[inline]
        fn revcmp(self, kmerlen: usize) -> $ty {
            let mut seq = self.swap_bytes() ^ dvm::<$ty>(2, 3);

            let left_nibbles = (seq & dvm::<$ty>(0xf, 0xff)) << 4;
            let right_nibbles = (seq & dvm::<$ty>(0xf0, 0xff)) >> 4;
            seq = left_nibbles | right_nibbles;

            let all_left_two_bits = (seq & dvm::<$ty>(0x3, 0xf)) << 2;
            let all_right_two_bits = (seq & dvm::<$ty>(0xc, 0xf)) >> 2;
            seq = all_left_two_bits | all_right_two_bits;

            seq >> (<$ty>::BITS - (kmerlen as u32) * 2)
        }
    }
    )*)
}
implement_revcmp!(usize, u64);

pub trait XmerHash {
    fn get_kmerlen(&self) -> usize;
    fn get_overbit(&self) -> usize;
    /// a reversible hash, no flip to exclude top bit.
    #[inline(always)]
    fn xmer_hash(&self, idx: usize, seed: usize) -> usize {
        let k = self.get_kmerlen();
        idx ^ ((idx & !seed & ((1 << k) - 1)) << k | ((idx >> k) & seed))
    }
    /// after compression, the hash is reversible but it may also require flipping bits. See extend_xmer()
    fn hash_and_compress(&self, seq: usize, x: usize) -> usize {
        self.compress_xmer(self.xmer_hash(seq, x))
    }
    /// unhash and uncompress to get the k-mer that was selected, sequence or reverse complement.
    fn unhash_and_uncompress_to_kmer(&self, hash: usize, x: usize) -> usize {
        // in compress_xmer() bits are flipped if the highest bit was set.
        let flip = self.xmer_hash(hash, x);
        if flip < flip.revcmp(self.get_kmerlen()) {
            // XXX: why is this not the inverse ??

            //then flipped, yes: mark.idx here !!
            let overbit = self.get_overbit();
            self.xmer_hash((overbit | (overbit - 1)) & !hash, x)
        } else {
            flip
        }
    }

    // same hash function, from Xmer.
    // zou de orientation bit achterwegen kunnen laten, dan zijn er in deze hot loop minder operaties..
    fn extend_xmer(&self, mark: &mut XmerLoc) -> Result<()> {
        let old_x = mark.p.x();
        // early extension, because it may fail if already extended at max.
        mark.p.extend()?;

        let hash = self.unhash_and_uncompress_to_kmer(mark.idx, old_x);
        // Now the compressed hash is undone. Hash it again but this time for the next extension.
        mark.idx = self.hash_and_compress(hash, mark.p.x());
        Ok(())
    }
    /// compress base_seq or xmer_hash from xmer module
    fn compress_xmer(&self, v: usize) -> usize {
        let overbit = self.get_overbit();
        if v & overbit == 0 {
            v
        } else {
            (overbit - 1) & !v
        }
    }
}

impl KmerConst {
    pub(crate) fn new(genomesize: u64, read_len: u16, seed: u16) -> Self {
        // bit width, required to store all (cumulative) genomic positions, is used as len

        // FIXME: met een priority bit per 2*16 kmers zou 2 basen = 4 bits minder ook moeten werken.
        let mut bitlen = genomesize.next_power_of_two().trailing_zeros() as u8;
        if (bitlen & 1) == 1 {
            // must be even.
            bitlen += 1
        }
        KmerConst::from_bitlen(bitlen, read_len, seed)
    }

    pub(crate) fn from_bitlen(bitlen: u8, read_len: u16, seed: u16) -> Self {
        // TODO: Do some sanity checks on sizes here
        let bitlen = usize::from(bitlen);
        let read_len = usize::from(read_len);
        let seed = usize::from(seed);

        let kmerlen = bitlen / 2;
        dbg_assert!(kmerlen > 0);
        dbg_assert!(read_len >= kmerlen);

        // e.g. with a RL 4 & KL 2: (0,1), (1,2), (2,3) => 3 kmers + 3x rc.
        let no_kmers = (read_len + 1 - kmerlen) * 2;

        dbg_init!(
            "read_len: {read_len}, kmerlen: {kmerlen}, no_kmers: {no_kmers} (including rc)\n--"
        );
        if !cfg!(debug_assertions) {
            eprintln!(
                "read_len: {read_len}, kmerlen: {kmerlen}, no_kmers: {no_kmers} (including rc)"
            );
        }
        let dna_topb2_shift = (u32::try_from(kmerlen).unwrap() << 1) - 2;
        let overbit = 1_usize
            .checked_shl(dna_topb2_shift + 1)
            .expect("k-mer shift");
        let xmer_mask = overbit - 1;
        KmerConst {
            no_kmers,
            kmerlen,
            bitlen,
            read_len,
            dna_topb2_shift,
            rc_mask: (1 << dna_topb2_shift) - 1,
            half_mask: (1 << kmerlen) - 1,
            overbit,
            xmer_mask,
            over_mask: overbit | xmer_mask,
            seed,
        }
    }
    pub(crate) fn get_kmer_boundaries(&self, pos: Position, mut contig: PosRange) -> PosRange {
        let read_len = Position::from(BasePos::from(self.read_len as u64));

        if pos > read_len {
            contig.bound_lower(pos - read_len)
        }
        contig.bound_upper(pos + Position::from_basepos((self.no_kmers / 2) as u64));
        contig
    }
}

impl XmerHash for KmerConst {
    fn get_kmerlen(&self) -> usize {
        self.kmerlen
    }
    fn get_overbit(&self) -> usize {
        self.overbit
    }
}
