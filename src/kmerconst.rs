use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{
    position::{BasePos, Position},
    xmer::xmer_hash,
};
use crate::rdbg::STAT_DB;
use anyhow::Result;
use num::{FromPrimitive, PrimInt};
use std::cmp;
use std::mem::size_of;

pub struct KmerConst {
    pub no_kmers: usize,
    pub kmerlen: usize,
    pub bitlen: usize,
    pub read_len: usize,
    pub extent: Vec<usize>,
}

pub trait RevCmp<T: PrimInt + FromPrimitive> {
    fn revcmp(self, kmerlen: u32) -> T;
}

/// create bitmask. e.g. dvm::<u32>(0xf0  0xff) => 0xf0_f0_f0_f0
#[inline]
fn dvm<T: PrimInt + FromPrimitive>(numerator: u32, divisor: u32) -> T {
    let base = T::max_value() / T::from_u32(divisor).unwrap();
    T::from_u32(numerator).unwrap() * base
}

#[macro_export]
macro_rules! implement_revcmp { ($($ty:ty),*) => ($(
    /// give twobit reverse complent for given kmerlen
    impl RevCmp<$ty> for $ty {
        #[inline]
        fn revcmp(self, kmerlen: u32) -> $ty {
            let mut seq = self.swap_bytes() ^ dvm::<$ty>(2, 3);

            let left_nibbles = (seq & dvm::<$ty>(0xf, 0xff)) << 4;
            let right_nibbles = (seq & dvm::<$ty>(0xf0, 0xff)) >> 4;
            seq = left_nibbles | right_nibbles;

            let all_left_two_bits = (seq & dvm::<$ty>(0x3, 0xf)) << 2;
            let all_right_two_bits = (seq & dvm::<$ty>(0xc, 0xf)) >> 2;
            seq = all_left_two_bits | all_right_two_bits;

            seq >> (size_of::<$ty>() * 8 - usize::try_from(kmerlen).unwrap() * 2)
        }
    }
    )*)
}
implement_revcmp!(usize, u64);

impl KmerConst {
    pub(crate) fn new(genomesize: usize, read_len: usize) -> Self {
        // bit width, required to store all (cumulative) genomic positions, is used as len

        let mut bitlen = genomesize.next_power_of_two().trailing_zeros() as usize;
        if (bitlen & 1) == 1 {
            // must be even.
            bitlen += 1
        }
        KmerConst::from_bitlen(bitlen, read_len)
    }

    pub(crate) fn from_bitlen(bitlen: usize, read_len: usize) -> Self {
        let kmerlen = bitlen / 2;

        // e.g. with a RL 4 & KL 2: (0,1), (1,2), (2,3) => 3 kmers.
        let no_kmers = read_len - kmerlen + 1;

        // generate all combinations of 2 kmer positions for extension, high and low bits.
        // overlap and combination duplicates are allowed, those will result respectively
        // in no xor-hash and a complemented xor-hash with the max for index.
        let extent: Vec<usize> = (0..kmerlen).collect();

        dbg_restart!("read_len: {}, kmerlen: {}\n--", read_len, kmerlen);
        if !cfg!(debug_assertions) {
            eprintln!("read_len: {}, kmerlen: {}", read_len, kmerlen);
        }
        KmerConst {
            no_kmers,
            kmerlen,
            bitlen,
            read_len,
            extent,
        }
    }

    // same hash function, from Xmer.
    pub(crate) fn get_next_xmer(
        &self,
        orig_hash: usize,
        mut p: ExtPosEtc,
    ) -> Result<(usize, ExtPosEtc)> {
        let old_x = p.x();
        p.extend()?;

        let k = self.kmerlen as u32;
        let overbit = 1 << (k * 2 - 1);

        // in get_hash_and_p() bits are flipped if the highest bit was set.
        //
        let mut hash = xmer_hash(orig_hash, old_x, k);

        if hash < hash.revcmp(k) {
            // XXX: why is this not the inverse ??

            //then flipped, yes: orig_hash here !!
            hash = xmer_hash(!orig_hash & (overbit | (overbit - 1)), old_x, k);
        }

        // set to idx for next extension; x is incremented:
        hash = xmer_hash(hash, p.x(), k);
        if (hash & overbit) == 0 {
            Ok((hash, p))
        } else {
            Ok(((overbit - 1) & !hash, p))
        }
    }

    pub(crate) fn get_kmer_boundaries(
        &self,
        pos: Position,
        contig: (Position, Position),
    ) -> (Position, Position) {
        let read_len = Position::from(BasePos::from(self.read_len));
        let lower = if pos > read_len {
            pos - read_len
        } else {
            Position::zero()
        };
        (
            cmp::max(lower, contig.0),
            cmp::min(pos + Position::from(BasePos::from(self.no_kmers)), contig.1),
        )
    }
}
