use std::mem::size_of;
use num_traits::PrimInt;
use num::{FromPrimitive,ToPrimitive,Unsigned};

#[derive(Copy, Clone)]
/// A kmer that dissociates index and strand orientation
pub struct Kmer<T> {
	pub dna: T,
	pub rc: T,
	max: T,
	topb2: T,
} //^-^\\

/// return true if dna sequence has the bit set for the first bit that differs between complements
pub fn test_template(dna: u64, rc: u64) -> bool {
	// The first bit that differs between complements, isolated with the '& wrapping_neg()',
	// the devbit, dictates the orientation. If none is set, a palindrome, use bit 1.
	let deviant = dna ^ rc;
	(dna & if deviant != 0 {deviant.wrapping_neg() & deviant} else {1}) != 0
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


macro_rules! implement_revcmp { ($($ty:ty),*) => ($(
	/// give twobit reverse complent for given kmerlen
	impl RevCmp<$ty> for $ty {
		#[inline]
		fn revcmp(self, kmerlen: usize) -> $ty {
			let mut seq = self.swap_bytes() ^ dvm::<$ty>(2, 3);
			seq = ((seq & dvm::<$ty>(0xf0, 0xff)) >> 4) | ((seq & dvm::<$ty>(0xf, 0xff)) << 4);
			seq = ((seq & dvm::<$ty>(0xc, 0xf)) >> 2) | ((seq & dvm::<$ty>(0x3, 0xf)) << 2);
			seq >> (size_of::<$ty>() * 8 - kmerlen * 2)
		}
	}
	)*)
}

implement_revcmp!(u8, u16, u32, u64, u128, usize);

impl<T> Kmer<T>
	where T: Unsigned + FromPrimitive + ToPrimitive
  {
	/// get a kmer for this length
	pub fn new(kmerlen: u32) -> Self {
		let bitlen = kmerlen * 2;
		let topb2 = bitlen - 2;
		Kmer {
			dna: <T>::zero(),
			rc: <T>::zero(),
			max: T::from_u64(1 << (topb2 + 1)).unwrap_or_else(|| panic!("{} twobits don't fit in a {}-bits kmer.", kmerlen, u64::max_value().count_ones())),
			topb2: T::from_u32(topb2).unwrap(),
		}
	}
	/// return a new kmer for given index, length and orientation.
	pub fn from_idx(index: usize, kmerlen: u32, ori: bool) -> Self {
		let bitlen = kmerlen * 2;
		let topb2 = bitlen - 2;
		let mut dna = usize::to_u64(&index).unwrap();
		let mut rc = dna.revcmp(kmerlen as usize);
		if !test_template(dna, rc) {
			let overbit = 1 << (topb2 + 1);
			let overmask = overbit | (overbit - 1);
			dna ^= overmask;
			rc ^= overmask;
		}
		let (dna, rc) = if ori {(dna, rc)} else {(rc, dna)};
		Kmer {
			dna: T::from_u64(dna).unwrap(),
			rc: T::from_u64(rc).unwrap(),
			max: T::from_u64(1 << (topb2 + 1)).unwrap_or_else(|| panic!("{} twobits don't fit in a {}-bits kmer.", kmerlen, u64::max_value().count_ones())),
			topb2: T::from_u32(topb2).unwrap(),
		}
	}
	/// adds twobit to kmer sequences, to .dna in the top two bits.
	pub fn add(&mut self, b2: u8) {
		debug_assert!(b2 <= 3);
		let topb2 = T::to_u64(&self.topb2).unwrap();
		let topless = (1 << topb2) - 1;
		let dna = T::to_u64(&self.dna).unwrap() >> 2;
		let rc = (T::to_u64(&self.rc).unwrap() & topless) << 2;
		self.dna = T::from_u64(dna | (u64::from(b2) << topb2)).unwrap();
		self.rc = T::from_u64(rc ^ 2 ^ u64::from(b2)).unwrap();
	}
	/// true if the kmer is from the template.
	pub fn is_template(&self) -> bool {
		test_template(T::to_u64(&self.dna).unwrap(), T::to_u64(&self.rc).unwrap())
	}
	/// return whether orientation needs to be changed
	pub fn update(&mut self, b2: u8, for_template: bool) -> bool {
		// 2 is added for next pos; orientation is set in first bit.
		if for_template {
			self.add(b2);
			self.is_template()
		} else {
			self.add(b2 ^ 2);
			!self.is_template()
		}
	}
	/// return an index specific per sequence but the same for the other orientation
	pub fn get_idx(&self, for_template: bool) -> usize {
		let seq = T::to_usize(if self.is_template() == for_template {&self.dna} else {&self.rc}).unwrap();

		// flipped if the top bit is set, to reduce size.
		let overbit = 1 << (T::to_u64(&self.topb2).unwrap() + 1);
		if (seq & overbit) == 0 {seq} else {(overbit - 1) & !seq}
	}
}

#[cfg(test)]
mod tests {
    use super::{Kmer,RevCmp};
    use std::{cmp,iter::once};
    use rand::{thread_rng, Rng};
    #[test]
    fn test_u64() {
        let mut kmer: Kmer<u64> = Kmer::new(32);
        for i in 0..32 {
            kmer.add(i & 3);
        }
        assert_eq!(kmer.dna, 0xE4E4E4E4E4E4E4E4);  // GTCAGTCAGTCAGTCA => 3210321032103210 (in 2bits)
        assert_eq!(kmer.rc, 0xB1B1B1B1B1B1B1B1);   // xor 0xaaaaaaaa and reverse per 2bit
        assert_eq!(kmer.is_template(), false); // first devbit is 1, ori is set in rc, so false
        assert_eq!(kmer.get_idx(true), 0x4E4E4E4E4E4E4E4E); // highest bit is set, so flipped.
    }
    #[test]
    fn test_u32() {
        let mut kmer: Kmer<u32> = Kmer::new(16);
        for i in 0..16 {
            kmer.add(i & 3);
        }
        assert_eq!(kmer.dna, 0xE4E4E4E4);
        assert_eq!(kmer.rc, 0xB1B1B1B1);
        assert_eq!(kmer.is_template(), false);
        assert_eq!(kmer.get_idx(true), 0x4E4E4E4E);
    }
    #[test]
    fn test_u8() {
        let mut kmer: Kmer<u8> = Kmer::new(4);
        for i in 0..4 {
            kmer.add(i & 3);
        }
        assert_eq!(kmer.dna, 0xE4);
        assert_eq!(kmer.rc, 0xB1);
        assert_eq!(kmer.is_template(), false);
        assert_eq!(kmer.get_idx(true), 0x4E);
    }
    #[test]
    fn test_usize() {
        let mut kmer: Kmer<usize> = Kmer::new(32);
        for i in 0..32 {
            kmer.add(i & 3);
        }
        assert_eq!(kmer.dna, 0xE4E4E4E4E4E4E4E4);
        assert_eq!(kmer.rc, 0xB1B1B1B1B1B1B1B1);
        assert_eq!(kmer.is_template(), false);
        assert_eq!(kmer.get_idx(true), 0x4E4E4E4E4E4E4E4E);
    }
    #[test]
    fn unique() {
        let mut seen = vec![false; 256];
        let mut kmer: Kmer<u8> = Kmer::new(4);
        for i in (0..255).chain(once(255)) {
            for j in 0..4 {
                kmer.add((i >> (j << 1)) & 3);
            }
            let x = (if kmer.is_template() {1} else {0}) | kmer.get_idx(true) << 1;
            assert!( ! seen[x], format!("0x{:x} already seen!", x));
            seen[x] = true;
        }
        assert_eq!(vec![true; 256], seen);
    }
    #[test]
    fn test_revcmp() {
        let mut rng = thread_rng();
        let kmerlen = rng.gen_range(2, 32);
        let mut kmer: Kmer<u64> = Kmer::new(kmerlen);
        for _ in 0..32 {
            kmer.add(rng.gen_range(0, 4));
        }
        assert_eq!(kmer.dna.revcmp(kmerlen as usize), kmer.rc);
    }
    #[test]
    fn test_from_idx() {
        let mut rng = thread_rng();
        let kmerlen = rng.gen_range(12, 32);
        let mut test_dna = 0;
        let mut test_rc = 0;
        let mut test_ori = false;
        let mut test_idx = 0xffffffffffffffff;
        let last = rng.gen_range(kmerlen+1, 102);
        let pick = rng.gen_range(kmerlen, cmp::max(last-1, kmerlen+1));

        let mut kmer: Kmer<u64> = Kmer::new(kmerlen);
        for i in 0..last {
            kmer.add(rng.gen_range(0, 4));
            if i == pick {
                test_dna = kmer.dna;
                test_rc = kmer.rc;
                test_idx = kmer.get_idx(true);
                test_ori = kmer.is_template();
            }
        }
        assert_ne!(test_idx, 0xffffffffffffffff);
        let kmer2 = Kmer::from_idx(test_idx, kmerlen, test_ori);

        assert_eq!(test_dna, kmer2.dna);
        assert_eq!(test_rc, kmer2.rc);

        let kmer3 = Kmer::from_idx(test_idx, kmerlen, !test_ori);

        assert_eq!(test_dna, kmer3.rc);
        assert_eq!(test_rc, kmer3.dna);
    }
    #[test]
    fn extra() {
        let mut kmer: Kmer<u64> = Kmer::new(4);
        for _ in 0..16 {
            kmer.add(1);
        }
        assert_eq!(kmer.dna, 0x55);
        assert_eq!(kmer.rc, 0xff);
    }
}
