use std::ops::BitOrAssign;
use std::ops::BitAnd;
use std::ops::BitAndAssign;
use num::FromPrimitive;
use num::ToPrimitive;
use num_traits::PrimInt;
use num_traits::NumAssign;
use std::fmt::LowerHex;

pub trait PriExtPosOri: PrimInt + BitAnd + NumAssign + BitOrAssign + BitAndAssign + FromPrimitive +
ToPrimitive + LowerHex {
	fn extension(&self) -> u64;
	fn x(&self) -> usize;
	fn extpos(&self) -> u64;
	fn same_ori(&self, p: u64) -> bool;
	fn top(&self) -> u64;
	fn pos(&self) -> u64;
	fn blacklist(&mut self);
	fn blacklisted(&self) -> bool;
	fn extend(&mut self);
	fn set_extension(&mut self, x: u64);
	fn priority(&self) -> u64;
	fn set_priority(&mut self);
	fn unset_priority(&mut self);
	fn clear_extension_and_priority(&mut self);
	fn is_same(&self, other: u64) -> bool;
}

impl PriExtPosOri for u64 {
	fn extension(&self) -> u64 {
		self & 0x7FFF000000000000
	}
	fn x(&self) -> usize {
		(self.extension() >> 48) as usize
	}
	fn extpos(&self) -> u64 {
		self & !(1 << 63)
	}
	fn same_ori(&self, p: u64) -> bool {
		(self & 1) == (p & 1)
	}
	fn top(&self) -> u64 {
		*self & !0xFFFFFFFFFFFF
	}
	fn pos(&self) -> u64 {
		*self & 0xFFFFFFFFFFFE
	}
	fn blacklist(&mut self) {
		if !self.blacklisted() {
			*self &= !0xFFFFFFFFFFFE;
			self.extend();
		}
	}
	fn blacklisted(&self) -> bool {
		self.pos() == 0
	}
	fn extend(&mut self) {
		*self += 1 << 48;
	}
	fn set_extension(&mut self, x: u64) {
		debug_assert!(x <= 0x7FFF);
		*self &= 0x8000FFFFFFFFFFFF;
		*self |= x << 48;
	}
	fn priority(&self) -> u64 {
		self & (1 << 63)
	}
	fn set_priority(&mut self) {
		*self |= 1 << 63;
	}
	fn unset_priority(&mut self) {
		*self &= 0x7FFFFFFFFFFFFFFF;
	}
	fn clear_extension_and_priority(&mut self) {
		*self &= 0xFFFFFFFFFFFF;
	}
	fn is_same(&self, other: u64) -> bool {
		*self == other
	}
}


/// for extended kmers, where position is midpoint
pub trait MidPos: PriExtPosOri {
	fn b2pos(&self) -> u64;
	fn is_replaceable_by(&self, new_entry: u64) -> bool;
	fn is_set_and_not(&self, other: u64) -> bool;
}

impl MidPos for u64 {
	/// The stored position is inbetween kmer or hash, and shifted. (pos << 1) - kmerlen + extension.
	/// this converts it to start of sequence. Zero-based so the first kmer is at 0.
	fn b2pos(&self) -> u64 {
		let p = *self & 0xFFFFFFFFFFFE;
		if p != 0 {
			let ext = 1 << self.x();
			debug_assert!(p >= ext, "{:x}, {:x}", p, ext);
			let half_ext = ext >> 1;
			if half_ext == 0 {
				p
			} else {
				p - half_ext
			}
		} else {
			0
		}
	}
	// *self == new_entry.top() ? blacklisted was extension below, but not new_entry's
	fn is_replaceable_by(&self, new_entry: u64) -> bool {
		*self <= new_entry.top()
	}
	fn is_set_and_not(&self, other: u64) -> bool {
		!self.blacklisted() && !self.is_same(other)
	}
}

#[derive(new, Clone, PartialEq)]
pub struct KmerLoc<T> {
	pub idx: usize,
	pub p: T,
}
impl<T: PriExtPosOri> KmerLoc<T> {
	pub fn reset(&mut self) {
		self.idx = usize::max_value();
		self.p = T::from_u64(0).unwrap();
	}

	pub fn is_set(&self) -> bool { self.idx != usize::max_value() }

	pub fn set(&mut self, idx: usize, p: T, x: usize) {
		self.idx = idx;
		self.p = p;
		self.p.set_extension(x as u64);
		self.p.set_priority();
	}

	pub fn priority(&self) -> u64 {
		self.p.to_u64().unwrap() & 0x8000000000000000
	}
	pub fn next(&mut self, ori: bool, is_template: bool) {
		// ori == true if kmer is for template, then we want 1 in self.p
		let p = self.p.to_u64().unwrap();
		if is_template {
			self.p += T::from_u64(if ori {3} else {2} - (1 & p)).unwrap();
		} else {
			self.p -= T::from_u64(if ori {1} else {2} + (1 & p)).unwrap();
		}
	}
}

#[cfg(test)]
mod tests {
    use super::KmerLoc;
    use rand::{thread_rng, Rng};
    #[test]
    fn forward() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, 100);
        let pick = rng.gen_range(20, 50);
        for _ in 0..pick {
            let ori = rng.gen_range(0, 2);
            kl.next(ori != 0, true);
            assert_eq!(ori, kl.p & 1);
        }
        assert_eq!(kl.p & !1, 100 + (pick << 1));
    }
    #[test]
    fn reverse() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, 100);
        let pick = rng.gen_range(20, 50);
        for _ in 0..pick {
            let ori = rng.gen_range(0, 2);
            kl.next(ori != 0, false);
            assert_eq!(ori, kl.p & 1);
        }
        assert_eq!(kl.p & !1, 100 - (pick << 1));
    }

}
