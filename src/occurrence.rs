use std::collections::VecDeque;

use kmer::Kmer;
pub use kmerconst::KmerConst;
use kmerloc::{KmerLoc,PriExtPosOri,MidPos};
use kmerstore::KmerStore;


// rename to Recurrance?
pub struct Occurrence<'a> {
	pub kc: &'a KmerConst,
	pub p: u64,
	pub is_template: bool, //
	pub kmer: Kmer<u64>,
	d: VecDeque<usize>,	   //
	pub mark: KmerLoc<u64>,
	pub i: usize,
	pub plim: u64,
}

impl<'a> Occurrence<'a> {
	pub fn new(ps: (u64, u64), kc: &'a KmerConst, ext: u64, is_template: bool) -> Self {
		let kmp = ext | if is_template { ps.0 } else { ps.1 };

		Occurrence {
			kc,
			p: kmp,
			is_template,
			kmer: Kmer::new(kc.kmerlen as u32),
			d: VecDeque::from(vec![0; kc.max_no_kmers]),
			mark: KmerLoc::new(usize::max_value(), ext),
			i: 0,
			plim: if is_template { ps.1 } else { ps.0 },
		}
	}


	pub fn extend(&mut self) -> bool {
		if self.p.x() < self.kc.ext_max() {
			self.p.extend();
			self.p.x() < self.kc.ext_max()
		} else {false}
	}

	pub fn try_extension_redefine_minimum(&mut self) -> bool {
		if self.extend() {
			self.set_next_mark();
			true
		} else {false}
	}

	pub fn all_kmers(&self) -> bool {
		self.i >= self.kc.readlen
	}

	pub fn hash_is_extreme(&mut self, hash: usize, x: usize) -> bool {
		!self.mark.is_set() || if self.is_template {
			(hash ^ self.kc.ext_domain(x)) < self.mark.idx
		} else {
			(hash ^ self.kc.ext_domain(x)) <= self.mark.idx
		}
	}

	pub fn get_ori_for_stored(&self, stored_at_index: u64) -> bool {
		if (stored_at_index & 1) == (self.mark.p & 1) {
			self.is_template
		} else {
			!self.is_template
		}
	}

	fn next(&mut self, ori: bool, is_template: bool) {

		self.i += 1;
		// ori == true if kmer is for template, then we want 1 in self.p
		if is_template {
			self.p += if ori {3} else {2} - (1 & self.p);
		} else {
			self.p -= if ori {1} else {2} + (1 & self.p);
		}
	}
	fn mark_is_leaving(&self) -> bool {
		self.p.pos() - ((self.kc.max_no_kmers as u64) << 1) == self.mark.p.pos() + 2
	}

	fn set_next_mark(&mut self) {

		self.mark.reset();
		let x = self.p.x();
		assert!(x < self.kc.ext_max());
		let ext = (1 << x) >> 1;
		let end_i = self.kc.max_no_kmers - (ext << 1);
		let base = self.i - self.kc.kmerlen;

		for i in 0..end_i {

			let d_i = (base + i - end_i) % self.kc.max_no_kmers;

			let mut hash = self.d[d_i];
			if x > 0 {
				let d_i2 = (base + i - self.kc.max_no_kmers) % self.kc.max_no_kmers;
				hash ^= self.d[d_i2];
			}
			if self.hash_is_extreme(hash, x) {
				let p = self.p - ((end_i + ext - i) << 1) as u64;
				self.mark.set(hash, p, x);
			}
		}
	}

	/// add b2 to kmer, move .p according to occ orientation
	/// return whether required nr of kmers for struct occurance were seen, since contig start.
	/// adds 2bit to stored sequence
	pub fn complete(&mut self, b2: u8) -> bool {

		let is_template = self.is_template;
		let ori_change = self.kmer.update(b2, is_template);

		self.next(ori_change, is_template);
		if self.i >= self.kc.kmerlen { // some kmers
			let t = (self.i - self.kc.kmerlen) % self.kc.max_no_kmers;

			self.d[t] = self.kmer.get_idx(is_template);
			//println!("[{}], idx:{}", self.i, self.d[t]);
			let x = self.p.x();
			let mut hash = self.d[t];
			if x > 0 {
				let t2 = self.kc.kmerlen + (1 << x);
				if self.i < t2 {
					return false;
				}
				hash ^= self.d[(self.i - t2) % self.kc.max_no_kmers];
			}
			if self.hash_is_extreme(hash, x) {
				let p = self.p - (if x == 0 {0} else {1 << x}) as u64;
				self.mark.set(hash, dbgf!(p, "{:#016x}"), x);
			}

			if self.all_kmers() {
				if self.i > self.kc.readlen && dbgx!(self.mark_is_leaving()) {
					// there is a leaving kmer and the extreme was popped.
					// need to reestablish new ext from all kmers and extensions.
					// but cannot be less than this extension:
					self.set_next_mark();
				}
				return true
			}
		}
		false
	}
}

#[cfg(test)]
mod tests {
	use super::KmerConst;
	use super::Occurrence;

	const READLEN: usize = 16;
	const SEQLEN: usize = 250;
	const EXTENT: u32 = 48;

	#[test]
	fn test_push_b2() {
		let c = KmerConst::new(READLEN, SEQLEN, EXTENT);
                let is_template = true;
		let mut occ = Occurrence::new((0, 100), &c, 0, is_template);
		for _ in 0..6 {
			occ.complete(1);
		}
		assert_eq!(occ.kmer.dna, 0x55);
		assert_eq!(occ.kmer.rc, 0xff);
		assert_eq!(occ.p, 0xc);

		occ.complete(2);
		assert_eq!(occ.kmer.dna, 0x95);
		assert_eq!(occ.kmer.rc, 0xfc);
		assert_eq!(occ.p, 0xf);
		assert_eq!(occ.kmer.get_idx(true), 0x6a);
	}
	#[test]
	fn test_push_b2_rc() {
		let c = KmerConst::new(READLEN, SEQLEN, EXTENT);
                let is_template = false;
		let mut occ = Occurrence::new((0, 16), &c, 0, is_template);
		occ.complete(2);
		for _ in 0..3 {
			occ.complete(1);
		}
		assert_eq!(occ.p, 0x9);
		assert_eq!(occ.kmer.dna, 0xfc);
		assert_eq!(occ.kmer.rc, 0x95);
		assert_eq!(occ.kmer.get_idx(true), 0x6a);
	}
	#[test]
	fn test_ori1() {
		let c = KmerConst::new(READLEN, SEQLEN, EXTENT);
                let is_template = true;
		let mut occ = Occurrence::new((0, 100), &c, 0, is_template);
		occ.kmer.dna = 0xaa;
		assert_eq!(occ.kmer.is_template(), true);
	}

	#[test]
	fn test_ori2() {
		let c = KmerConst::new(READLEN, SEQLEN, EXTENT);
                let is_template = true;
		let mut occ = Occurrence::new((0, 100), &c, 0, is_template);
		occ.kmer.rc = 0xaa;
		assert_eq!(occ.kmer.is_template(), false);
	}

	#[test]
	fn occurrence() {
		let c = KmerConst::new(READLEN, SEQLEN, EXTENT);
                let is_template = true;
		let mut occ = Occurrence::new((0, 100), &c, 0, is_template);
		for _ in 0..8 {
			occ.complete(0);
		}
		for _ in 0..8 {
			occ.complete(1);
		}
		for _ in 0..8 {
			occ.complete(2);
		}
		for _ in 0..8 {
			occ.complete(3);
		}
	}
}

