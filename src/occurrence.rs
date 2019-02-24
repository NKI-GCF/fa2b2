use std::collections::VecDeque;

use kmer::Kmer;
pub use kmerconst::{KmerConst,afstand};
use kmerloc::{KmerLoc,PriExtPosOri,MidPos};
use kmerstore::KmerStore;


// rename to Recurrance?
pub struct Occurrence<'a> {
	pub kc: &'a KmerConst,
	pub p: u64,
	pub kmer: Kmer<u64>,
	d: VecDeque<usize>,	   //
	pub mark: KmerLoc<u64>,
	pub i: usize,
	pub plim: u64,
}

impl<'a> Occurrence<'a> {
	pub fn new(ps: (u64, u64), kc: &'a KmerConst, ext: u64) -> Self {
		let kmp = ext | ps.0;

		Occurrence {
			kc,
			p: kmp,
			kmer: Kmer::new(kc.kmerlen as u32),
			d: VecDeque::from(vec![0; kc.max_no_kmers]),
			mark: KmerLoc::new(usize::max_value(), ext),
			i: 0,
			plim: ps.1,
		}
	}


	pub fn extend(&mut self) -> bool {
		if self.p.x() < self.kc.ext_max {
			self.p.extend();
			self.p.x() < self.kc.ext_max
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
		(hash ^ self.kc.ext_domain(x)) < self.mark.idx
	}

	fn next(&mut self, ori: bool) {

		self.i += 1;
		// ori == true if kmer is for template, then we want 1 in self.p
		self.p += if ori {3} else {2} - (1 & self.p);
	}
	fn mark_is_leaving(&self) -> bool {
		self.p.pos() - ((self.kc.max_no_kmers as u64) << 1) == self.mark.p.pos() + 2
	}

	fn set_next_mark(&mut self) {

		self.mark.reset();
		let x = self.p.x();
		let end_i = (self.kc.max_no_kmers - afstand(x, self.kc.kmerlen)) + 1;
		let ext = (1 << x) >> 1;
		debug_assert!(end_i != 0, "{:#x}, ext_max:{}", self.p, self.kc.ext_max);
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
		debug_assert!(self.mark.is_set());
	}

	/// add b2 to kmer, move .p according to occ orientation
	/// return whether required nr of kmers for struct occurance were seen, since contig start.
	/// adds 2bit to stored sequence
	pub fn complete(&mut self, ks: &KmerStore<u64>, b2: u8, n: usize) -> bool {

		let ori_change = self.kmer.update(b2, true);

		self.next(ori_change);
		if self.i >= self.kc.kmerlen { // some kmers
			let t = (self.i - self.kc.kmerlen) % self.kc.max_no_kmers;

			debug_assert!(t < self.d.len(), "i:{} kml:{} t:{}, dlen:{}", self.i, self.kc.kmerlen, t, self.d.len());
			self.d[t] = self.kmer.get_idx(true);
			//println!("[{}], idx:{}", self.i, self.d[t]);
			let x_start = if n == 0 {0} else {self.p.x()};

			for x in x_start..(self.p.x() + 1) {
				let t2 = self.kc.kmerlen + (if x == 0 {0} else {1 << x});
				if self.i < t2 {
					return false;
				}
				let hash = self.d[t] ^ if x == 0 {0} else {
					self.d[(self.i - t2) % self.kc.max_no_kmers]
				};
				let mut p = self.p.with_ext_prior(x);
				p -= (if x == 0 {0} else {1 << x}) as u64;
				if self.hash_is_extreme(hash, x) && ks.kmp[hash].is_replaceable_by(p) {
					self.mark.idx = hash;
					self.mark.p = p;
					break;
				}
			}

			if self.all_kmers() {
				if self.i > self.kc.readlen && dbgx!(self.mark_is_leaving()) {
					debug_assert!(n == 0);
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
	use super::KmerStore;
	use super::Occurrence;

	const READLEN: usize = 16;
	const SEQLEN: usize = 250;

	#[test]
	fn test_push_b2() {
		let c = KmerConst::new(READLEN, SEQLEN);
		let ks = KmerStore::<u64>::new(c.bitlen);
		let mut occ = Occurrence::new((0, 100), &c, 0);
		for _ in 0..6 {
			occ.complete(&ks, 1, 0);
		}
		assert_eq!(occ.kmer.dna, 0x55);
		assert_eq!(occ.kmer.rc, 0xff);
		assert_eq!(occ.p, 0xc);

		occ.complete(&ks, 2, 0);
		assert_eq!(occ.kmer.dna, 0x95);
		assert_eq!(occ.kmer.rc, 0xfc);
		assert_eq!(occ.p, 0xf);
		assert_eq!(occ.kmer.get_idx(true), 0x6a);
	}
	#[test]
	fn occurrence() {
		let c = KmerConst::new(READLEN, SEQLEN);
		let ks = KmerStore::<u64>::new(c.bitlen);
		let mut occ = Occurrence::new((0, 100), &c, 0);
		for _ in 0..8 {
			occ.complete(&ks, 0, 0);
		}
		for _ in 0..8 {
			occ.complete(&ks, 1, 0);
		}
		for _ in 0..8 {
			occ.complete(&ks, 2, 0);
		}
		for _ in 0..8 {
			occ.complete(&ks, 3, 0);
		}
	}
}

