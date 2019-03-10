use std::collections::VecDeque;

use kmer::{Kmer,test_template};
pub use kmerconst::{KmerConst,afstand};
use kmerloc::{KmerLoc,PriExtPosOri,MidPos};
use kmerstore::KmerStore;


// rename to Recurrance?
pub struct Occurrence<'a> {
	pub kc: &'a KmerConst,
	pub p: u64,
	d: VecDeque<Kmer<u64>>,	   //
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
			d: VecDeque::from(vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers]),
			mark: KmerLoc::new(usize::max_value(), ext),
			i: 0,
			plim: ps.1,
		}
	}
	pub fn kmer(&self) -> Kmer<u64> {
		self.d[(self.i - self.kc.kmerlen) % self.kc.no_kmers]
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

	fn complete_kmer(&mut self, b2: u8) -> bool{

		// ori == true if kmer is for template, then we want 1 in self.p
		let new_idx = match self.i.checked_sub(self.kc.kmerlen).map(|i| i % self.kc.no_kmers) {
			Some(old_idx) => {
				let new_idx = if old_idx + 1 != self.kc.no_kmers {old_idx + 1} else {0};
				self.d[new_idx] = self.d[old_idx];
				new_idx
			},
			None => 0
		};
		self.p += if self.d[new_idx].update(b2, true) {3} else {2} - (1 & self.p);
		self.i += 1;
		self.i >= self.kc.kmerlen
	}
	fn mark_is_leaving(&self) -> bool {
		self.p.pos() - ((self.kc.no_kmers as u64) << 1) == self.mark.p.pos()
	}

	fn set_next_mark(&mut self) {

		self.mark.reset();
		let x = self.p.x();
		let offs = afstand(x, self.kc.kmerlen);
		let end_i = self.kc.no_kmers - offs;

		debug_assert!(end_i != 0, "{:#x}, ext_max:{}", self.p, self.kc.ext_max);
		let base = self.i - self.kc.kmerlen;

		for i in 0..end_i {

			let d_i = base.wrapping_sub(i) % self.kc.no_kmers;

			let mut kmer = self.d[d_i];
			if x > 0 {
				let d_i2 = base.wrapping_sub(i+offs) % self.kc.no_kmers;
				kmer.hash(self.d[d_i2]);
			}
			let idx = kmer.get_idx(true);
			if self.hash_is_extreme(idx, x) {
				let p = self.p - ((end_i + offs - i) << 1) as u64;
				self.mark.set(idx, p, x);
			}
		}
		debug_assert!(self.mark.is_set());
	}

	/// add b2 to kmer, move .p according to occ orientation
	/// return whether required nr of kmers for struct occurance were seen, since contig start.
	/// adds 2bit to stored sequence
	pub fn complete(&mut self, ks: &KmerStore<u64>, b2: u8, n: usize) -> bool {

		if self.complete_kmer(b2) {

			//println!("[{}], idx:{}", self.i, self.d[t]);
			let x_start = if n == 0 {0} else {self.p.x()};

			for x in x_start..(self.p.x() + 1) {
				let afs = afstand(x, self.kc.kmerlen);
				let t2 = self.kc.kmerlen + afs;
				if self.i < t2 {
					return false;
				}
				let mut kmer = self.kmer();
				if x != 0 {
					kmer.hash(self.d[(self.i - afs) % self.kc.no_kmers]);
				}
				let hash = kmer.get_idx(true);
				let mut p = self.p.with_ext_prior(x) - (afs << 1) as u64;
				if self.hash_is_extreme(hash, x) && dbgx!(ks.kmp[hash].is_replaceable_by(p) || ks.kmp[hash].is_same(p)) {
					self.mark.idx = hash;
					if test_template(kmer.dna, kmer.rc) {
						p |= 1;
					}
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
		let kc = KmerConst::new(READLEN, SEQLEN);
		let ks = KmerStore::<u64>::new(kc.bitlen);
		let mut occ = Occurrence::new((0, 100), &kc, 0);
		for _ in 0..6 {
			occ.complete(&ks, 1, 0);
		}
		let mut kmer = occ.kmer();
		assert_eq!(kmer.dna, 0x55);
		assert_eq!(kmer.rc, 0xff);
		assert_eq!(occ.p, 0xc);

		occ.complete(&ks, 2, 0);
		kmer = occ.kmer();
		assert_eq!(kmer.dna, 0x95);
		assert_eq!(kmer.rc, 0xfc);
		assert_eq!(occ.p, 0xf);
		assert_eq!(kmer.get_idx(true), 0x6a);
	}
	#[test]
	fn occurrence() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let ks = KmerStore::<u64>::new(kc.bitlen);
		let mut occ = Occurrence::new((0, 100), &kc, 0);
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

