use std::collections::VecDeque;

use kmer::{Kmer,test_template};
pub use kmerconst::{KmerConst,afstand};
use kmerloc::{KmerLoc,PriExtPosOri};
use rdbg::STAT_DB;

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
		if self.p.x() + 1 < self.kc.ext_max {
			self.p.extend();
			true
		} else {false}
	}

	pub fn all_kmers(&self) -> bool {
		self.i >= self.kc.readlen
	}

	fn hash_is_extreme(&mut self, hash: usize, p: u64) -> bool {
		let xh = hash ^ self.kc.ext_domain(p.x());
		xh < self.mark.idx || (xh == self.mark.idx && p < self.mark.p)
	}

	fn complete_kmer(&mut self, b2: u8) -> bool {

		// ori == true if kmer is for template, then we want 1 in self.p
		let new_idx = match self.i.checked_sub(self.kc.kmerlen).map(|i| i % self.kc.no_kmers) {
			Some(old_idx) => {
				let new_idx = if old_idx + 1 != self.kc.no_kmers {old_idx + 1} else {0};
				self.d[new_idx] = self.d[old_idx];
				new_idx
			},
			None => 0
		};
		dbg_assert!(new_idx < self.d.len(), "new_idx >= self.d.len(), {} >= {}", new_idx, self.d.len());
		self.p += if self.d[new_idx].update(b2, true) {3} else {2} - (1 & self.p);
		self.i += 1;
		self.i >= self.kc.kmerlen
	}
	fn mark_is_leaving(&self, new_xmer: bool) -> bool {

		dbgf!(self.p.pos(), "{} - ({} << 1) == {} ?", self.kc.no_kmers, self.mark.p.pos());
		self.p.pos() - ((self.kc.no_kmers as u64) << 1) == self.mark.p.pos() || {

			dbg_assert!(new_xmer || self.kc.readlen != self.kc.kmerlen);
			false
		}
	}

	pub fn set_next_mark(&mut self) -> bool {

		self.mark.reset();
		let x = self.p.x();
		let afs = afstand(x, self.kc.kmerlen);
		if dbgx!(self.kc.no_kmers <= afs) {
			return false;
		}
		for i in 0..(self.kc.no_kmers - afs) {
			let _ = self.set_if_extreme(i, x);
		}
		dbg_assert!(self.mark.is_set());
		true
	}

	fn set_if_extreme(&mut self, i: usize, x: usize) -> bool
	{
		let base = self.i - self.kc.kmerlen;
		let afs = afstand(x, self.kc.kmerlen);
		let d_i = base.wrapping_sub(i) % self.kc.no_kmers;
		let mut kmer = self.d[d_i];
		if x > 0 {
			let d_i2 = base.wrapping_sub(afs + i) % self.kc.no_kmers;
			dbg_print!("dna:{:#x}, dna2:{:#x}", self.d[d_i2].dna, kmer.dna);
			kmer.hash(self.d[d_i2]);
		}
		let hash = kmer.get_idx(true);
		let mut p = self.p.with_ext_prior(x) - ((afs + i) << 1) as u64;
		if test_template(kmer.dna, kmer.rc) {
			p |= 1;
		}
		if self.hash_is_extreme(hash, p) {
			self.mark.set(hash, p, x);
			true
		} else {
			false
		}
	}

	/// add b2 to kmer, move .p according to occ orientation
	/// return whether required nr of kmers for struct occurance were seen, since contig start.
	/// adds 2bit to stored sequence
	pub fn complete(&mut self, b2: u8, x_start: usize) -> bool {

		if self.complete_kmer(b2) {
			let mut new_xmer = false;

			for x in x_start..(self.p.x() + 1) {
				let afs = afstand(x, self.kc.kmerlen);
				if self.i < self.kc.kmerlen + afs {
					return false;
				}
				if self.set_if_extreme(0, x) {
					dbgf!(self.mark.p, "{:#x}");
					new_xmer = true;
					break;
				}
			}
			dbg_assert!(self.mark.is_set(), "or return false?");

			if self.all_kmers() {
				if self.i > self.kc.readlen && dbgx!(self.mark_is_leaving(new_xmer)) {
					// there is a leaving kmer and the extreme was popped.
					// need to reestablish new ext from all kmers and extensions.
					// but cannot be less than this extension:
					return self.set_next_mark();
				}
				return true
			}
		}
		false
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	const READLEN: usize = 16;
	const SEQLEN: usize = 250;

	#[test]
	fn test_push_b2() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut occ = Occurrence::new((0, 100), &kc, 0);
		for _ in 0..6 {
			occ.complete(1, 0);
		}
		let mut kmer = occ.kmer();
		dbg_assert_eq!(kmer.dna, 0x55);
		dbg_assert_eq!(kmer.rc, 0xff);
		dbg_assert_eq!(occ.p, 0xc);

		occ.complete(2, 0);
		kmer = occ.kmer();
		dbg_assert_eq!(kmer.dna, 0x95);
		dbg_assert_eq!(kmer.rc, 0xfc);
		dbg_assert_eq!(occ.p, 0xf);
		dbg_assert_eq!(kmer.get_idx(true), 0x6a);
	}
	#[test]
	fn occurrence() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut occ = Occurrence::new((0, 100), &kc, 0);
		for _ in 0..8 {
			occ.complete(0, 0);
		}
		for _ in 0..8 {
			occ.complete(1, 0);
		}
		for _ in 0..8 {
			occ.complete(2, 0);
		}
		for _ in 0..8 {
			occ.complete(3, 0);
		}
	}
}

