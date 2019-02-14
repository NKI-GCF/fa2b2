extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate num_traits;
extern crate num;

use std::cmp;
use std::slice::Iter;

use kmerconst::KmerConst;
use occurrence::Occurrence;
use kmerloc::{PriExtPosOri,MidPos};
use kmerstore::KmerStore;

pub struct KmerIter<'a> {
	n_stretch: u64,
	goffs: u64,
	occ: Vec<Occurrence<'a>>,
	ks: &'a mut KmerStore<u64>,
} //^-^\\

impl<'a> KmerIter<'a> {
	pub fn new(ks: &'a mut KmerStore<u64>, kc: &'a KmerConst) -> Self {
		let occ = vec![Occurrence::new((0, u64::max_value()), kc, 0, true)];
		KmerIter { n_stretch: 0, goffs: 0, occ, ks }
	}

	fn finalize_n_stretch_if_pending(&mut self) {
		if self.n_stretch > 0 {
			eprintln!("added new contig. Ns:{}", self.n_stretch);
			self.ks.offset_contig(self.n_stretch);
			self.goffs += self.n_stretch;
			self.n_stretch = 0;
		}
	}

	/// if N, insert contig (once) and count stretch, false
	/// else store and update occurrance kmers, complete ? true : false
	fn complete_occurance_or_contig(
		&mut self,
		n: usize,
		b2: u8,
	) -> bool {
		if let Some(occ) = self.occ.get_mut(n) {
			let p = occ.p.pos();
			if b2 < 4 {
				if let Some(qb) = self.ks.b2.get_mut(p as usize >> 3) {
					*qb |= b2 << (p & 6);
				}
				self.ks.p_max = (p + 2) & !1;
				if occ.complete(b2) {
					return true;
				}
			} else if occ.i != 0 {
				self.goffs += occ.i as u64;
				eprintln!("started N-stretch at {}.", p);
				self.ks.push_contig(p, self.goffs);

				// clear all except orientation and position to rebuild at the start of a new contig.
				//assert_eq!(occ.i, 1);
				occ.d.clear();
				occ.i = 0;
				self.n_stretch = 1;
				return false;
			} else {
				self.n_stretch += 1;
				return false;
			}
		}
		self.finalize_n_stretch_if_pending();
		false
	}

	/// when rebuilding eq for recurrent kmer, and extending take into account contig boundaries
	/// for that site
	fn get_plimits(&self, p: u64) -> (u64, u64) {

		// binary search; limit endp to end of contig
		debug_assert!(p != 0);
		let i = self.ks.get_contig(p);

		debug_assert!(i < self.ks.contig.len());
		debug_assert!(self.ks.contig[i].twobit <= p);

		let x = self.ks.get_twobit_before(i).unwrap_or(0);
		let kc = self.occ[0].kc;
		let addl = ((kc.max_no_kmers + kc.kmerlen) << 1) as u64;
		let p_max = self.ks.p_max;
		//eprintln!("p:{:x} > x:{:x} + addl:{:x} ?",  p, x, addl);
		let left = if p > x + addl { p - addl } else { x };
		let right = cmp::min(
			p + (kc.readlen << 1) as u64,
			self.ks.get_twobit_after(i).unwrap_or(p_max),
		);
		(left, right)
	}

	fn rebuild_kmer_stack_with_extension(&self, stored_at_index: &mut u64, is_template: bool) -> Occurrence<'a> {

		stored_at_index.extend();

		let plimits = self.get_plimits(stored_at_index.b2pos());

		Occurrence::new(plimits, self.occ[0].kc, stored_at_index.extension(), is_template)
	}

	fn search_occ_for_pos(&self, original_n: usize, stored_at_index: u64) -> usize {

		assert!(original_n == self.occ.len() - 1);

		for n in 0..original_n {
			if dbgf!(stored_at_index == self.occ[n].mark.p,
					"{:#?}: {:#016x} == {:#016x}?", stored_at_index, self.occ[n].mark.p) {
				return dbgx!(n);
			}
		}
		original_n
	}

	/// when a min_pos can replace an existing entry in ks.kmp, the existing entry
	/// needs extension, and to be added to occ.
	fn next_after_replacement(&mut self, stored_at_index: &mut u64, original_n: usize) -> usize {

		let original_stored_ori = self.occ[original_n].get_ori_for_stored(*stored_at_index);

		let mut n = self.search_occ_for_pos(original_n, *stored_at_index);
		if self.occ[dbgx!(n)].mark.p != *stored_at_index {

			let next = self.rebuild_kmer_stack_with_extension(stored_at_index, dbgx!(original_stored_ori));

			if n != 0 {
				// position is written below, this self.occ[n] is done.
				// Could add it to recurring kmer_stacks, though.
				self.occ[n] = next;
			} else {
				// 0th is never overwritten.
				self.occ.push(next);
				n += 1;
			}
		}
		n
	}

	fn next_b2(&self, seq: &mut Iter<u8>, n: usize) -> Option<u8> {
		if n == 0 {
			seq.next().map(|c| (c >> 1) & 7u8)
		} else {
			//dbgf!(self.occ[n].mark.p, "{:#016x}");
			self.ks.b2_for_p(self.occ[n].p.pos())
		}
	}

	pub fn markcontig<T: MidPos>(
		&mut self,
		seq: &mut Iter<u8>,
	) -> u64 {

		let mut n = 0;
		self.goffs = 0;
		self.ks.push_contig(self.occ[n].p.pos(), self.goffs);

		while let Some(b2) = self.next_b2(seq, n) {

			if !self.complete_occurance_or_contig(n, b2) {
				continue;
			}
			if n != 0 { dbgf!(self.occ[n].mark.p, "{:#016x}"); }
			loop {

				let (min_index, min_pos) = (self.occ[n].mark.idx, self.occ[n].mark.p);

				let mut stored_at_index = self.ks.kmp[min_index];

				if dbgf!(stored_at_index.is_replaceable_by(min_pos),
					"{:#?}\nstored_at_index:{:#016x}, min_pos:{:#016x}\n", stored_at_index, min_pos) {

					self.ks.kmp[min_index] = min_pos;
					dbgf!(self.ks.kmp[min_index], "{:#x} [{:#x}]", min_index);

					if dbgx!(stored_at_index.is_set_and_not(min_pos)) {

						n = self.next_after_replacement(&mut stored_at_index, n);

					} else if dbgx!(n > 0) {
						let recurrent = self.occ.pop();
						// TODO add to another stack for fast lookup for multimappers.
						n -= 1;
						continue;
					}

					break; // position written (done) or added next_stack which requires extension.
				}
				if dbgx!(self.occ[n].try_extension_redefine_minimum()) {
					dbgf!(self.occ[n].mark.p, "{:x}, n:{}", n);

					if dbgx!(stored_at_index.extension() == min_pos.extension()) {
						self.ks.kmp[min_index].blacklist();
						// both stored and current require extension

						// keep before asignment of n:
						let stored_ori = self.occ[n].get_ori_for_stored(stored_at_index);

						// search whether stored is already in self.occ (TODO enable binary search?)
						n = self.search_occ_for_pos(n, stored_at_index);

						if dbgf!(stored_at_index != self.occ[n].mark.p,
								 "{:#?}\nn:{}: {:x}, {:x}", n, stored_at_index, self.occ[n].mark.p) {

							let next_stack = self.rebuild_kmer_stack_with_extension(&mut stored_at_index, dbgx!(stored_ori));

							n = self.occ.len();
							self.occ.push(next_stack);
							// has to be rebuilt
						}
						break; // next_stack .
					}
				} else {
					self.ks.kmp[min_index].blacklist();
					if dbgx!(n == 0) { // never pop 0th. 0th needs to be renewed when done.
						break;
					}
					let blacklisted = self.occ.pop();
					// TODO: add blacklisted in unmappable regions.
					n -= 1;
					if dbgx!(n == 0 && self.occ[n].mark.p == self.ks.kmp[self.occ[n].mark.idx]) {
						// test already done
						break;
					}
				}
			}
		}
		self.finalize_n_stretch_if_pending();
		self.occ[0].p.pos() as u64
	}
}

#[cfg(test)]
mod tests {
	use super::KmerConst;
	use super::KmerIter;
	use super::KmerStore;
	const READLEN: usize = 16;
	const SEQLEN: usize = 250;
	const EXTENT: u32 = 48;

	#[test]
	fn test_16n() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new(&mut ks, &kc);
			let seq1: Vec<u8> = b"NNNNNNNNNNNNNNNN"[..].to_owned();
			kmi.markcontig::<u64>(&mut seq1.iter());
		}
		assert_eq!(ks.contig.len(), 1);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 16);
	}
	#[test]
	fn test_1n() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new(&mut ks, &kc);
			let seq1: Vec<u8> = b"N"[..].to_owned();
			kmi.markcontig::<u64>(&mut seq1.iter());
		}
		assert_eq!(ks.contig.len(), 1);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 1);
	}
	#[test]
	fn test_1n1c1n() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new(&mut ks, &kc);
			let seq1: Vec<u8> = b"NCN"[..].to_owned();
			kmi.markcontig::<u64>(&mut seq1.iter());
		}
		assert_eq!(ks.contig.len(), 2);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 1);
		assert_eq!(ks.contig[1].twobit, 2);
		assert_eq!(ks.contig[1].genomic, 3);

	}
	#[test]
	fn test_17_c() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		assert_eq!(ks.kmp.len(), 128);
		{
			let mut kmi = KmerIter::new(&mut ks, &kc);

			// XXX ik verwacht hier (17 x C!!) een multimapper, geen positie, maar geblacklist.
			let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCCC"[..].to_owned();
			kmi.markcontig::<u64>(&mut seq1.iter());
		}
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		let unresolveable = (1 << 63) | ((kc.ext_max() as u64) << 48);
		assert_eq!(ks.kmp[0], unresolveable, "test_17_c(): {:x}", ks.kmp[0]);
	}
	#[test]
	fn test_1n18c1n() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new(&mut ks, &kc);

			// XXX ik verwacht hier (17 x C!!) een multimapper, geen positie, maar geblacklist.
			let seq1: Vec<u8> = b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned();
			kmi.markcontig::<u64>(&mut seq1.iter());
		}
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		let unresolveable = (1 << 63) | ((kc.ext_max() as u64) << 48);
		assert_eq!(ks.kmp[0], unresolveable);
	}
	#[test]
	fn test_1n16c() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new(&mut ks, &kc);

			// XXX ik verwacht hier (17 x C!!) een multimapper, geen positie, maar geblacklist.
			let seq1: Vec<u8> = b"NCCCCCCCCCCCCCCCC"[..].to_owned();
			kmi.markcontig::<u64>(&mut seq1.iter());
		}
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		let first_pos = (1 << 63) | ((kc.kmerlen as u64)<< 1);
		assert_eq!(ks.kmp[0], first_pos);
	}
	#[test]
	fn test_16at1a() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new(&mut ks, &kc);

			// XXX ik verwacht hier (17 x C!!) een multimapper, geen positie, maar geblacklist.
			let seq1: Vec<u8> = b"ATATATATATATATATA"[..].to_owned();
			kmi.markcontig::<u64>(&mut seq1.iter());
		}
		for i in 0..ks.kmp.len() {
			if i != 0x0 && i != 0x22 {
				assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
			}
		}
		let unresolveable = (1 << 63) | ((kc.ext_max() as u64) << 48);
		assert_eq!(ks.kmp[0], unresolveable);
	}
}
