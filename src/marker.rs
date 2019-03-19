extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate num_traits;
extern crate num;

use std::cmp;
use std::slice::Iter;

use occurrence::Occurrence;
use kmerloc::{PriExtPosOri,MidPos};
use kmerstore::KmerStore;
use rdbg::STAT_DB;

pub struct KmerIter<'a> {
	n_stretch: u64,
	goffs: u64,
	pub(super) occ: &'a mut Vec<Occurrence<'a>>,
	pub(super) ks: &'a mut KmerStore<u64>,
} //^-^\\

impl<'a> KmerIter<'a> {
	pub fn new(ks: &'a mut KmerStore<u64>, occ: &'a mut Vec<Occurrence<'a>>) -> Self {
		KmerIter { n_stretch: 0, goffs: 0, occ, ks }
	}

	fn finalize_n_stretch_if_pending(&mut self) {
		if self.n_stretch > 0 {
			dbg_print!("added new contig. Ns:{}", self.n_stretch);
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
				if occ.complete(b2, 0) {
					return true;
				}
			} else if occ.i != 0 {
				self.goffs += occ.i as u64;
				dbg_print!("started N-stretch at {}.", p);
				self.ks.push_contig(p, self.goffs);

				// clear all except orientation and position to rebuild at the start of a new contig.
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
		dbg_assert!(p != 0);
		let i = self.ks.get_contig(p);

		dbg_assert!(i < self.ks.contig.len());
		dbg_assert!(self.ks.contig[i].twobit <= p);

		let contig_start = self.ks.get_twobit_before(i).unwrap_or(0);
		let kc = self.occ[0].kc;
		let p_max = self.ks.p_max;

		let p_rl = (kc.readlen << 1) as u64;

		let left = if p >= contig_start + p_rl {
			p - p_rl
		} else {
			contig_start
		};
		dbgf!(left, "{:?}, p:{}, contig_start:{}, p_rl:{}", p, contig_start, p_rl);

		let right = cmp::min(
			p + ((kc.readlen - kc.kmerlen) << 1) as u64,
			self.ks.get_twobit_after(i).unwrap_or(p_max),
		);
		dbgf!(right, "{:?}, p_max:{}", self.ks.get_twobit_after(i).unwrap_or(p_max));

		(left, right)
	}

	fn rebuild_kmer_stack(&mut self, min_index: usize, new_val: u64) -> Occurrence<'a> {

		let mut stored = self.ks.kmp[min_index];

		let plimits = self.get_plimits(stored.pos());
		let mut occ = Occurrence::new(plimits, self.occ[0].kc, stored.extension());

		while {
			let p = occ.p.pos();
			let b2 = self.ks.b2_for_p(p).unwrap();

			dbg_print!("=> b2 {:x} <=", b2);
			let x = occ.p.x();
			let _ = occ.complete(b2, x);
			if occ.p.pos() <= stored.pos() {
				if occ.p.pos() == stored.pos() {
					self.ks.kmp[min_index] = new_val;
				}
				true
			} else {
				if occ.p.pos() == stored.pos() + 2 {
					stored = occ.mark.p;
				}
				occ.mark.p == stored && occ.p.pos() < occ.plim
			}
		} {}

		occ
	}

	fn search_occ_for_pos(&self, original_n: usize, stored_at_index: u64) -> usize {


		for n in 0..original_n {
			if dbgf!(stored_at_index == self.occ[n].mark.p,
					"{:#?}: {:#016x} == {:#016x}?", stored_at_index, self.occ[n].mark.p) {
				return dbgx!(n);
			}
		}
		original_n
	}
	fn add_newstack(&mut self, next_stack: Occurrence<'a>, n: usize) {
		if n < self.occ.len() {
			self.occ[n] = next_stack;
		} else {
			self.occ.push(next_stack);
		}
	}

	/// when the stored entry in ks.kmp is replaced or blacklisted, the occurence for the stored
	/// needs to be extended to find, there, a minpos for the extended kmers.
	/// the entry is added to the stack, return the index.
	fn get_next_for_extension(&mut self, min_index: usize, new_val: u64, original_n: usize,
						   is_replaced: bool) -> usize {

		let stored_at_index = self.ks.kmp[min_index];
		let mut n = self.search_occ_for_pos(original_n, stored_at_index);
		if dbgx!(stored_at_index != self.occ[n].mark.p) {

			let next_stack = self.rebuild_kmer_stack(min_index, new_val);

			if dbgx!(is_replaced && n != 0) {
				// position is written below, this self.occ[n] is done.
				// Could add it to recurring kmer_stacks, though.
				if next_stack.p.pos() < next_stack.plim {
					self.occ[n] = next_stack;
				} else if n != 0 {
					n -= 1;
				}
			} else if next_stack.p.pos() < next_stack.plim {
				// 0th is never overwritten.
				n += 1;
				self.add_newstack(next_stack, n);
			}
		}
		n
	}

	fn next_b2(&self, seq: &mut Iter<u8>, n: usize) -> Option<u8> {
		if n == 0 {
			seq.next().map(|c| (c >> 1) & 7u8)
		} else {
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
			dbg_print!("=> twobit {:x} (n: {}) <=", b2, n);

			if !self.complete_occurance_or_contig(n, b2) {
				continue;
			}
			loop {
				dbg_print!("---------[ n: {} ]-------------", n);

				let (min_index, min_pos) = (self.occ[n].mark.idx, self.occ[n].mark.p);
				dbg_assert!(min_index < self.ks.kmp.len(), "{:x}, {:x}", min_pos, self.occ[n].p);
				let mut stored_at_index = self.ks.kmp[min_index];

				if dbgx!(stored_at_index.is_replaceable_by(min_pos)) {
					dbg_print!("[{:#x}] (={:#016x}) <= {:#016x}", min_index, stored_at_index, min_pos);

					self.ks.kmp[min_index] = min_pos;

					if dbgx!(stored_at_index.is_set_and_not(min_pos)) {

						n = self.get_next_for_extension(min_index, min_pos, n, true);

					} else {
						self.ks.kmp[min_index] = min_pos;
						if n > 0 {
							// We're done on a recurrence when observing an already set value, not mark.p
							if dbgf!(stored_at_index.is_same(min_pos) || self.occ[n].p.pos() == self.occ[n].plim,
									"{:?}, {:#x}, {:#x}", self.occ[n].p, self.occ[n].plim) {
								n -= 1;
								continue;
							}
						}
					}

					break; // position written (done) or added next_stack which requires extension.
				}
				if dbgx!(self.occ[n].extend() && self.occ[n].set_next_mark()) {
					dbgf!(self.occ[n].mark.p, "{:x}");

					if dbgx!(stored_at_index.extension() == min_pos.extension()) {
						let mut blacklist = self.ks.kmp[min_index];
						blacklist.blacklist();
						dbgf!(min_index, "[{:#x}] is blacklist'd ({:#x})", self.ks.kmp[min_index]);

						// both stored and current require extension. Stored is handled now.
						n = self.get_next_for_extension(min_index, blacklist, n, false);

						break;
					}
				} else {
					self.ks.kmp[min_index].blacklist();
					dbgf!(min_index, "[{:#x}] is blacklisted ({:#x})", self.ks.kmp[min_index]);
					if dbgx!(n == 0) { // 0th needs to be renewed when done.
						break;
					}
					n -= 1;
					if dbgx!(n == 0 && self.occ[n].mark.p == self.ks.kmp[self.occ[n].mark.idx]) {
						// test already done
						break;
					}
				}
			}
		}
		dbgx!(self.finalize_n_stretch_if_pending());
		self.occ[0].p.pos() as u64
	}
}

#[cfg(test)]
mod tests {
	use kmerconst::KmerConst;
	use super::*;
	const READLEN: usize = 16;
	const SEQLEN: usize = 250;

	fn process<'a>(occ: &'a mut Vec<Occurrence<'a>>, ks: &'a mut KmerStore<u64>, seq: Vec<u8>) {
		let mut kmi = KmerIter::new(ks, occ);
		kmi.markcontig::<u64>(&mut seq.iter());
	}

	#[test]
	fn test_16n() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"NNNNNNNNNNNNNNNN"[..].to_owned());
		}
		dbg_assert_eq!(ks.contig.len(), 1);
		dbg_assert_eq!(ks.contig[0].twobit, 0);
		dbg_assert_eq!(ks.contig[0].genomic, 16);
	}
	#[test]
	fn test_1n() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"N"[..].to_owned());
		}
		dbg_assert_eq!(ks.contig.len(), 1);
		dbg_assert_eq!(ks.contig[0].twobit, 0);
		dbg_assert_eq!(ks.contig[0].genomic, 1);
	}
	#[test]
	fn test_1n1c1n() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"NCN"[..].to_owned());
		}
		dbg_assert_eq!(ks.contig.len(), 2);
		dbg_assert_eq!(ks.contig[0].twobit, 0);
		dbg_assert_eq!(ks.contig[0].genomic, 1);
		dbg_assert_eq!(ks.contig[1].twobit, 2);
		dbg_assert_eq!(ks.contig[1].genomic, 3);

	}
	#[test]
	fn test_17c() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"CCCCCCCCCCCCCCCCC"[..].to_owned());
		}
		dbg_assert_eq!(ks.kmp.len(), 128);
		for i in 0..ks.kmp.len() {
			dbg_assert!(ks.kmp[i].blacklisted(), "[{}], {:x}", i, ks.kmp[i]);
		}
	}
	#[test]
	fn test_1n18c1n() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned());
		}
		for i in 0..ks.kmp.len() {
			dbg_assert!(ks.kmp[i].blacklisted(), "[{}], {:x}", i, ks.kmp[i]);
		}
	}
	#[test]
	fn test_1n16c() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"NCCCCCCCCCCCCCCCC"[..].to_owned());
		}
		for i in 1..ks.kmp.len() {
			dbg_assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		let first_pos = (1 << 63) | ((kc.kmerlen as u64)<< 1);
		dbg_assert_eq!(ks.kmp[0], first_pos);
	}
	#[test]
	fn test_18at() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"ATATATATATATATATAT"[..].to_owned());
		}
		for i in 0..ks.kmp.len() {
			dbg_assert!(ks.kmp[i].blacklisted(), "[{}], {:x}", i, ks.kmp[i]);
		}
	}
	#[test]
	fn test_reconstruct1() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		let ks_kmp_len = ks.kmp.len();
		let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
		let mut kmi = KmerIter::new(&mut ks, &mut occ);
		let seq_vec =
			b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..].to_owned();
		let mut seq = seq_vec.iter();
		kmi.markcontig::<u64>(&mut seq);
		for hash in 0..ks_kmp_len {
			let mut p = kmi.ks.kmp[hash];
			if !p.blacklisted() {
				let new_stack = kmi.rebuild_kmer_stack(hash, p);
				kmi.add_newstack(new_stack, 1);
				while let Some(b2) = kmi.next_b2(&mut seq, 1) {
					if kmi.complete_occurance_or_contig(1, b2) {
						break;
					}
				}
				let mark_p = kmi.occ[1].mark.p;
				dbg_print!("testing: [{:#x}]: {:#x} == (stored p:){:#x}", hash, mark_p, p);
				dbg_assert_eq!(mark_p, p);
			}
		}
	}
	#[test]
	fn test_reconstruct_simple() { // all mappable.
		let kc = KmerConst::new(4, 15);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		let ks_kmp_len = ks.kmp.len();
		let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
		let mut kmi = KmerIter::new(&mut ks, &mut occ);
		let seq_vec = b"GGAACCTTCAGAGTG"[..].to_owned();
		let mut seq = seq_vec.iter();
		kmi.markcontig::<u64>(&mut seq);
		for hash in 0..ks_kmp_len {
			let mut p = kmi.ks.kmp[hash];
			if !p.blacklisted() {
				let new_stack = kmi.rebuild_kmer_stack(hash, p);
				kmi.add_newstack(new_stack, 1);
				while let Some(b2) = kmi.next_b2(&mut seq, 1) {
					if kmi.complete_occurance_or_contig(1, b2) {
						break;
					}
				}
				let mark_p = kmi.occ[1].mark.p;
				dbg_print!("testing: [{:#x}]: {:#x} == (stored p:){:#x}", hash, mark_p, p);
				dbg_assert_eq!(mark_p, p);
			}
		}
	}
	#[test]
	fn test_reconstruct_gs4_rl1to4_all() { // all mappable.
		let seqlen: usize = 8;

		let mut bitlen = seqlen.next_power_of_two().trailing_zeros() as usize;
		if (bitlen & 1) == 1 { // must be even.
			bitlen += 1
		}
		let kmerlen = bitlen / 2;

		for rl in kmerlen..=seqlen {
			for gen in 0..=4_usize.pow(seqlen as u32) {
				let kc = KmerConst::new(rl, seqlen);
				let mut ks = KmerStore::<u64>::new(kc.bitlen);
				let ks_kmp_len = ks.kmp.len();
				let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
				let mut kmi = KmerIter::new(&mut ks, &mut occ);
				let seq_vec:Vec<_> = (0..seqlen).map(|i| match (gen >> (i << 1)) & 3 {
				0 => 'A', 1 => 'C', 2 => 'T', 3 => 'G', _ => dbg_panic!("here")}).collect();
				dbg_print!("sequence: {:?}", seq_vec);

				let vv: Vec<u8> = seq_vec.iter().map(|&c| c as u8).collect();
				let mut seq = vv.iter();
				kmi.markcontig::<u64>(&mut seq);
				dbg_print!("-- testing hashes --");
				for hash in 0..ks_kmp_len {
					let mut p = kmi.ks.kmp[hash];
					if !p.blacklisted() {
						dbg_print!("hash: [{:#x}]: p: {:#x}", hash, p);
						let _ = kmi.rebuild_kmer_stack(hash, p);
					}
				}
			}
		}
	}

}
