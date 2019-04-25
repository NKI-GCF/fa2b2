use std::collections::VecDeque;

use kmer::{Kmer,test_template};
pub use kmerconst::KmerConst;
use kmerloc::{KmerLoc,PriExtPosOri};
use rdbg::STAT_DB;
use std::fmt;

// rename to Recurrance?
pub struct Occurrence<'a> {
	pub kc: &'a KmerConst,
	pub p: u64,
	d: VecDeque<Kmer<u64>>,	   // misschien is deze on the fly uit ks te bepalen?
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

	/// hierin vinden .i & .p increments and kmer .d[] update plaats. 1+ kmers ? true.
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
		//dbgf!(self.p, "p_add:{} i:{}, d[i]:{:#x} rc:{:#x} p:{:#x}", new_idx, self.d[new_idx].dna, self.d[new_idx].rc, self.p);
		self.i += 1;
		self.i >= self.kc.kmerlen
	}

	/// add b2 to kmer, move .p according to occ orientation
	/// return whether required nr of kmers for struct occurance were seen, since contig start.
	/// adds 2bit to stored sequence
	pub fn complete(&mut self, b2: u8, x_start: usize) -> bool {

		if !self.complete_kmer(b2) { //XXX: hierna is self.i al geincrement.
			return false;
		}
		let x_end = self.p.x() + 1;
		//dbgf!(x_start, "{}..{}", x_end);

		for x in x_start..x_end {
			if self.i < self.kc.kmerlen + self.kc.afstand(x) || self.set_if_extreme(0, x) {
				break;
			}
		}
		self.mark.is_set() && self.all_kmers()
	}

	pub fn set_mark_after_extension_if_possible(&mut self, x: usize) {
		self.mark.reset();
		let _ = self.i < self.kc.kmerlen + self.kc.afstand(x) || self.set_if_extreme(0, x);
	}

	/// geef de current kmer.
	pub fn kmer(&self) -> Kmer<u64> {
		// Todo: on the fly uit ks.
		self.d[(self.i - self.kc.kmerlen) % self.kc.no_kmers]
	}

	/// extend positie (als kmer/hash niet replaceble was); true indien mogelijk.
	// waar wordt dit ongedaan gemaakt??
	pub fn extend(&mut self) -> bool {
		if self.p.x() + 1 < self.kc.ext_max {
			self.p.extend();
			true
		} else {false}
	}

	/// is occ complete? call na complete_kmer() - self.i increment.
	pub fn all_kmers(&self) -> bool {
		self.i >= self.kc.readlen
	}

	/// true als voor de occurance de hash een nieuw minimum/optimum is.
	fn hash_is_extreme(&mut self, hash: usize, p: u64) -> bool {
		let xh = hash ^ self.kc.ext_domain(p.x());
		//dbgf!(p < self.mark.p, "{:?}, p:{:x}, mark.p:{:#x}", p, self.mark.p))
		xh < self.mark.idx || (xh == self.mark.idx && p < self.mark.p)
	}

	/// is de minimum/optimum leaving?
	pub fn mark_is_leaving(&self) -> bool {
		//let afs = self.kc.afstand(self.mark.p.x());
		let p_pos = self.p.pos();
		let dist = (self.kc.no_kmers << 1) as u64;
		//dbgf!(p_pos, "{} - (({} + {}) << 1) == {} ?", self.kc.no_kmers, afs, self.mark.p.pos());
		p_pos == self.mark.p.pos() + dist
	}

	/// bepaal van alle kmers het nieuwe minimum/optimum (na leaving mark of extensie)
	pub fn set_next_mark(&mut self) -> bool {

		self.mark.reset();
		let x = self.p.x();
		let afs = self.kc.afstand(x);
		//dbg_print!("no_kmers: {}, x:{}, afstand: {}", self.kc.no_kmers, x, afs);
		if self.kc.no_kmers <= afs {
			return dbgx!(false);
		}
		for i in 0..(self.kc.no_kmers - afs) {
			let _ = self.set_if_extreme(i, x);
		}
		//dbg_print!("mark was set: [{:#x}]:{:#x}", self.mark.idx, self.mark.p);
		dbg_assert!(self.mark.is_set());
		true
	}

	/// voor een offset i en extensie x, maak de kmer/hash en zet mark + return true als extreme.
	fn set_if_extreme(&mut self, i: usize, x: usize) -> bool
	{
		let base = self.i - self.kc.kmerlen;
		let afs = self.kc.afstand(x);
		let d_i = base.wrapping_sub(i) % self.kc.no_kmers;
		let mut kmer = self.d[d_i];
		//dbg_print!("base:{}, i:{} d_i:{} idx:{:#x}, x:{}", base, i, d_i, self.d[d_i].get_idx(true), x);
		if x > 0 {
			let d_i2 = base.wrapping_sub(afs + i) % self.kc.no_kmers;
			dbg_print!("dna:{:#x}, dna2:{:#x}", self.d[d_i2].dna, kmer.dna);
			//dbgf!(base.wrapping_sub(afs + i) % self.kc.no_kmers, "[{}] ^ [{}] (afs: {})", i % self.kc.no_kmers, afs);
			kmer.hash(self.d[d_i2]);
		}
		let hash = kmer.get_idx(true);
		let mut p = self.p.with_ext(x) - ((afs + i) << 1) as u64;
		//dbg_print!("p:{:#x}, afs:{:#x}, i:{:#x} dna:{:#x}, rc:{:#x}", p, afs, i, kmer.dna, kmer.rc);
		if test_template(kmer.dna, kmer.rc) {
			p |= 1;
		}
		if self.hash_is_extreme(hash, p) {
			//dbgf!(hash, "[{:x}] = {:x} = {:x} - (({:x} + {:x}) << 1) | 1?", p, self.p, afs, i);
			self.mark.set(hash, p, x);
			//dbgf!(self.mark.set(hash, p, x), "{:?} [{:#x}]: {:#x}", self.mark.idx, self.mark.p);
			true
		} else {
			false
		}
	}
}

impl<'a> fmt::Display for Occurrence<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		let p = self.p.pos() as usize;
		let mp = self.mark.p.pos() as usize;
		let n = self.kc.kmerlen + self.p.x();
		let o = " ".repeat(((p >> 1) - self.kc.readlen) * 5);
		let r = (p - mp) >> 1;
		if r == 0 {
			let x = self.kc.readlen - n;
			let s = if x != 0 {" ".repeat(x << 2) + "|"} else {String::from("")};
			write!(f, "{2}<{3}{: ^1$x}>", self.mark.idx, n << 2, o, s)
		} else if r + self.kc.kmerlen == self.kc.readlen {
			let x = self.kc.readlen - n;
			let s = if x != 0 {String::from("|") + &" ".repeat(x << 2)} else {String::from("")};
			write!(f, "{2}<{: ^1$x}{3}>", self.mark.idx, n << 2, o, s)
		} else {
			//let l = self.kc.readlen - r - n;
			//let ls = if o {" ".repeat(o) + "|"} else {String::from("")};
			//let rs = if l {String::from("|") + &" ".repeat(l << 2)} else {String::from("")};
			//write!(f, "{2}<{3}|{: ^1$x}|{4}>", self.mark.idx, n << 2, o, ls, rs)
			write!(f, "")
		}
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

