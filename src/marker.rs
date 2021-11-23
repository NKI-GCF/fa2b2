extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate num;
extern crate num_traits;

use std::cmp;
use std::slice::Iter;

use crate::kmerloc::{MidPos, PriExtPosOri};
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::{anyhow, Result};

pub struct KmerIter<'a> {
    n_stretch: u64,
    goffs: u64,
    pub(super) occ: &'a mut Vec<Scope<'a>>,
    pub(super) ks: &'a mut KmerStore<u64>,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(ks: &'a mut KmerStore<u64>, occ: &'a mut Vec<Scope<'a>>) -> Self {
        KmerIter {
            n_stretch: 0,
            goffs: 0,
            occ,
            ks,
        }
    }

    /// indien we een n-stretch aan het verlengen waren, sluit deze af.
    // XXX de n_stretch logica hier lijkt borked.
    fn finalize_n_stretch_if_pending(&mut self) {
        if self.n_stretch > 0 {
            dbg_print!("added new contig. Ns:{}", self.n_stretch);
            self.ks.offset_contig(self.n_stretch);
            self.goffs += self.n_stretch;
            self.n_stretch = 0;
        }
    }

    /// if N, insert contig (once) and count stretch => false
    /// else store and update scope kmers, complete ? true : false
    fn complete_scope_or_contig(&mut self, n: usize, b2: u8) -> bool {
        let mut ret = false;
        let p = self.occ[n].p.pos();
        if b2 < 4 {
            if let Some(qb) = self.ks.b2.get_mut(p as usize >> 3) {
                *qb |= b2 << (p & 6);
            }
            if self.occ[n].complete(b2, 0) {
                if dbgx!(self.occ[n].mark_is_leaving()) {
                    ret = self.occ[n].set_next_mark();
                } else if n == 0 {
                    let mark_idx = self.occ[n].mark.idx;
                    self.occ[n].p.clear_extension();
                    ret = self.occ[n].set_next_mark();

                    // een poging om deze branch overbodig te maken.
                    dbg_assert_eq!(mark_idx, self.occ[n].mark.idx);
                } else {
                    ret = true;
                }
            } else {
                self.finalize_n_stretch_if_pending();
            }
            self.ks.p_max = p + 4; // == self.occ[n].p.pos() + 2;
        } else if self.occ[n].i != 0 {
            self.goffs += self.occ[n].i as u64;
            dbg_print!("started N-stretch at {}.", p);
            self.ks.push_contig(p, self.goffs);

            // clear all except orientation and position to rebuild at the start of a new contig.
            //dbg_assert_eq!(self.occ[n].i, 1);
            self.occ[n].i = 0;
            self.n_stretch = 1;
        } else {
            self.n_stretch += 1;
        }
        ret
    }

    /// when rebuilding occ for recurrent kmer, and extending take into account contig boundaries
    /// for that site
    fn get_comtig_limits(&self, p: u64) -> (u64, u64) {
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
        dbgf!(
            left,
            "{:?}, p:{}, contig_start:{}, p_rl:{}",
            p,
            contig_start,
            p_rl
        );

        let right = cmp::min(
            p + ((kc.readlen - kc.kmerlen) << 1) as u64,
            self.ks.get_twobit_after(i).unwrap_or(p_max),
        );
        dbgf!(
            right,
            "{:?}, p_max:{}",
            self.ks.get_twobit_after(i).unwrap_or(p_max)
        );

        (left, right)
    }

    /// rebuild occ until occ.mark.p reaches stored position.
    fn rebuild_kmer_stack(&self, min_index: usize) -> Result<Scope<'a>> {
        let stored = self.ks.kmp[min_index];
        let x = stored.x();

        let comtig_limits = self.get_comtig_limits(stored.pos());
        let mut occ = Scope::new(comtig_limits, self.occ[0].kc, stored.extension());

        while {
            let p = occ.p.pos();
            let b2 = self.ks.b2_for_p(p)?;
            dbg_print!("=> b2 {:x}, p: {:#x} <=", b2, p);
            if occ.complete(b2, x) && occ.mark_is_leaving() {
                // komt voor aangezien start p soms te vroeg is. e.g. 2/3:C<AA|A>AA.
                occ.set_next_mark();
            }
            occ.mark.p != stored
        } {
            dbg_assert!(occ.p.pos() < occ.plim, "{:#x}, {:#x}", occ.mark.p, stored);
        }
        Ok(occ)
    }

    /// continue occ rebuild
    fn extend_kmer_stack(&self, occ: &mut Scope<'a>) {
        let x = occ.p.x();
        occ.set_mark_after_extension_if_possible(x);
        while occ.p.pos() < occ.plim {
            let p = occ.p.pos();
            let b2 = self.ks.b2_for_p(p).unwrap();
            dbg_print!("-> b2 {:x}, p: {:#x} <-", b2, p);
            if occ.complete(b2, x) && occ.mark_is_leaving() {
                let mark_is_set = occ.set_next_mark();
                dbg_assert!(mark_is_set);
                if self.ks.kmp[occ.mark.idx] == occ.mark.p {
                    // herspoord = klaar
                    break;
                }
            }
        }
    }

    fn search_occ_for_pos(&self, original_n: usize, stored_at_index: u64) -> usize {
        // XXX: not understood, but this seems to occur.
        //dbg_assert!(original_n == self.occ.len() - 1, "{}, {}", original_n, self.occ.len() - 1);

        for n in 0..original_n {
            if dbgf!(
                stored_at_index == self.occ[n].mark.p,
                "{:#?}: {:#x} == {:#x}?",
                stored_at_index,
                self.occ[n].mark.p
            ) {
                return dbgx!(n);
            }
        }
        original_n
    }
    fn add_newstack(&mut self, next_stack: Scope<'a>, n: usize) {
        if n < self.occ.len() {
            self.occ[n] = next_stack;
        } else {
            self.occ.push(next_stack);
        }
    }

    /// when the stored entry in ks.kmp is replaced or blacklisted, the occurence for the stored
    /// needs to be extended to find, there, a minpos for the extended kmers.
    /// the entry is added to the stack, return the index.
    fn get_next_for_extension(
        &mut self,
        min_index: usize,
        new_p: u64,
        original_n: usize,
    ) -> Result<usize> {
        let stored = self.ks.kmp[min_index];
        let mut n = self.search_occ_for_pos(original_n, stored);
        if dbgx!(stored != self.occ[n].mark.p) {
            let mut next_stack = self.rebuild_kmer_stack(min_index)?;
            dbg_assert!(next_stack.mark.idx == min_index);
            let is_blacklisted = new_p.no_pos();
            dbgf!(
                self.ks.kmp[min_index] = new_p,
                "{:?}, kmp[{:#x}] = {:#x}",
                min_index,
                new_p
            );
            self.extend_kmer_stack(&mut next_stack);

            if dbgx!(!is_blacklisted && n != 0) {
                // position is written below, this self.occ[n] is done.
                // Could add it to recurring kmer_stacks, though.
                if next_stack.p.pos() < next_stack.plim {
                    self.occ[n] = next_stack;
                } else {
                    n -= 1;
                }
            } else if dbgx!(next_stack.p.pos() < next_stack.plim) {
                // 0th is never overwritten.
                n += 1;
                self.add_newstack(next_stack, n);
            }
        }
        Ok(n)
    }

    fn next_b2(&self, seq: &mut Iter<u8>, n: usize) -> Result<u8> {
        //dbgf!(self.occ[n].p, "{:#x} n:{}", n);
        if n == 0 {
            seq.next()
                .map(|c| (c >> 1) & 7u8)
                .ok_or(anyhow!("End of seq"))
        } else {
            self.ks.b2_for_p(self.occ[n].p.pos())
        }
    }

    pub fn markcontig<T: MidPos>(&mut self, seq: &mut Iter<u8>) -> Result<u64> {
        let mut n = 0;
        self.goffs = 0;
        self.ks.push_contig(self.occ[n].p.pos(), self.goffs);

        while let Ok(b2) = self.next_b2(seq, n) {
            dbg_print!("=> twobit {:x} (n: {}) <=", b2, n);

            if !self.complete_scope_or_contig(n, b2) {
                continue;
            }
            //dbgf!(self.occ[n].p, "{:#x}, mark.p:{:#x}", self.occ[n].mark.p);
            loop {
                dbg_print!("---------[ n: {} ]-------------", n);

                let (min_index, mut min_pos) = (self.occ[n].mark.idx, self.occ[n].mark.p);
                dbg_assert!(
                    min_index < self.ks.kmp.len(),
                    "{:x}, {:x}",
                    min_pos,
                    self.occ[n].p
                );
                let stored_at_index = self.ks.kmp[min_index];

                if dbgx!(stored_at_index.is_replaceable_by(min_pos)) {
                    //dbg_print!("[{:#x}] (={:#x}) <= {:#x}", min_index, stored_at_index, min_pos);
                    dbg_print!("{}", self.occ[n]);

                    if dbgx!(stored_at_index.is_set_and_not(min_pos)) {
                        let orig_n = n;
                        n = self.get_next_for_extension(min_index, min_pos, n)?;
                        if n == orig_n {
                            continue;
                        }
                    } else {
                        self.ks.kmp[min_index] = min_pos;
                        if n > 0 {
                            // We're done on a recurrence when observing an already set value, not mark.p
                            if dbgf!(
                                stored_at_index.is_same(min_pos)
                                    || self.occ[n].p.pos() == self.occ[n].plim,
                                "{:?}, {:#x}, {:#x}",
                                self.occ[n].p,
                                self.occ[n].plim
                            ) {
                                // TODO could add [n] to another stack for fast lookup for multimappers.
                                n -= 1;
                                //dbg_print!("n:{} x:{:#x} p:{:#x} mark.p:{:#x}", n,
                                //self.occ[n].mark.idx, self.occ[n].p, self.occ[n].mark.p);
                                continue;
                            }
                            //dbg_assert!(self.occ[n].p.pos() <= self.occ[n].plim);
                        }
                    }

                    break; // position written (done) or added next_stack which requires extension.
                } else {
                    dbg_print!("{}\t\t<!>", self.occ[n]);
                }
                //dbgf!(self.occ[n].mark.p, "{:x}");
                if dbgx!(self.occ[n].extend()) {
                    let _long_enough = self.occ[n].set_next_mark();
                    dbg_assert!(_long_enough);

                    if dbgx!(stored_at_index.extension() == min_pos.extension()) {
                        // both stored and current require extension. Stored is handled now, new is postponed.
                        min_pos.blacklist();
                        let orig_n = n;
                        n = self.get_next_for_extension(min_index, min_pos, n)?;
                        dbg_print!("{}\t\t<0>", self.occ[n]);
                        if dbgx!(n == orig_n) {
                            self.occ[n].p.clear_extension();
                            let _ = self.occ[n].set_next_mark();
                            dbg_assert_eq!(self.occ[n].mark.idx, min_index);
                            dbgf!(self.occ[n].mark.idx, "\n[{:#x}]");
                            continue;
                        }
                        break;
                    }
                } else {
                    dbg_print!("{}\t\t<XXX>", self.occ[n]);
                    self.ks.kmp[min_index].blacklist();
                    dbgf!(
                        min_index,
                        "[{:#x}] is blacklisted ({:#x})",
                        self.ks.kmp[min_index]
                    );
                    if dbgx!(n == 0) {
                        // 0th needs to be renewed when done.
                        //dbgf!(self.occ[n].mark.p, "mark.p:{:x}, p:{:x}", self.occ[n].p);
                        break;
                    }
                    //let blacklisted = self.occ.pop();
                    // TODO: could add [n] to blacklisted in unmappable regions.
                    n -= 1;
                    if dbgx!(n == 0 && self.occ[n].mark.p == self.ks.kmp[self.occ[n].mark.idx]) {
                        // test already done
                        break;
                    }
                }
            }
        }
        dbgx!(self.finalize_n_stretch_if_pending());
        Ok(self.occ[0].p.pos() as u64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::KmerConst;
    const READLEN: usize = 16;
    const SEQLEN: usize = 250;

    fn process<'a>(occ: &'a mut Vec<Scope<'a>>, ks: &'a mut KmerStore<u64>, seq: Vec<u8>) {
        let mut kmi = KmerIter::new(ks, occ);
        kmi.markcontig::<u64>(&mut seq.iter());
    }

    #[test]
    fn test_16n() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
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
            let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
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
            let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
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
            let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
            process(&mut occ, &mut ks, b"CCCCCCCCCCCCCCCCC"[..].to_owned());
        }
        dbg_assert_eq!(ks.kmp.len(), 128);
        for i in 0..ks.kmp.len() {
            dbg_assert!(ks.kmp[i].no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
    }
    #[test]
    fn test_1n18c1n() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
            process(&mut occ, &mut ks, b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned());
        }
        for i in 0..ks.kmp.len() {
            dbg_assert!(ks.kmp[i].no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
    }
    #[test]
    fn test_1n16c() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
            process(&mut occ, &mut ks, b"NCCCCCCCCCCCCCCCC"[..].to_owned());
        }
        for i in 1..ks.kmp.len() {
            dbg_assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
        }
        let first_pos = (1 << 63) | ((kc.kmerlen as u64) << 1);
        dbg_assert_eq!(ks.kmp[0], first_pos);
    }
    #[test]
    fn test_18at() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
            process(&mut occ, &mut ks, b"ATATATATATATATATAT"[..].to_owned());
        }
        for i in 0..ks.kmp.len() {
            dbg_assert!(ks.kmp[i].no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
    }
    #[test]
    fn test_reconstruct1() {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        let ks_kmp_len = ks.kmp.len();
        let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
        let mut kmi = KmerIter::new(&mut ks, &mut occ);
        let seq_vec = b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..]
            .to_owned();
        let mut seq = seq_vec.iter();
        kmi.markcontig::<u64>(&mut seq);
        for hash in 0..ks_kmp_len {
            let mut p = kmi.ks.kmp[hash];
            if !p.no_pos() {
                let new_stack = kmi.rebuild_kmer_stack(hash).unwrap();
                dbg_assert_eq!(new_stack.mark.p, p);
            }
        }
    }
    #[test]
    fn test_reconstruct_simple() {
        // all mappable.
        let kc = KmerConst::new(4, 15);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        let ks_kmp_len = ks.kmp.len();
        let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
        let mut kmi = KmerIter::new(&mut ks, &mut occ);
        let seq_vec = b"GGAACCTTCAGAGTG"[..].to_owned();
        let mut seq = seq_vec.iter();
        kmi.markcontig::<u64>(&mut seq);
        for hash in 0..ks_kmp_len {
            let mut p = kmi.ks.kmp[hash];
            if !p.no_pos() {
                let new_stack = kmi.rebuild_kmer_stack(hash).unwrap();
                kmi.add_newstack(new_stack, 1);
                while let Ok(b2) = kmi.next_b2(&mut seq, 1) {
                    if kmi.complete_scope_or_contig(1, b2) {
                        break;
                    }
                }
                let mark_p = kmi.occ[1].mark.p;
                dbg_print!(
                    "testing: [{:#x}]: {:#x} == (stored p:){:#x}",
                    hash,
                    mark_p,
                    p
                );
                dbg_assert_eq!(mark_p, p);
            }
        }
    }
    #[test]
    fn test_reconstruct_gs4_rl1to4_all() {
        // all mappable.
        let seqlen: usize = 6; //8;

        let mut bitlen = seqlen.next_power_of_two().trailing_zeros() as usize;
        if (bitlen & 1) == 1 {
            // must be even.
            bitlen += 1
        }
        let kmerlen = bitlen / 2;

        for rl in kmerlen..=seqlen {
            for gen in 0..=4_usize.pow(seqlen as u32) {
                let kc = KmerConst::new(rl, seqlen);
                let mut ks = KmerStore::<u64>::new(kc.bitlen);
                let ks_kmp_len = ks.kmp.len();
                let mut occ: Vec<Scope> = vec![Scope::new((0, u64::max_value()), &kc, 0)];
                let mut kmi = KmerIter::new(&mut ks, &mut occ);
                let seq_vec: Vec<_> = (0..seqlen)
                    .map(|i| match (gen >> (i << 1)) & 3 {
                        0 => 'A',
                        1 => 'C',
                        2 => 'T',
                        3 => 'G',
                        _ => dbg_panic!("here"),
                    })
                    .collect();
                dbg_print!("[{:#x}] sequence:\n{:?}", gen, seq_vec);

                let vv: Vec<u8> = seq_vec.iter().map(|&c| c as u8).collect();
                let mut seq = vv.iter();
                kmi.markcontig::<u64>(&mut seq);
                dbg_print!("-- testing hashes --");
                for hash in 0..ks_kmp_len {
                    let mut p = kmi.ks.kmp[hash];
                    if !p.no_pos() {
                        dbg_print!("hash: [{:#x}]: p: {:#x}", hash, p);
                        dbg_assert_eq!(kmi.rebuild_kmer_stack(hash).unwrap().mark.p, p);
                    }
                }
            }
        }
    }
}
