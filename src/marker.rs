extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate num;
extern crate num_traits;

use crate::kmerconst::KmerConst;
use std::slice::Iter;

use crate::kmerloc::{MidPos, PriExtPosOri};
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::{anyhow, ensure, Result};

pub struct KmerIter<'a> {
    n_stretch: u64,
    goffs: u64,
    pub(super) scp: Scope<'a>,
    pub(super) ks: &'a mut KmerStore<u64>,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(ks: &'a mut KmerStore<u64>, kc: &'a KmerConst) -> Self {
        let scp = Scope::new((0, u64::max_value()), &kc, 0);
        KmerIter {
            n_stretch: 0,
            goffs: 0,
            scp,
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

    fn get_scope(&mut self, seq: &mut Iter<u8>, past: Option<&mut Scope<'a>>) -> Result<bool> {
        if let Some(scp) = past {
            return self.ks.b2_for_p(scp.p).and_then(|b2| {
                dbg_print!("=> twobit {:x} (past) <=", b2);
                if scp.complete(b2, 0) {
                    Ok(true)
                } else {
                    Err(anyhow!("cannot complete "))
                }
            });
        }
        while let Some(b2) = seq.next().map(|c| (c >> 1) & 0x7) {
            dbg_print!("=> twobit {:x} (head) <=", b2);
            let p = self.scp.p;
            if b2 < 4 {
                // new sequence is also stored, to enable lookup later.
                if let Some(qb) = self.ks.b2.get_mut(p.byte_pos()) {
                    *qb |= b2 << (p & 6);
                }
                self.ks.p_max = p.pos() + 4;
                if self.scp.complete(b2, 0) {
                    return Ok(true);
                }
            } else {
                if self.scp.i != 0 {
                    dbg_print!("started N-stretch at {}.", p);
                    self.goffs += self.scp.i as u64;
                    self.ks.push_contig(p, self.goffs);

                    // clear all except orientation and position to rebuild at the start of a new contig.
                    self.scp.i = 0;
                    self.n_stretch = 0;
                }
                self.n_stretch += 1;
            }
        }
        Ok(false) // end of contig
    }

    /// When rebuilding and exteniding scp for recurrent kmer, mind contig boundaries
    fn get_contig_limits(&self, p: u64) -> (u64, u64) {
        let (contig_start, contig_end) = self.ks.get_contig_start_end_for_p(p);
        self.scp.kc.get_kmer_boundaries(p, contig_start, contig_end)
    }

    /// rebuild scp until scp.mark.p reaches stored position.
    fn rebuild_kmer_stack(&self, stored: u64) -> Result<Scope<'a>> {
        let x = stored.x();

        let contig_limits = self.get_contig_limits(stored.pos());
        let mut scp = Scope::new(contig_limits, self.scp.kc, stored.extension());

        while {
            let b2 = self.ks.b2_for_p(scp.p).unwrap();
            dbg_print!("=> b2 {:x}, p: {:#x} <=", b2, scp.p);
            if scp.complete(b2, x) && scp.mark_is_leaving() {
                // komt voor aangezien start p soms te vroeg is. e.g. 2/3:C<AA|A>AA.
                scp.set_next_mark();
            }
            scp.mark.p != stored
        } {
            ensure!(scp.p.pos() < scp.plim, "{:#x}, {:#x}", scp.mark.p, stored);
        }
        Ok(scp)
    }

    pub fn markcontig<T: MidPos>(&mut self, seq: &mut Iter<u8>) -> Result<u64> {
        self.goffs = 0;
        self.ks.push_contig(self.scp.p.pos(), self.goffs);
        let mut past = None;

        while self.get_scope(seq, past.as_mut())? {
            loop {
                dbg_print!("---------[ past:{:?} ]-------------", past.is_some());
                let (min_index, min_pos) = if let Some(scp) = past.as_ref() {
                    (scp.mark.idx, scp.mark.p)
                } else {
                    (self.scp.mark.idx, self.scp.mark.p)
                };
                dbg_assert!(
                    min_index < self.ks.kmp.len(),
                    "{:x}, {:x}",
                    min_pos,
                    self.scp.p
                );
                dbg_assert!(min_index < self.ks.kmp.len(), "{:x}", min_pos);
                let stored_at_index = self.ks.kmp[min_index];

                if dbgx!(stored_at_index.is_replaceable_by(min_pos)) {
                    //dbg_print!("[{:#x}] (={:#x}) <= {:#x}", min_index, stored_at_index, min_pos);
                    if let Some(scp) = past {
                        dbg_print!("{}", scp);
                    } else {
                        dbg_print!("{}", self.scp);
                    }
                    if dbgx!(stored_at_index.is_set_and_not(min_pos)) {
                        let mut new_scope = self.rebuild_kmer_stack(stored_at_index)?;
                        new_scope.extend_kmer_stack(self.ks);
                        self.ks.kmp[min_index] = min_pos;
                        past = Some(new_scope);
                        // go back and extend
                    } else {
                        // If already set it was minpos. Then leave dupbit state, else it was zero.
                        self.ks.kmp[min_index] |= min_pos;
                        past = None;
                        break;
                    }
                } else if dbgx!(stored_at_index.extension() == min_pos.extension()) {
                    self.ks.kmp[min_index].set_dup();
                } else {
                    dbg_print!("\t\t<!>");
                }
                if let Some(scp) = past.as_mut() {
                    if dbgx!(!scp.extend()) {
                        break;
                    }
                    // includes check for readlength
                    ensure!(scp.set_next_mark());
                } else if dbgx!(!self.scp.extend()) {
                    // TODO: could add to blacklisted in unmappable regions.
                    break;
                }
            }
            if let Some(scp) = past.as_ref() {
                if dbgx!(scp.p.pos() >= scp.plim) {
                    // extension requires bases, but we're at contig limit, just leave it.
                    past = None;
                }
            }
        }
        dbgx!(self.finalize_n_stretch_if_pending());
        Ok(self.scp.p.pos() as u64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::KmerConst;
    use anyhow::Result;
    const READLEN: usize = 16;
    const SEQLEN: usize = 250;

    fn process<'a>(ks: &'a mut KmerStore<u64>, kc: &'a KmerConst, seq: Vec<u8>) -> Result<u64> {
        let mut kmi = KmerIter::new(ks, kc);
        kmi.markcontig::<u64>(&mut seq.iter())
    }

    #[test]
    fn test_16n() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NNNNNNNNNNNNNNNN"[..].to_owned())?;
        }
        dbg_assert_eq!(ks.contig.len(), 1);
        dbg_assert_eq!(ks.contig[0].twobit, 0);
        dbg_assert_eq!(ks.contig[0].genomic, 16);
        Ok(())
    }
    #[test]
    fn test_1n() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"N"[..].to_owned())?;
        }
        dbg_assert_eq!(ks.contig.len(), 1);
        dbg_assert_eq!(ks.contig[0].twobit, 0);
        dbg_assert_eq!(ks.contig[0].genomic, 1);
        Ok(())
    }
    #[test]
    fn test_1n1c1n() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCN"[..].to_owned())?;
        }
        dbg_assert_eq!(ks.contig.len(), 2);
        dbg_assert_eq!(ks.contig[0].twobit, 0);
        dbg_assert_eq!(ks.contig[0].genomic, 1);
        dbg_assert_eq!(ks.contig[1].twobit, 2);
        dbg_assert_eq!(ks.contig[1].genomic, 3);
        Ok(())
    }
    #[test]
    fn test_17c() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"CCCCCCCCCCCCCCCCC"[..].to_owned())?;
        }
        dbg_assert_eq!(ks.kmp.len(), 128);
        for i in 0..ks.kmp.len() {
            dbg_assert!(ks.kmp[i].no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
        Ok(())
    }
    #[test]
    fn test_1n18c1n() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned())?;
        }
        for i in 0..ks.kmp.len() {
            dbg_assert!(ks.kmp[i].no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
        Ok(())
    }
    #[test]
    fn test_1n16c() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCC"[..].to_owned())?;
        }
        for i in 1..ks.kmp.len() {
            dbg_assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
        }
        let first_pos = (1 << 63) | ((kc.kmerlen as u64) << 1);
        dbg_assert_eq!(ks.kmp[0], first_pos);
        Ok(())
    }
    #[test]
    fn test_18at() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"ATATATATATATATATAT"[..].to_owned())?;
        }
        for i in 0..ks.kmp.len() {
            dbg_assert!(ks.kmp[i].no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
        Ok(())
    }
    #[test]
    fn test_reconstruct1() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        let ks_kmp_len = ks.kmp.len();
        let mut kmi = KmerIter::new(&mut ks, &kc);
        let seq_vec = b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..]
            .to_owned();
        let mut seq = seq_vec.iter();
        kmi.markcontig::<u64>(&mut seq)?;
        for hash in 0..ks_kmp_len {
            let p = kmi.ks.kmp[hash];
            if !p.no_pos() {
                let new_stack = kmi.rebuild_kmer_stack(p).unwrap();
                dbg_assert_eq!(new_stack.mark.p, p);
            }
        }
        Ok(())
    }
    #[test]
    fn test_reconstruct_gs4_rl1to4_all() -> Result<()> {
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
                let mut kmi = KmerIter::new(&mut ks, &kc);
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
                kmi.markcontig::<u64>(&mut seq)?;
                dbg_print!("-- testing hashes --");
                for hash in 0..ks_kmp_len {
                    let p = kmi.ks.kmp[hash];
                    if !p.no_pos() {
                        dbg_print!("hash: [{:#x}]: p: {:#x}", hash, p);
                        dbg_assert_eq!(kmi.rebuild_kmer_stack(p).unwrap().mark.p, p);
                    }
                }
            }
        }
        Ok(())
    }
}
