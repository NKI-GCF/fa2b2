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
    //steekproefi: u32,
    n_stretch: u64,
    goffs: u64,
    repetitive: Option<u64>,
    pub(super) scp: Scope<'a>,
    pub(super) ks: &'a mut KmerStore<u64>,
    pub(super) kc: &'a KmerConst,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(ks: &'a mut KmerStore<u64>, kc: &'a KmerConst) -> Self {
        let scp = Scope::new((0, u64::max_value()), &kc, 0);
        KmerIter {
            //steekproefi: 10_000,
            n_stretch: 0,
            goffs: 0,
            repetitive: None,
            scp,
            ks,
            kc,
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

    // false
    fn build_scope(&mut self, seq: &mut Iter<u8>, past: Option<&mut Scope<'a>>) -> Result<bool> {
        if let Some(scp) = past {
            return self.ks.b2_for_p(scp.p).and_then(|b2| {
                dbg_print!("=> twobit {:x} (past) <=", b2);
                scp.complete_and_update_mark(b2, 0)
            });
            // false: cannot complete scope (contig end?) or find next mark after leaving mark.
        }
        while let Some(b2) = seq.next().map(|c| (c >> 1) & 0x7) {
            dbg_print!("=> twobit {:x} (head) <=", b2);
            let p = self.scp.p;
            if b2 < 4 {
                self.finalize_n_stretch_if_pending();
                // new sequence is also stored, to enable lookup later.
                if let Some(qb) = self.ks.b2.get_mut(p.byte_pos()) {
                    *qb |= b2 << (p & 6);
                }
                self.ks.p_max = p.pos() + 4;
                // ext cannot be unset downstream: scp funcs also used for scope rebuild.
                self.scp.p.set_extension(0);
                return self.scp.complete_and_update_mark(b2, 0);
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
        Err(anyhow!("end of contig."))
    }

    /// When rebuilding and exteniding scp for recurrent kmer, mind contig boundaries
    fn get_contig_limits(&self, p: u64) -> (u64, u64) {
        let p = p.pos();
        let (contig_start, contig_end) = self.ks.get_contig_start_end_for_p(p);
        self.kc.get_kmer_boundaries(p, contig_start, contig_end)
    }

    /// rebuild scp until scp.mark.p reaches stored position.
    fn rebuild_scope(&self, stored: u64) -> Result<Option<Scope<'a>>> {
        ensure!(stored.pos() != PriExtPosOri::no_pos());
        let x = stored.x();

        let contig_limits = self.get_contig_limits(stored);
        let mut scp = Scope::new(contig_limits, self.kc, stored.extension());

        loop {
            let b2 = self.ks.b2_for_p(scp.p).unwrap();
            dbg_print!("=> b2 {:x}, p: {:#x} <=", b2, scp.p);
            if scp.complete_and_update_mark(b2, x)? {
                // komt voor aangezien start p soms te vroeg is. e.g. 2/3:C<AA|A>AA.
                scp.set_next_mark();
            }
            if scp.mark.p == stored {
                break;
            }
            if scp.p.pos() >= scp.plim.1 {
                return Ok(None);
            }
        }
        Ok(Some(scp))
    }

    fn set_idx_pos(&mut self, idx: usize, p: u64) {
        /*self.steekproefi -= 1;
        if self.steekproefi == 0 {
            self.steekproefi = 10_000;
            // XXX waarom geeft dit telkens een p met dupbit set en ext == 3 ???
            eprintln!("[{:#x}] = {:#x}", idx, p);
        }*/
        self.ks.kmp[idx].set(p);
    }

    pub fn markcontig<T: MidPos>(&mut self, seq: &mut Iter<u8>) -> Result<u64> {
        self.goffs = 0;
        self.ks.push_contig(self.scp.p.pos(), self.goffs);
        let mut past = None;

        loop {
            match self.build_scope(seq, past.as_mut()) {
                Ok(true) => {}
                Ok(false) => continue,
                Err(e) => {
                    // fixme: use thiserror?
                    if e.to_string() == "end of contig." {
                        break;
                    }
                    if past.is_some() {
                        dbg_print!("{} for past (p:{:x})", e, past.unwrap().mark.p);
                        past = None;
                    } else {
                        dbg_print!("{} (p:{:x})", e, self.scp.mark.p);
                    }
                    continue;
                }
            }
            loop {
                /*dbg_print!(
                    "---------[ past:{:?} ]-------------",
                    self.scp[1].p.is_set()
                );*/
                let (min_idx, min_p) = if let Some(scp) = past.as_ref() {
                    (scp.mark.idx, scp.mark.p)
                } else {
                    (self.scp.mark.idx, self.scp.mark.p)
                };
                dbg_assert!(min_idx < self.ks.kmp.len(), "{:x}, {:x}", min_idx, min_p);
                let stored_p = self.ks.kmp[min_idx];

                if dbgx!(stored_p.is_replaceable_by(min_p)) {
                    //dbg_print!("[{:#x}] (={:#x}) <= {:#x}", min_idx, stored_p, min_p);
                    // dbg_print!("{}", self.get_scp());
                    if dbgx!(stored_p.is_set_and_not(min_p)) {
                        if let Some(mut new_scope) = self.rebuild_scope(stored_p)? {
                            new_scope.extend_kmer_stack(self.ks);
                            dbgx!(self.set_idx_pos(min_idx, min_p));
                            past = Some(new_scope);
                        } else {
                            past = None;
                            break;
                        }
                        //dbgx!(self.ks.kmp[min_idx].set(min_p));
                        // go back and extend
                    } else {
                        // If already set it was min_p. Then leave dupbit state.
                        if self.ks.kmp[min_idx].is_no_pos() {
                            self.set_idx_pos(min_idx, min_p);
                        }
                        past = None;
                        break;
                    }
                } else if dbgx!(stored_p.extension() == min_p.extension()) {
                    // If a kmer occurs multiple times within an extending readlength, only
                    // the first gets a position. During mapping this rule also should apply.
                    let check_linked = if let Some(scp) = past.as_mut() {
                        scp.downstream_on_contig(stored_p)
                    } else {
                        self.scp.downstream_on_contig(stored_p)
                    };
                    if check_linked {
                        if let Some(mut new_scope) = self.rebuild_scope(stored_p)? {
                            if new_scope.in_linked_scope(self.ks, min_p, min_idx)? {
                                past = None;
                                break;
                            }
                        }
                    }
                    self.ks.kmp[min_idx].set_dup();
                } else {
                    dbg_print!("\t\t<!>");
                }
                if let Some(scp) = past.as_mut() {
                    if dbgx!(!scp.extend()) {
                        break;
                    }
                    ensure!(scp.set_next_mark());
                } else if dbgx!(!self.scp.extend()) {
                    break;
                }
                // includes check for readlength
            }
            if let Some(scp) = past.as_ref() {
                if dbgx!(scp.p.pos() >= scp.plim.1) {
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
        let mut observed_first = false;
        for i in 1..ks.kmp.len() {
            dbg_assert!(ks.kmp[i].is_no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
        let first_pos = (kc.kmerlen as u64) << 1;
        dbg_assert_eq!(ks.kmp[0], first_pos);
        Ok(())
    }
    #[test]
    fn test_1n18c1n() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned())?;
        }
        for i in 1..ks.kmp.len() {
            dbg_assert!(ks.kmp[i].is_no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
        let first_pos = (kc.kmerlen as u64) << 1;
        dbg_assert_eq!(ks.kmp[0], first_pos);
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
            dbg_assert!(ks.kmp[i].is_no_pos(), "[{}], {:x}", i, ks.kmp[i]);
        }
        let first_pos = (kc.kmerlen as u64) << 1;
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
            dbg_assert!(
                ks.kmp[i].is_no_pos() || !ks.kmp[i].is_dup(),
                "[{}], {:x}",
                i,
                ks.kmp[i]
            );
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
            if p.is_set() {
                let scope = kmi.rebuild_scope(p).unwrap();
                dbg_assert_eq!(scope.mark.p, p);
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
                    if p.is_set() {
                        dbg_print!("hash: [{:#x}]: p: {:#x}", hash, p);
                        dbg_assert_eq!(kmi.rebuild_scope(p).unwrap().mark.p, p);
                    }
                }
            }
        }
        Ok(())
    }
}
