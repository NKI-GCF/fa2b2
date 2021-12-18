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
    period: Option<u64>,
    pub(super) scp: [Scope<'a>; 2],
    pub(super) ks: &'a mut KmerStore<u64>,
    pub(super) kc: &'a KmerConst,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(ks: &'a mut KmerStore<u64>, kc: &'a KmerConst) -> Self {
        let scp = [
            Scope::new((0, u64::max_value()), kc, 0),
            Scope::new((0, u64::max_value()), kc, PriExtPosOri::no_pos()),
        ];
        KmerIter {
            //steekproefi: 10_000,
            n_stretch: 0,
            goffs: 0,
            period: None,
            scp,
            ks,
            kc,
        }
    }

    /// sluit een n-stretch af, die we aan het verlengen waren.
    fn finalize_n_stretch(&mut self) {
        dbg_print!("added new contig. Ns:{}", self.n_stretch);
        self.ks.offset_contig(self.n_stretch);
        self.goffs += self.n_stretch;
        self.n_stretch = 0;
    }

    /// Mark / store recurring xmers. Skippable if repetitive on contig. Else mark as dup.
    fn recurrance_is_skippable(&mut self, p: u64, min_idx: usize) -> bool {
        let scp = self.get_scp();
        let pos = p.pos();
        if pos >= scp.plim.0 {
            let dist = dbgx!(scp.mark.p.pos() - pos);
            if let Some(period) = self.period {
                // Repeats can be more complex than a regular repetition. e.g.
                // if a transposon is inserted + inserted inverted, and this is the
                // repeat, then an xmers could recur at irregular intervals.
                // but then the xmer out of line is extended and gets its own period,
                // so it resolves that way too.
                if dist % period == 0 {
                    // XXX A xmer could recur exacly on period but not be a(n exact) repeat.
                    self.ks.extend_repetitive(min_idx, dist as u32);
                    return true;
                }
            } else {
                self.period = Some(dist);
                return true;
            }
        } else {
            self.period = None;
        }
        self.ks.kmp[min_idx].set_dup();
        false
    }

    /// When rebuilding and exteniding scp for recurrent kmer, mind contig boundaries
    fn get_contig_limits(&self, p: u64) -> (u64, u64) {
        let p = p.pos();
        let (contig_start, contig_end) = self.ks.get_contig_start_end_for_p(p);
        self.kc.get_kmer_boundaries(p, contig_start, contig_end)
    }

    /// rebuild scp until scp.mark.p reaches stored position.
    fn rebuild_scope(&mut self, stored: u64) -> Result<bool> {
        ensure!(stored.pos() != PriExtPosOri::no_pos());
        let x = stored.x();
        let contig_limits = self.get_contig_limits(stored);
        let scp = &mut self.scp[1];

        scp.set(contig_limits, stored.extension());

        loop {
            let b2 = self.ks.b2_for_p(scp.p).unwrap();
            dbg_print!("=> b2 {:x}, p: {:#x} <=", b2, scp.p);
            if scp.complete_and_update_mark(b2, x)? {
                // komt voor aangezien start p soms te vroeg is. e.g. 2/3:C<AA|A>AA.
                scp.set_next_mark()?;
            }
            if scp.mark.p == stored {
                break;
            }
            if scp.p.pos() >= scp.plim.1 {
                scp.p.clear();
                return Ok(false);
            }
        }
        Ok(true)
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

    fn next_past_b2(&mut self) -> Option<u8> {
        if self.scp[1].p.is_set() {
            self.ks.b2_for_p(self.scp[1].p).ok()
        } else {
            None
        }
    }

    fn get_scp(&mut self) -> &mut Scope<'a> {
        &mut self.scp[if self.scp[1].p.is_set() { 1 } else { 0 }]
    }

    pub fn markcontig<T: MidPos>(&mut self, seq: &mut Iter<u8>) -> Result<u64> {
        self.goffs = 0;
        self.ks.push_contig(self.scp[0].p.pos(), self.goffs);

        loop {
            let mut test;
            if self.scp[1].p.is_set() {
                test = self.ks.b2_for_p(self.scp[1].p).and_then(|b2| {
                    dbg_print!("=> twobit {:x} (past) <=", b2);
                    self.scp[1].complete_and_update_mark(b2, 0)
                });
                // false: cannot complete scope (contig end?) or find next mark after leaving mark.
            } else {
                test = Err(anyhow!("end of contig."));
                while let Some(b2) = seq.next().map(|c| (c >> 1) & 0x7) {
                    dbg_print!("=> twobit {:x} (head) <=", b2);
                    let p = self.scp[0].p;
                    if b2 < 4 {
                        // new sequence is also stored, to enable lookup later.
                        if let Some(qb) = self.ks.b2.get_mut(p.byte_pos()) {
                            *qb |= b2 << (p & 6);
                        }
                        self.ks.p_max = p.pos() + 4;
                        if self.n_stretch > 0 {
                            self.finalize_n_stretch();
                            // If a repetition ends in an N-strectch, thereafter offset to period
                            // may differ or the repetition could be different or entirely gone.
                            self.period = None;
                        } else if self
                            .period
                            .and_then(|d| self.ks.b2_for_p(self.scp[0].p - d).ok())
                            .filter(|&old_b2| old_b2 == b2)
                            .is_some()
                        {
                            //everything that needs to be done.
                            self.scp[0].increment(b2);
                            continue;
                        } else {
                            self.period = None;
                        }
                        // scp funcs also used for scope rebuild, therefore ext is set here.
                        self.scp[0].p.set_extension(0);
                        test = self.scp[0].complete_and_update_mark(b2, 0);
                        break;
                    }
                    if self.scp[0].i != 0 {
                        dbg_print!("started N-stretch at {}.", p);
                        self.goffs += self.scp[0].i as u64;
                        self.ks.push_contig(p, self.goffs);

                        // clear all except orientation and position to rebuild at the start of a new contig.
                        self.scp[0].i = 0;
                        self.n_stretch = 0;
                    }
                    self.n_stretch += 1;
                }
            }
            match test {
                Ok(true) => {}
                Ok(false) => continue,
                Err(e) => {
                    // fixme: use thiserror?
                    if e.to_string() == "end of contig." {
                        break;
                    }
                    dbg_print!("{} (p:{:x})", e, self.get_scp().mark.p);
                    if self.scp[1].p.is_set() {
                        self.scp[1].p.clear();
                    }
                    continue;
                }
            }
            loop {
                dbg_print!(
                    "---------[ past:{:?} ]-------------",
                    self.scp[1].p.is_set()
                );
                let (min_idx, min_p) = (self.get_scp().mark.idx, self.get_scp().mark.p);
                dbg_assert!(min_idx < self.ks.kmp.len(), "{:x}, {:x}", min_idx, min_p);
                let stored_p = self.ks.kmp[min_idx];

                if dbgx!(stored_p.is_replaceable_by(min_p)) {
                    //dbg_print!("[{:#x}] (={:#x}) <= {:#x}", min_idx, stored_p, min_p);
                    // dbg_print!("{}", self.get_scp());
                    if dbgx!(stored_p.is_set_and_not(min_p)) {
                        if !self.rebuild_scope(stored_p)? {
                            self.scp[1].p.clear();
                            break;
                        }
                        self.scp[1].extend_kmer_stack(self.ks)?;
                        //dbgx!(self.ks.kmp[min_idx].set(min_p));
                        dbgx!(self.set_idx_pos(min_idx, min_p));
                        // go back and extend
                    } else {
                        // If already set it was min_p. Then leave dupbit state.
                        if self.ks.kmp[min_idx].is_no_pos() {
                            self.set_idx_pos(min_idx, min_p);
                        }
                        self.scp[1].p.clear();
                        break;
                    }
                } else if dbgx!(stored_p.extension() == min_p.extension()) {
                    // If a kmer occurs multiple times within an extending readlength, only
                    // the first gets a position. During mapping this rule also should apply.
                    if self.recurrance_is_skippable(stored_p, min_idx) {
                        break;
                    }
                } else {
                    dbg_print!("\t\t<!>");
                }
                if dbgx!(!self.get_scp().extend()) {
                    break;
                }
                if self.scp[1].p.is_set() {
                    // includes check for readlength
                    self.scp[1].set_next_mark()?;
                }
            }
            if self.scp[1].p.is_set() && dbgx!(self.scp[1].p.pos() >= self.scp[1].plim.1) {
                // extension requires bases, but we're at contig limit, just leave it.
                self.scp[1].p.clear();
            }
        }
        if self.n_stretch > 0 {
            dbgx!(self.finalize_n_stretch());
        }
        Ok(self.scp[0].p.pos() as u64)
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
        let first_pos = 1 | (kc.kmerlen as u64) << 1;
        let mut seen = 0;
        for i in 1..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert!(ks.kmp[i] == first_pos, "[{:x}], {:x}", i, ks.kmp[i]);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
        Ok(())
    }
    #[test]
    fn test_1n18c1n() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned())?;
        }
        let first_pos = 1 | (kc.kmerlen as u64) << 1;
        let mut seen = 0;
        for i in 1..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert!(ks.kmp[i] == first_pos, "[{:x}], {:x}", i, ks.kmp[i]);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
        Ok(())
    }
    #[test]
    fn test_1n16c() -> Result<()> {
        let kc = KmerConst::new(READLEN, SEQLEN);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCC"[..].to_owned())?;
        }
        let mut seen = 0;
        let first_pos = 1 | (kc.kmerlen as u64) << 1;
        for i in 1..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert!(ks.kmp[i] == first_pos, "[{:x}], {:x}", i, ks.kmp[i]);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
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
        let mut seen = 0;
        for hash in 0..ks_kmp_len {
            let p = kmi.ks.kmp[hash];
            if p.is_set() {
                kmi.rebuild_scope(p)?;
                dbg_assert_eq!(kmi.scp[1].mark.p, p);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 9); // not 100% sure how many this should be
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
                        kmi.rebuild_scope(p)?;
                        dbg_assert_eq!(kmi.scp[1].mark.p, p);
                    }
                }
            }
        }
        Ok(())
    }
}
