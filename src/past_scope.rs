use crate::kmer::Kmer;
use crate::kmer::TwoBit;
use crate::kmerconst::KmerConst;
use crate::kmerloc::{ExtPosEtc, KmerLoc};
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::{ensure, Result};
use std::fmt;

pub struct PastScope<'a> {
    kc: &'a KmerConst,
    p: u64,
    i: usize,
    mod_i: usize,
    plim: (u64, u64),
    period: u64,
    mark: KmerLoc,
    d: Vec<Kmer<u64>>, // misschien is deze on the fly uit ks te bepalen?
    z: Vec<usize>,
}

impl<'a> PastScope<'a> {
    pub fn new<T: ExtPosEtc + fmt::LowerHex>(
        ks: &KmerStore,
        kc: &'a KmerConst,
        p: &T,
        idx: usize,
    ) -> Result<Self> {
        let pos = p.pos();
        ensure!(pos != ExtPosEtc::no_pos());
        let plim = ks.get_contig_start_end_for_p(pos);
        let bound = kc.get_kmer_boundaries(pos, plim);
        dbg_print!("{:x}, {:x}", bound.0, bound.1);

        let mut scp = PastScope {
            kc,
            p: bound.0 | p.extension(),
            i: 0,
            mod_i: 0,
            plim,
            period: 0,
            mark: KmerLoc::new(usize::max_value(), p.extension()),
            d: vec![Kmer::new(kc.kmerlen as u32, bound.0); kc.no_kmers],
            z: (0..kc.no_kmers).collect(),
        };
        let x = scp.p.x();
        let bin = kc.get_kmers(x);

        loop {
            if scp.increment(ks.b2_for_p(scp.p, false).unwrap()) {
                // we weten extension op voorhand.
                let base = scp.get_i() - kc.kmerlen;
                if (scp.set_if_optimum(x, base, bin) && scp.all_kmers()) || scp.remark(false)? {
                    if p.same_pos_and_ext(scp.mark.p) {
                        break;
                    }
                    if scp.mark.get_idx() == idx {
                        dbg_print!(
                            "idx {:x} observed but for {:#x}, not {:#x}",
                            idx,
                            scp.mark.p,
                            p
                        );
                    }
                    // XXX ik zou een assertion hier logischer vinden
                    /*assert!(
                        scp.p.pos() < bound.1,
                        "kmer {:x} not observed for {:x} !!",
                        idx,
                        p
                    );*/
                    if scp.p.pos() >= bound.1 {
                        dbg_print!("kmer {:x} not observed for {:x} !!", idx, p);
                        scp.p.clear();
                        break;
                    }
                }
            }
        }
        Ok(scp)
    }
    fn is_before_end_of_contig(&self) -> bool {
        self.p.as_pos() < self.plim.1
    }
    fn is_on_contig(&self, pos: u64) -> bool {
        pos >= self.plim.0 && pos < self.plim.1
    }
}

impl<'a> Scope for PastScope<'a> {
    fn get_mark(&self) -> Option<&KmerLoc> {
        if self.mark.is_set() {
            Some(&self.mark)
        } else {
            None
        }
    }
    fn get_kc(&self) -> &KmerConst {
        self.kc
    }
    fn get_p(&self) -> u64 {
        self.p
    }
    fn get_i(&self) -> usize {
        self.i
    }
    fn get_d(&self, i: usize) -> &Kmer<u64> {
        &self.d[i]
    }
    fn is_repetitive(&self) -> bool {
        self.period != 0
    }
    fn clear_p_extension(&mut self) {
        self.p.clear_extension();
    }

    /// add twobit to k-mers, update k-mer vec, increment pos and update orientation
    /// true if we have at least one kmer.
    fn increment(&mut self, b2: TwoBit) -> bool {
        // XXX: function is hot
        if self.i >= self.kc.kmerlen {
            let old_d = self.d[self.mod_i];
            self.mod_i += 1;
            if self.mod_i == self.kc.no_kmers {
                self.mod_i = 0;
            }
            self.d[self.mod_i] = old_d;
        }
        // first bit is strand bit, set according to kmer orientation bit.
        self.p.set_ori(self.d[self.mod_i].update(b2));
        self.p.incr_pos();
        self.i += 1;
        self.i >= self.kc.kmerlen
    }
    fn extend_p(&mut self) {
        self.p.extend();
    }
    fn mark_reset(&mut self) {
        self.mark.reset();
    }
    fn set_mark(&mut self, idx: usize, p: u64, x: usize) {
        dbg_print!("{:<30}<P>", format!("[{:x}] = {:x} | x({})", idx, p, x));
        self.mark.set(idx, p, x);
    }
    fn set_period(&mut self, period: u64) {
        self.period = period;
    }
}

impl<'a> fmt::Display for PastScope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.unshift_pos() as usize;
        let mp = self.mark.p.unshift_pos() as usize;
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat((p - self.kc.venster) * 5);
        let r = p - mp;
        if r == 0 {
            let x = self.kc.venster - n;
            let s = if x != 0 {
                " ".repeat(x << 2) + "|"
            } else {
                String::from("")
            };
            write!(f, "{2}<{3}{: ^1$x}>", self.mark.get_idx(), n << 2, o, s)
        } else if r + self.kc.kmerlen == self.kc.venster {
            let x = self.kc.venster - n;
            let s = if x != 0 {
                String::from("|") + &" ".repeat(x << 2)
            } else {
                String::from("")
            };
            write!(f, "{2}<{: ^1$x}{3}>", self.mark.get_idx(), n << 2, o, s)
        } else {
            //let l = self.kc.venster - r - n;
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
    use crate::kmerconst::KmerConst;
    use crate::marker::KmerIter;
    use anyhow::Result;
    use noodles_fasta as fasta;
    const SEQLEN: usize = 250;

    #[test]
    fn test_reconstruct1() -> Result<()> {
        let kc = KmerConst::new(SEQLEN);
        let mut ks = KmerStore::new(kc.bitlen, 10_000);
        let mut kmi = KmerIter::new(&mut ks, &kc);
        let seq_vec = b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..]
            .to_owned();
        kmi.ks.p_max = (seq_vec.len() as u64).as_pos();
        let definition = fasta::record::Definition::new("test", None);
        let sequence = fasta::record::Sequence::from(seq_vec);
        kmi.markcontig::<u64>(fasta::Record::new(definition, sequence))?;
        let mut seen = 0;
        for hash in 0..kmi.ks.kmp.len() {
            let p = kmi.ks.kmp[hash];
            if p.is_set() {
                dbg_print!("---[ {:#x} p:{:x} ]---", hash, p);
                let scp = PastScope::new(&kmi.ks, &kc, &p, hash)?;
                dbg_assert_eq!(scp.mark.p, p.rep_dup_masked(), "[{}]: {:x}", seen, hash);
                seen += 1;
            }
        }
        dbg_assert_eq!(
            seen,
            23,
            "XXX: the number of seen kmers could change, though"
        );
        Ok(())
    }
    #[test]
    fn test_reconstruct_gs4_all() -> Result<()> {
        // all mappable.
        let seqlen: usize = 8;
        let kc = KmerConst::new(seqlen);

        for gen in 0..=4_usize.pow(seqlen as u32) {
            let mut ks = KmerStore::new(kc.bitlen, 10_000);
            let mut kmi = KmerIter::new(&mut ks, &kc);
            kmi.ks.p_max = (seqlen as u64).as_pos();
            let seq_vec: Vec<_> = (0..seqlen)
                .map(|i| match (gen >> (i << 1)) & 3 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'T',
                    3 => 'G',
                    _ => unreachable!(),
                })
                .collect();
            dbg_print!("-- k: {} rl: {} {:#x} seq:", kc.kmerlen, kc.venster, gen);
            dbg_print!("{:?}", seq_vec);

            let vv: Vec<u8> = seq_vec.into_iter().map(|c| c as u8).collect();
            let definition = fasta::record::Definition::new("test", None);
            let sequence = fasta::record::Sequence::from(vv);
            kmi.markcontig::<u64>(fasta::Record::new(definition, sequence))?;
            for hash in 0..kmi.ks.kmp.len() {
                let p = kmi.ks.kmp[hash];
                if p.is_set() {
                    dbg_print!("hash: [{:#x}]: p: {:#x}", hash, p);
                    let scp = PastScope::new(&kmi.ks, &kc, &p, hash)?;
                    dbg_assert_eq!(
                        scp.mark.p,
                        p.rep_dup_masked(),
                        "reps: {}",
                        kmi.ks.repeat.len()
                    );
                }
            }
        }
        Ok(())
    }
}
