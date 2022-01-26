use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::{
    extended_position::ExtPosEtc,
    extension::Extension,
    position::{PosRange, Position},
};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use crate::xmer_location::XmerLoc;
use anyhow::{ensure, Result};

pub struct PastScope<'a> {
    kc: &'a KmerConst,
    scp: Scope<'a>,
    plim: PosRange,
    period: Position,
}

impl<'a> PastScope<'a> {
    pub(crate) fn new(kc: &'a KmerConst) -> Result<Self> {
        Ok(PastScope {
            kc,
            scp: Scope::new(kc)?,
            plim: Default::default(),
            period: Default::default(),
        })
    }
    fn rebuild(&mut self, ks: &KmerStore, p: ExtPosEtc, idx: usize) -> Result<XmerLoc> {
        let pos = Position::from(p);
        ensure!(pos != Position::default());
        self.plim = ks.get_contig_start_end_for_pos(pos);
        let bound = self.kc.get_kmer_boundaries(pos, self.plim);
        dbg_print!("{}", bound);
        let extension = Extension::from(p);

        let mut p = ExtPosEtc::from((extension, bound.lower()));
        let mut mark = XmerLoc {
            p: ExtPosEtc::from(extension),
            idx: Default::default(),
        };

        //fixme: make this a proper iterator
        for b2 in ks.bit_slice(bound)?.windows(2) {
            if let Some(median_xmer) = self.scp.updated_median_xmer(b2) {
                if p.same_pos_and_ext(median_xmer.p) {
                    mark.p = p;
                    return Ok(mark);
                }
                if median_xmer.get_idx() == idx {
                    dbg_print!(
                        "idx {:x} observed but for {}, not {}",
                        idx,
                        median_xmer.p,
                        p
                    );
                }
                p.incr_pos();
            }
        }
        dbg_print!("kmer {:x} not observed for {} !!", idx, p);
        Ok(mark)
    }
    fn is_on_contig(&self, pos: Position) -> bool {
        self.plim.has_in_range(pos)
    }
}

/*impl<'a> fmt::Display for PastScope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.unshift_pos();
        let mp = self.mark.p.unshift_pos();
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat((p - self.kc.read_len) * 5);
        let r = p - mp;
        if r == 0 {
            let x = self.kc.read_len - n;
            let s = if x != 0 {
                " ".repeat(x << 2) + "|"
            } else {
                String::from("")
            };
            write!(f, "{2}<{3}{: ^1$x}>", self.mark.get_idx(), n << 2, o, s)
        } else if r + self.kc.kmerlen == self.kc.read_len {
            let x = self.kc.read_len - n;
            let s = if x != 0 {
                String::from("|") + &" ".repeat(x << 2)
            } else {
                String::from("")
            };
            write!(f, "{2}<{: ^1$x}{3}>", self.mark.get_idx(), n << 2, o, s)
        } else {
            //let l = self.kc.read_len - r - n;
            //let ls = if o {" ".repeat(o) + "|"} else {String::from("")};
            //let rs = if l {String::from("|") + &" ".repeat(l << 2)} else {String::from("")};
            //write!(f, "{2}<{3}|{: ^1$x}|{4}>", self.mark.idx, n << 2, o, ls, rs)
            write!(f, "")
        }
    }
}*/

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::KmerConst;
    use crate::marker::KmerIter;
    use crate::new_types::position::BasePos;
    use anyhow::Result;
    use noodles_fasta as fasta;
    const SEQLEN: usize = 250;
    const READLEN: u16 = 6;

    #[test]
    fn test_reconstruct1() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        let mut kmi = KmerIter::new(&mut ks, &kc, vec![]);
        let seq_vec = b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..]
            .to_owned();
        let definition = fasta::record::Definition::new("test", None);
        let sequence = fasta::record::Sequence::from(seq_vec);
        kmi.markcontig(&kc, fasta::Record::new(definition, sequence))?;
        let mut scp = PastScope::new(&kc);
        let mut seen = 0;
        for hash in 0..kmi.ks.kmp.len() {
            let p = kmi.ks.kmp[hash];
            if p.is_set() {
                dbg_print!("---[ {:#x} p:{:?} ]---", hash, p);
                let mark = scp.rebuild(&mut kmi.ks, p, hash)?;
                dbg_assert_eq!(mark.p, p.rep_dup_masked(), "[{}]: {:x}", seen, hash);
                seen += 1;
            }
        }
        dbg_assert_eq!(
            seen,
            42,
            "XXX: the number of seen kmers could change, though"
        );
        Ok(())
    }
    #[test]
    fn test_reconstruct_gs4_all() -> Result<()> {
        // all mappable.
        let seqlen: usize = 8;
        let kc = KmerConst::new(seqlen, 2, 0);
        let mut scp = PastScope::new(&kc);

        for gen in 0..=4_usize.pow(seqlen as u32) {
            let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
            let mut kmi = KmerIter::new(&mut ks, &kc, vec![]);
            let seq_vec: Vec<_> = (0..seqlen)
                .map(|i| match (gen >> (i << 1)) & 3 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'T',
                    3 => 'G',
                    _ => unreachable!(),
                })
                .collect();
            dbg_print!("-- k: {} rl: {} {:#x} seq:", kc.kmerlen, kc.read_len, gen);
            dbg_print!("{:?}", seq_vec);

            let vv: Vec<u8> = seq_vec.into_iter().map(|c| c as u8).collect();
            let definition = fasta::record::Definition::new("test", None);
            let sequence = fasta::record::Sequence::from(vv);
            kmi.markcontig(&kc, fasta::Record::new(definition, sequence))?;
            for hash in 0..kmi.ks.kmp.len() {
                let p = kmi.ks.kmp[hash];
                if p.is_set() {
                    dbg_print!("hash: [{:#x}]: p: {:?}", hash, p);
                    let mark = scp.rebuild(&mut kmi.ks, p, hash)?;
                    dbg_assert_eq!(mark.p, p.rep_dup_masked(), "reps: {}", kmi.ks.repeat.len());
                }
            }
        }
        Ok(())
    }
}
