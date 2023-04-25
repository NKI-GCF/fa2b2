// (c) Roel Kluin, 2023, GPL v3

use crate::kmerconst::{KmerConst, XmerHash};
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

// The past scope is a struct currently only in use for tests. It serves to revisit stored twobit
// reference sequence and retrieve the k-mer from that site. Similar to the alignment, except for
// the source. In earlier revisions it was used to retrieve different k-mers from the site, but
// this has become obsolete. It's cheaper to keep on hashing a single chosen k-mer.
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
        let mut pos = ExtPosEtc::from((extension, bound.lower()));
        let mut mark = XmerLoc::default();
        mark.p = ExtPosEtc::from(extension);
        self.scp.reset();

        //fixme: make this a proper iterator
        for b2 in ks.bit_slice(bound)?.chunks(2) {
            //dbg_print!("{b2}");
            if let Some(mut median_xmer) = self.scp.updated_median_xmer(b2) {
                median_xmer.idx = self.kc.hash_and_compress(median_xmer.idx, median_xmer.p.x());
                dbg_print!("{}", median_xmer);
                if pos.same_pos_and_ext(median_xmer.p) {
                    return Ok(median_xmer);
                }
                if median_xmer.get_idx() == idx {
                    dbg_print!(
                        "idx {:x} observed but for {}, not {}",
                        idx,
                        median_xmer.p,
                        pos
                    );
                }
            }
            pos.incr_pos();
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
    use crate::index::multi_thread;
    use crate::kmerconst::KmerConst;
    use crate::marker::KmerIter;
    use crate::new_types::position::BasePos;
    use crate::past_scope::tests::fasta::Reader;
    use anyhow::Result;
    use noodles_fasta as fasta;
    use std::io;
    const GENOME_SIZE: u64 = 250;
    const READLEN: u16 = 6;

    #[test]
    fn test_reconstruct_minimal() -> Result<()> {
        // indexing
        let kc = KmerConst::new(32, 4, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0x2a)?;
        let record = b">test\nAAAAAA";
        multi_thread(&mut ks, kc, Reader::new(io::BufReader::new(&record[..])), 4)?;

        // Testing whether for every stored xmer, if we go back to the sequence at that site,
        // We again retrieve that same xmer.
        let kc = KmerConst::new(32, 4, 0);
        let mut pscp = PastScope::new(&kc)?;
        let p = ExtPosEtc::from_basepos::<u64>(3);
        let mut test = XmerLoc { idx: 0x2a, p: ExtPosEtc::from_basepos::<u64>(3) };
        test.p.set_ori(true);
        let mark = pscp.rebuild(&mut ks, p, test.p.as_u64() as usize)?;

        dbg_assert_eq!(mark, test);
        Ok(())
    }

    #[test]
    fn test_reconstruct_10() -> Result<()> {
        // indexing
        let kc = KmerConst::new(32, 4, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        let record = b">test\nNNNNNNGCGATATTCT";
        multi_thread(&mut ks, kc, Reader::new(io::BufReader::new(&record[..])), 4)?;

        // Testing whether for every stored xmer, if we go back to the sequence at that site,
        // We again retrieve that same xmer.
        let kc = KmerConst::new(32, 4, 0);
        let mut pscp = PastScope::new(&kc)?;
        let mut seen = 0;
        for hash in 0..ks.kmp.len() {
            let p = ks.kmp[hash];
            if p.is_set() {
                dbg_print!("---[ {:#x} p:{:?} ]---", hash, p);
                let mark = pscp.rebuild(&mut ks, p, hash)?;
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
    fn test_reconstruct1() -> Result<()> {
        // indexing
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        let record =
            b">test\nGCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC";
        multi_thread(&mut ks, kc, Reader::new(io::BufReader::new(&record[..])), 4)?;

        // Testing whether for every stored xmer, if we go back to the sequence at that site,
        // We again retrieve that same xmer.
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let mut pscp = PastScope::new(&kc)?;
        let mut seen = 0;
        for hash in 0..ks.kmp.len() {
            let p = ks.kmp[hash];
            if p.is_set() {
                dbg_print!("---[ {:#x} p:{:?} ]---", hash, p);
                let mark = pscp.rebuild(&mut ks, p, hash)?;
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
        let genome_size: u64 = 8;
        let kc = KmerConst::new(genome_size, 2, 0);
        let mut pscp = PastScope::new(&kc)?;

        for gen in 0..=4_usize.pow(genome_size as u32) {
            let kc = KmerConst::new(genome_size, 2, 0);
            let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
            let seq_vec: Vec<_> = (0..genome_size)
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
            let record = [b">test\n", &vv[..]].concat();
            multi_thread(&mut ks, kc, Reader::new(io::BufReader::new(&record[..])), 2)?;
            for hash in 0..ks.kmp.len() {
                let p = ks.kmp[hash];
                if p.is_set() {
                    dbg_print!("hash: [{:#x}]: p: {:?}", hash, p);
                    let mark = pscp.rebuild(&mut ks, p, hash)?;
                    dbg_assert_eq!(mark.p, p.rep_dup_masked(), "reps: {}", ks.repeat.len());
                }
            }
        }
        Ok(())
    }
}
