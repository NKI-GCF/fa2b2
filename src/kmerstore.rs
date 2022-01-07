use ahash::AHashMap;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering::{Equal, Greater, Less};

use crate::kmerloc::PriExtPosOri;
use crate::rdbg::STAT_DB;

use anyhow::{anyhow, ensure, Result};

#[derive(Serialize, Deserialize)]
pub struct Contig {
    pub twobit: u64,
    pub genomic: u64, // if 0: next contig
}

type ContigRng = (u64, u64);
pub type Repeat = (u32, u32);

#[derive(Serialize, Deserialize)]
pub struct KmerStore<T> {
    pub p_max: u64,
    pub opt: u64,
    pub repetition_max_dist: u64,
    pub b2: Vec<u8>,
    pub kmp: Vec<T>, // position + strand per k-mer.
    pub contig: Vec<Contig>,
    pub repeat: AHashMap<u64, Repeat>,
}

impl<T: PriExtPosOri> KmerStore<T> {
    pub fn new(bitlen: usize, repetition_max_dist: u64) -> Self {
        let shift = bitlen - 2;
        KmerStore {
            p_max: 0,
            opt: 0,
            repetition_max_dist: repetition_max_dist.as_pos(),
            b2: vec![0; 1 << shift],                  // sequence (4 per u8).
            kmp: vec![T::no_pos(); 1 << (shift + 1)], // kmer positions
            //kmp,
            contig: Vec::new(), // contig info
            repeat: AHashMap::new(),
        }
    }
    pub fn push_contig(&mut self, p: u64, goffs: u64) {
        self.contig.push(Contig {
            twobit: p.pos(), // TODO: remaining bits could mean something about contig.
            genomic: goffs,
        });
    }
    /// Adjust offset for contig.
    pub fn offset_contig(&mut self, offset: u64) {
        if let Some(ctg) = self.contig.last_mut() {
            ctg.genomic += offset;
        }
    }
    pub fn set_kmp(&mut self, min_idx: usize, min_p: u64) {
        self.kmp[min_idx].set(min_p);
    }
    pub fn get_bitlen(&self) -> usize {
        self.b2.len().trailing_zeros() as usize + 2
    }

    /// binary search contig lower boundary
    fn get_contig(&self, p: u64) -> usize {
        let mut size = self.contig.len();
        let mut base = 0;
        while size > 0 {
            size /= 2;
            let mid = base + size;
            base = match (self.contig[mid].twobit).cmp(&p) {
                Less => mid,
                Greater | Equal => base,
            };
        }
        base
    }
    pub fn get_contig_start_end_for_p(&self, p: u64) -> ContigRng {
        let i = self.get_contig(p);

        dbg_assert!(i < self.contig.len());
        dbg_assert!(self.contig[i].twobit <= p);
        (
            self.get_twobit_before(i).unwrap_or(0),
            self.get_twobit_after(i).unwrap_or(self.p_max),
        )
    }
    fn get_twobit_after(&self, i: usize) -> Result<u64> {
        ensure!(
            i + 1 < self.contig.len(),
            "get_twobit_after(): End of contig"
        );
        Ok(self.contig[i + 1].twobit)
    }
    fn get_twobit_before(&self, i: usize) -> Result<u64> {
        ensure!(i != 0, "get_twobit_before(): Start of contig");
        Ok(self.contig[i - 1].twobit)
    }
    pub fn b2_for_p(&self, p: u64, tag: &str) -> Result<u8> {
        let pos = p.pos();
        ensure!(
            pos < self.p_max || tag == "(R)",
            "running into sequence head"
        );
        self.b2
            .get(p.byte_pos())
            .map(|x| {
                let b2 = (x >> (p & 6)) & 3;
                dbg_print!("{:<30}{}", format!("=> b2 {:x}, p: {:#x}", b2, p), tag);
                b2
            })
            .ok_or_else(|| anyhow!("stored pos past contig? {:#}", p))
    }
    pub fn extend_repetitive(&mut self, min_pos: u64, dist: u32) {
        let repeat = self.repeat.entry(min_pos).or_insert((dist, 0));
        repeat.1 = dist;
    }

    pub fn replace_repetitive(&mut self, old_pos: u64, new_pos: u64, dist: u32) {
        let old_repeat = self.repeat.remove(&old_pos).expect("oldpos not stored");

        let repeat = self.repeat.entry(new_pos).or_insert(old_repeat);
        repeat.1 += dist;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rdbg::STAT_DB;
    #[test]
    fn it_works() {
        //let mut ks = KmerStore::new(32); // allocates 16 gig of data
        let mut ks = KmerStore::<u64>::new(2, 10_000);
        let mut p = 0;
        let goffs = p & !1;

        ks.push_contig(p, goffs);
        ks.offset_contig(10000); // simulate N-stretch of 10000
        p += 64; // 64 Nt's == one readlength
        ks.push_contig(p, goffs);
        ks.offset_contig(64);

        ks.push_contig(p, goffs);
        ks.offset_contig(10000); // another N-stretch of 10000
        ks.push_contig(p, goffs);
        let i = ks.get_contig(11016);
        dbg_assert_eq!(ks.contig[i].twobit, 64);
    }
}
