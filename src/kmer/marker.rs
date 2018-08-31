extern crate bio;
extern crate flate2;
extern crate arrayvec;
extern crate rustc_serialize;
extern crate bincode;
extern crate clap;
extern crate bit_reverse;
extern crate kmer;
extern crate statdeq;

use flate2::{Compression,write::ZlibEncoder,read::ZlibDecoder};

use std::{io,cmp::{self,Ordering::*},fs::File,collections::VecDeque};

use bincode::{SizeLimit,rustc_serialize::{decode_from,encode_into}};

use self::kmer::Kmer;
use self::statdeq::*; // custom logging (NB_*)

#[derive(RustcEncodable, RustcDecodable)]
pub struct Contig {
    twobit: u64,
    genomic: u64,    // if 0: next contig
}

#[derive(RustcEncodable, RustcDecodable)]
pub struct KmerStore {
    contig: Vec<Contig>,
    pub kmp: Vec<u64>,	    // position + strand per k-mer.
    b2: Vec<u8>,
}

impl KmerStore {
    pub fn new(bitlen: usize) -> Self {
        let shift = bitlen - 2;
        let bufsz = 1 << (shift + 1);
        assert_eq!(shift & 1, 0);
        KmerStore {
            contig: Vec::new(),       // contig info
            kmp: vec![0; bufsz],      // kmer positions
            b2: vec![0; 1 << shift],  // sequence (4 per u8).
        }
    }

    pub fn write(&mut self, path: &str) -> io::Result<()> {

	let writer = io::BufWriter::new(File::create(path).unwrap());
        let mut encoder = ZlibEncoder::new(writer, Compression::Fast);
        encode_into(self, &mut encoder, SizeLimit::Infinite).unwrap();

        Ok(())
    }
    pub fn read(&mut self, path: &str) -> io::Result<()> {
        let f = try!(File::open(path));
        let reader = io::BufReader::new(f);

        let mut decoder = ZlibDecoder::new(reader);
        *self = decode_from(&mut decoder, SizeLimit::Infinite).unwrap();

        Ok(())
    }
    fn update_contig(&mut self, eq: &mut ExtQueue, goffs: u64, pos_mask: u64) {
        if !eq.d.is_empty() || (eq.loc.p == 0 && self.contig.is_empty()){
            let p = eq.loc.p & pos_mask;
            self.contig.push(Contig {
                twobit: p, // TODO: bits could mean something about contig.
                genomic: p - goffs,
            });
            eq.clear();
        }
        if let Some(ctg) = self.contig.last_mut() {
            ctg.genomic += 2;
        }
    }
    fn get_contig(&self, p: u64) -> usize {
        let mut size = self.contig.len();
        let mut base = 0;
        while size > 0 {
            size /= 2;
            let mid = base + size;
            base = match (self.contig[mid].twobit).cmp(&p) {
                Less    => mid,
                Greater => base,
                Equal   => {
                    NB_panic!("kmer offset at contig boundary!?");
                    0
                },
            };
        }
        base
    }
    fn get_twobit_after(&self, i: usize) -> u64 {
        if i + 1 < self.contig.len() {
            self.contig[i+1].twobit
        } else {
            (self.b2.len() << 3) as u64
        }
    }
    fn get_twobit_before(&self, i: usize) -> u64 {
        if i != 0 {self.contig[i+1].twobit} else {0}
    }
    fn b2_for_p(&self, p: u64) -> u8 {
        (self.b2[p as usize >> 3] >> (p & 6)) & 3
    }
}

#[derive(Copy, Clone)]
pub struct KmerLoc {
    idx: usize,
    p: u64,
}
impl KmerLoc {
    fn new(idx: usize, p: u64) -> Self {
        KmerLoc {idx, p}
    }
    /*fn copy(&self) -> KmerLoc {
        KmerLoc { idx: self.idx, p: self.p}
    }*/
    fn next(&mut self, step_and_strand: u64) {
        self.p += step_and_strand ^ (self.p & 1);
    }
}

fn test_pop(marked: &VecDeque<KmerLoc>, parti: &KmerLoc) -> bool{
    let front = marked.front().unwrap();
    front.idx == parti.idx && front.p == parti.p
}
pub struct ExtQueue {
    loc: KmerLoc,
    is_template: bool,
    max_no_kmers: usize,
    ext_max: usize,
    kmer: Kmer<u64>,
    d: VecDeque<KmerLoc>,
    marked: Vec<VecDeque<KmerLoc>>,
}

impl ExtQueue {
    pub fn new(p: u64, readlen: usize, kmerlen: usize, is_template: bool) -> Self {
        let max_no_kmers = readlen - kmerlen;
        let mut marked = vec![VecDeque::with_capacity(max_no_kmers)];
        let mut ext_max = 1;
        while max_no_kmers < (1 << ext_max) {
            marked.push(VecDeque::with_capacity(max_no_kmers - (1 << ext_max)));
            ext_max += 1;
        }
        ExtQueue {
            loc: KmerLoc::new(0, p),
            is_template,
            max_no_kmers,
            ext_max,
            kmer: Kmer::new(kmerlen as u32),
            d: VecDeque::with_capacity(max_no_kmers),
            marked, // marked; either min or max, dependent on extension (resp. 0 or 1+)
        }
    }
    fn next(&mut self, b2: u8) {
        let flip_ori = if self.is_template {self.kmer.add(b2)} else { ! self.kmer.add(b2 ^ 2)};
        // 2 is added for next pos; orientation is set in first bit.
        self.loc.next(if flip_ori {3} else {2});
    }
    fn reestablish_marked(&mut self, x: usize, loc: KmerLoc) {
        let marked = &mut self.marked[x];
        loop {
            match marked.back() {
                None => break,
                Some(back) => {
                    if (x & 1) == 0 {
                        if back.idx < loc.idx || (back.idx == loc.idx && back.p < loc.p) {
                            break;
                        }
                    } else if back.idx > loc.idx || (back.idx == loc.idx && back.p > loc.p) {
                        break;
                    }
                },
            }
            marked.pop_back();
        }
        marked.push_back(loc);
    }

    fn update(&mut self, x: usize, leaving: Option<KmerLoc>) -> bool {
        let offs = 1 << x;
        if self.d.len() <= offs {
            if self.d.len() == offs {NB!("ExtQueue is complete [{}]", offs);}
            return false;
        }
        let loc = self.get_hash_loc(&self.loc, offs);
        self.reestablish_marked(x, loc);

        if let Some(parti) = leaving {
            let hash = self.get_hash_loc(&parti, offs);
            let marked = &mut self.marked[x];
            if test_pop(&marked, &hash) {
                marked.pop_front();
            }
        }
        true
    }

    fn is_p(&self, x: usize, p: u64, pos_mask: u64) -> bool {
        match self.marked[x].front() {
            Some(marked) => {
                assert!(marked.p != 0);
                NB!("ol.p:{}, p:{}", (marked.p & pos_mask) >> 1, (p & pos_mask) >> 1);
                (p & pos_mask & !1) == (marked.p & pos_mask & !1)
            },
            None => false,
        }
    }

    fn get_hash_loc(&self, loc: &KmerLoc, offs: usize) -> KmerLoc {
        let hash = loc.idx ^ self.d[offs].idx;
        let half = offs / 2;
        let ori = if hash != 0 {
            if (hash & hash.wrapping_neg()) != 0 {1} else {0}
        } else {
            self.d[half].p & 1
        };
        KmerLoc::new(hash, (loc.p - (half << 1) as u64) ^ ((loc.p ^ ori) & 1))
    }

    fn all_kmers_complete(&self) -> bool { self.d.len() == self.max_no_kmers }

    fn clear(&mut self) {
        self.d.clear();
        for marked in &mut self.marked {
            marked.clear();
        }
    }
}

pub struct KmerIter {
    pos_mask: u64, // : all these are used only once.
    pos_ori_bitcount: u32, //
    priority_shft: u32, //
    kmerlen: usize,
    readlen: usize,
    pub ks: KmerStore,
} //^-^\\

impl KmerIter {
    pub fn new(readlen: usize, bitlen: usize, pos_ori_bitcount: u32) -> Self {
        NB_new!(500);

	KmerIter {
            pos_mask: (1 << pos_ori_bitcount) - 1,
            pos_ori_bitcount,
            priority_shft: 63,
            kmerlen: bitlen / 2,
            readlen,
            ks: KmerStore::new(bitlen),
        }
    }
    fn pos(&self, p: u64) -> u64 {(p & self.pos_mask) >> 1}
    fn b2pos(&self, p: u64) -> u64 {p & self.pos_mask & !1}

    fn get_plimits(&self, p: u64) -> (u64, u64) {

        // binary search; limit endp to end of contig
        let i = self.ks.get_contig(p);
        let add = (self.readlen << 1) as u64;
        println!("contig[{}].twobit:{}, p:{}, add:{}", i, self.ks.contig[i].twobit, p, add); //XXX:1

        if self.ks.contig[i].twobit < p {
            let x = self.ks.get_twobit_before(i);
            (if p > x + add {p - add} else {x}, cmp::min(p + add, self.ks.get_twobit_after(i)))
        } else {
            let x = self.ks.get_twobit_before(i);
            (if p > x + add {p - add} else {x}, cmp::min(p + add, self.ks.contig[i].twobit))
        }
    }
    fn relocate(&mut self, is_template: bool, p: u64) {
        let in_use = p & (1 << self.priority_shft);
        let plimits = self.get_plimits(self.b2pos(p));
        let plim;
        NB!("plimits: ({}, {})", self.b2pos(plimits.0), self.b2pos(plimits.1));
        let mut eq = if is_template {
            plim = plimits.1;
            ExtQueue::new(plimits.0, self.readlen, self.kmerlen, true)
        } else {
            plim = plimits.0;
            ExtQueue::new(plimits.1, self.readlen, self.kmerlen, false)
        };

        // rebuild Extqueue up to where repeat kmer occurs
        let x = (p & !(1 << self.priority_shft)) >> self.pos_ori_bitcount;
        NB!("relocate x:{}", x);
        NB_assert!(x == 0); //remove later: relocated for extended kmer is unlikely in beginning
        let mut already_adapted = false;

        for i in 0..self.readlen {
            //println!("p:{}, eq.loc.p:{}", self.pos(p), self.pos(eq.loc.p));
            let b2 = self.ks.b2_for_p(eq.loc.p);
            eq.next(b2);
            if i >= self.kmerlen {
                if !already_adapted && eq.is_p(x as usize, p, self.pos_mask) {
                    // XXX waarom groeit extensie niet ??
                    let n = (eq.loc.p & !(1 << self.priority_shft)) >> self.pos_ori_bitcount;
                    NB!("0x{:08x}: required extension incremented to {}", eq.loc.idx, n + 1);
                    self.ks.kmp[eq.loc.idx] = in_use | ((n + 1) << self.pos_ori_bitcount);
                    already_adapted = true;
                }
                let _ = self.progress(&mut eq);
            }
        }
        NB_assert!(eq.all_kmers_complete());
        println!("complete"); //XXX:2
        NB_panic!("FIXME: will stack overflow after here");
//        NB_assert!((eq.loc.p & in_use) == in_use);
        // make sure this position can no longer be used for this extension
        self.ks.kmp[eq.loc.idx] = (eq.loc.p & !self.pos_mask) + (x << self.pos_ori_bitcount);

        // progress() will extend kmers that were invalidated
        while (eq.loc.p & self.pos_mask & !1) != plim {
            let _ = self.progress(&mut eq);
            let b2 = self.ks.b2_for_p(eq.loc.p);
            eq.next(b2);
        }
    }

    fn update_kmer_entry(&mut self, new: &KmerLoc, x: u64, priority: u64, is_template: bool) -> bool {

        // for a minimum / maximum only: set priority bit.
        let ext_bits = (x << self.pos_ori_bitcount) | priority;
        loop {
            let former = self.ks.kmp[new.idx];
            let former_had_priority = (former & (1 << self.priority_shft)) != 0;
            let pos = self.pos(new.p);
            let olpos = self.pos(former);

            if former == (new.p | ext_bits) {
                //NB!("0x{:08x}: already set to {}", new.idx, pos);
                return true;
            }
            let ori = if (former & 1) == (new.p & 1) {is_template} else {!is_template};
            if former > ext_bits {
                if (former & !self.pos_mask) == ext_bits {
                    NB!("0x{:08x}: previous entry has same ext({}); p:{} <=> {}",
                        new.idx, x, olpos, pos);
                    self.relocate(ori, former & !1);
                    // progress() -> update_kmer_entry() -> relocate() -> progress()
                }
                // if there is one recurring, there may be more, therefore we
                // store the entire eq. TODO: eq by eq extension: problem: orientation.
                // rev_relocate, afh van orientatie van oldp & 1?
                NB!("0x{:08x}: refuse overriding higher priority | extension entry", new.idx);
                return false; // irreplaceable
            }
            // former == ext_bits: marked as unusable for previous, available for this extension.

            if !former_had_priority {
                if former == 0 {
                    NB!("0x{:08x}: new address for {}", new.idx, pos);
                } else {
                    NB!("0x{:08x}: accept override non-priority entry by {}", new.idx, pos);
                }
                break;
            }
            // in use
            if priority == 0 {
                NB!("0x{:08x}: refuse overriding priority entry for non-priorty call", new.idx);
                return false;
            }
            NB!("extend previous entry (this is retaken). p:{} <=> {}", olpos, pos);
            self.relocate(ori, former & !1);
        }
        if let Some(e) = self.ks.kmp.get_mut(new.idx) {
            *e = new.p | ext_bits;
            return true;
        }
        false
    }

    // for each added kmer the minimum and extension maxima are updated.
    // returns whether the re was a k-mer not yet in use.
    fn progress(&mut self, eq: &mut ExtQueue) -> bool {

        eq.loc.idx = eq.kmer.get_idx(true);
        let in_use = 1 << self.priority_shft;

        let leaving = if eq.all_kmers_complete() {
            let parti = eq.d.pop_front().unwrap();
            let marked = &mut eq.marked[0];
            if test_pop(&marked, &parti) {
                marked.pop_front();
            }
            Some(parti)
        } else {
            None
        };
        eq.d.push_back(eq.loc);
        let loc = eq.loc;
        eq.reestablish_marked(0, loc);
        let mut already_stored = self.update_kmer_entry(
            eq.marked[0].front().unwrap(), 0, in_use, eq.is_template);

        for x in 1..eq.ext_max {
            // TODO: also implement non-priority entries (not min/max).
            //
            //if eq.loc.p >= (self.ks.b2.len() << 3) as u64 { // new sequence
            //    // bij relocation zijn posities al geschreven of invalidated.
            //    let _ = self.update_kmer_entry(eq, x as u64, 0, eq.is_template);
            //}
            if !eq.update(x, leaving) {
                break;
            }
            if !already_stored && eq.all_kmers_complete() {
                let new = eq.marked[0].front().unwrap();
                already_stored = self.update_kmer_entry(new, x as u64, in_use, eq.is_template);
            }
        }
        if !already_stored {
            // TODO:
            // If still not unique, include all kmers wrapping sum hashes for minima within
            // 64bp, until 1 << 16; to enable paired-end matching.
            //
            // Maybe, not every, but only for neighbouring unmappable sites?
        }
        already_stored
    }

    // with 16 Nt we can make an entry, but minima / maxima are determined only with
    // 64 Nts. Minima for kmers with no extension; for extensions 2, 4, 8, 16, 32 and 48
    // a hashed kmer maximum is used as index, instead. (with 48 there is one entry).
    // If still not unique, include all kmers wrapping sum hashes for minima within
    // 64bp, until 1 << 16 This is to enable paired-end matching.
    pub fn markcontig(&mut self, seq: &[u8]) -> u64 {

        // position to b2 offset: Ns and contig offset excluded.
        let mut eq = ExtQueue::new(0, self.readlen, self.kmerlen, true);
        let goffs = eq.loc.p & !1;
        let mut i = 0;

	for c in seq {
            // store sequence in twobit
            let mut b2 = 0;
            if let Some(qb) = self.ks.b2.get_mut(eq.loc.p as usize >> 3) {
                b2 = (*c >> 1) & 0x7; // convert ascii to 2bit
                if b2 < 4 {
                    let b2_shft = eq.loc.p & 6;
                    eq.next(b2); // adds to dna/rc of kmer and increments location
                    *qb |= b2 << b2_shft;
                }
            }
            if b2 < 4 {
                //println!("i:{}", i);
                if i >= self.kmerlen {
                    let _ = self.progress(&mut eq);
                }
                i += 1;
            } else {
                self.ks.update_contig(&mut eq, goffs, self.pos_mask ^ 1);
                i = 0;
            }
        }
        eq.loc.p
    }
}

#[cfg(test)]
mod tests {
    use super::ExtQueue;
    use super::KmerIter;
    const READLEN: usize = 16;
    const KMERLEN: usize = 4;
    const BITLEN: usize = KMERLEN * 2;
    const BUFSZ: u64 = 1 << BITLEN - 1;

    #[test]
    fn test_kmi() {
        let mut kmi = KmerIter::new(READLEN, BITLEN, 48);
        let seq: Vec<u8> = b"CCCCCC"[..].to_owned();
        let p = kmi.markcontig(&b"CCCCCC"[..].to_owned());
        assert_eq!(p, 0xc);
    }


    #[test]
    fn test_push_b2() {
        let mut eq = ExtQueue::new(0, READLEN, KMERLEN, true);
        for _ in 0..6 {
            eq.next(1);
        }
        assert_eq!(eq.kmer.dna, 0x55);
        assert_eq!(eq.kmer.rc, 0xff);
        assert_eq!(eq.loc.p, 0xc);

        eq.next(2);
        assert_eq!(eq.kmer.dna, 0x95);
        assert_eq!(eq.kmer.rc, 0xfc);
        assert_eq!(eq.loc.p, 0xf);
        assert_eq!(eq.get_idx(BUFSZ), 0x6a);

    }
    #[test]
    fn test_push_b2_rc() {
        let mut eq = ExtQueue::new(0, READLEN, KMERLEN, true);
        eq.next(0);
        for _ in 0..3 {
            eq.next(3);
        }
        assert_eq!(eq.loc.p, 0x8);
        assert_eq!(eq.kmer.dna, 0xfc);
        assert_eq!(eq.kmer.rc, 0x95);
        assert_eq!(eq.get_idx(BUFSZ), 0x6a);
    }
    #[test]
    fn test_ori1() {
        let mut eq = ExtQueue::new(0, READLEN, KMERLEN, true);
        eq.kmer.dna = 0xaa;
        assert_eq!(eq.kmer.get_ori(), true);
    }

    #[test]
    fn test_ori2() {
        let mut eq = ExtQueue::new(0, READLEN, KMERLEN, true);
        eq.kmer.rc = 0xaa;
        assert_eq!(eq.kmer.get_ori(), false);
    }
}

