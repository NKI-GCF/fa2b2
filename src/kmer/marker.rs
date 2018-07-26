extern crate bio;
extern crate flate2;
extern crate arrayvec;
extern crate rustc_serialize;
extern crate bincode;
extern crate clap;
extern crate bit_reverse;

use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;

use std::io;
use std::fs::File;

use bincode::SizeLimit;
use bincode::rustc_serialize::{decode_from,encode_into};

use std::collections::VecDeque;

use std::cmp::Ordering::*;

#[derive(RustcEncodable, RustcDecodable)]
pub struct Contig {
    twobit: u64,
    genomic: u64,    // if 0: next contig
}


#[derive(Copy, Clone)]
struct Kmer {
    dna: u64,
    rc: u64,
    pow2: u64,
    mask: u64,
} //^-^\\

// determine orientation reduced kmer (ork) for maximum deviant rot 8 within 64bp.
// only store pos for this.

impl Kmer {
    pub fn new(pow2: u64) -> Self {
	    Kmer {
            dna: 0,
            rc: 0,
            pow2,
            mask: pow2 - 1,
        }
    }
    fn add_to_seq(&mut self, b2: u64) -> u64 {
        // Add a nucleotide to the sequence. A Nt is added to template(0) in the top two bits.
        self.dna = (self.dna >> 2) | (b2 * self.pow2);
        self.rc = ((self.rc & self.mask) << 2) ^ 2 ^ b2;

        // 2 is added for next pos; orientation is set in first bit.
        if self.get_ori() {3} else {2}
    }
    // Orientation is true if the deviant bit is set for the template.
    fn get_ori(&self) -> bool {
        let deviant = self.dna ^ self.rc;
        // The first bit that differs between complements, isolated with the '& wrapping_neg()',
        // the devbit, dictates the orientation. If none is set, a palindrome, use bit 1.
        let devbit = if deviant != 0 {deviant & deviant.wrapping_neg()} else {1};
        (self.dna & devbit) != 0
    }
    fn get_strand(&self, is_template: bool) -> u64 {
        if self.get_ori() == is_template {self.dna} else {self.rc}
    }
}

#[derive(Copy, Clone)]
pub struct KmerLoc {
    idx: u64,
    p: u64,
}
impl KmerLoc {
    fn new(idx: u64, p: u64) -> Self {
        KmerLoc { idx, p}
    }
    fn copy(&self) -> KmerLoc {
        KmerLoc { idx: self.idx, p: self.p}
    }
    fn next(&mut self, step_and_strand: u64) {
        self.p += step_and_strand ^ (self.p & 1);
    }
}

pub struct ExtQueue {
    loc: KmerLoc,
    i: usize,
    kmer: Kmer,
    d: Vec<KmerLoc>,
    outlier: Vec<VecDeque<KmerLoc>>,
}

impl ExtQueue {
    pub fn new(p: u64, max_no_kmers: usize, pow2: u64) -> Self {
        let mut outlier = vec![VecDeque::with_capacity(max_no_kmers)];
        for x in 1..6 {
            outlier.push(VecDeque::with_capacity(max_no_kmers - (1 << x)));
        }
        ExtQueue {
            loc: KmerLoc::new(0, p),
            i: 0,
            kmer: Kmer::new(pow2),
            d: Vec::with_capacity(max_no_kmers),
            outlier, // outlier; either min or max, dependent on extension (resp. 0 or 1+)
        }
    }
    fn next(&mut self, is_template: bool, b2: u8) {
        if is_template {
            self.loc.next(self.kmer.add_to_seq(b2 as u64));
        } else {
            self.loc.next(1 - self.kmer.add_to_seq((b2 ^ 2) as u64));
        }
    }
    fn update(&mut self, loc: &KmerLoc, x: usize) {
        if let Some(outlier) = self.outlier.get_mut(x) {
            // XXX: this seems like it may not ever occur.
            let is_same = match outlier.front() {
                Some(front) => front.idx == loc.idx && front.p == loc.p,
                None => false,
            };
            if is_same {
                outlier.pop_front();
            }
        }
        if x == 0 {
            let capacity = self.d.capacity();
            self.d[self.i % capacity] = *loc;
        }

        if let Some(outlier) = self.outlier.get_mut(x) {
            while outlier.len() > 0 {
                if let Some(back) = outlier.back() {
                    if x != 0 {
                        if back.idx < loc.idx || (back.idx == loc.idx && back.p < loc.p) {
                            break;
                        }
                    } else {
                        if back.idx > loc.idx || (back.idx == loc.idx && back.p > loc.p) {
                            break;
                        }
                    }
                }
                outlier.pop_back();
            }
            outlier.push_back(*loc);
        }
    }

    fn is_p(&self, p: u64) -> bool {
        match self.outlier[0].front() {
            Some(outlier) => outlier.p == p,
            None => false,
        }
    }

    fn get_hash_loc(&self, offs: usize) -> KmerLoc {
        if offs == 1 { // bij extensie 0: geen hashing.
            return self.loc;
        }
        let hash = self.loc.idx ^ self.d[(self.i - offs) % self.d.capacity()].idx;
        let half = offs / 2;
        let ori = if hash != 0 {
            if (hash & hash.wrapping_neg()) != 0 {1} else {0}
        } else {
            self.d[(self.i - half) % self.d.capacity()].p & 1
        };
        KmerLoc::new(hash, (self.loc.p - (half << 1) as u64) ^ ((self.loc.p ^ ori) & 1))
    }

    fn all_kmers_complete(&self) -> bool { self.d.len() == self.d.capacity() }

    fn clear(&mut self) {
        self.i = 0;
        self.d.clear();
        for outlier in &mut self.outlier {
            outlier.clear();
        }
    }
    fn get_idx(&mut self, bufsz: u64, is_template: bool) -> u64 {
        let seq = self.kmer.get_strand(is_template);

        // flipped if the top bit is set, to reduce size.
        if (seq & bufsz) == 0 {seq} else {(bufsz - 1) & !seq}
    }
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
    fn new_contig(&mut self, eq: &mut ExtQueue, goffs: u64) {
        if eq.d.len() != 0 {
            eq.clear();
            eq.loc.p &= !1;
            if let Some(genomic) = eq.loc.p.checked_sub(goffs) {
                self.contig.push(Contig { twobit: eq.loc.p, genomic});
            }
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
                Equal   => panic!("kmer offset at contig boundary!?"),
            };
        }
        base
    }
}

pub struct KmerIter {
    pos_mask: u64, // : all these are used only once.
    bufsz: u64,
    pub p: u64,
    pos_ori_bitcount: u32, //
    in_use_shft: u32, //
    max_no_kmers: usize,
    pub ks: KmerStore,
} //^-^\\

impl KmerIter {
    pub fn new(readlen: usize, bitlen: usize, pos_ori_bitcount: u32) -> Self {
        // FIXME: different kmer lengths not supported currently. requires kmer extension
        // windows to be adapted)
        assert_eq!(bitlen, 16);
        // to try to make use of the entire u32 domain.

	    KmerIter {
            pos_mask: (1 << pos_ori_bitcount) - 1,
            bufsz: 1 << (bitlen - 1),
            p: 0,
            pos_ori_bitcount,
            in_use_shft: 63,
            max_no_kmers: readlen - (bitlen / 2),
            ks: KmerStore::new(bitlen),
        }
    }

    fn get_plimit(&mut self, is_template: bool, p: u64) -> u64 {
        // binary search; limit endp to end of contig
        let mut base = self.ks.get_contig(p);
        if is_template {
            let default = p + ((self.max_no_kmers as u64) << 1);
            if self.ks.contig[base].twobit < p {
                base += 1;
                if base == self.ks.contig.len() {
                    let end_of_seq = (self.ks.b2.len() << 3) as u64;
                    return if end_of_seq > default {default} else {end_of_seq};
                }
            }
            let end_of_contig = self.ks.contig[base].twobit;
            if default > end_of_contig {default} else {end_of_contig}
        } else {
            let default = p - ((self.max_no_kmers as u64) << 1);
            if self.ks.contig[base].twobit > p {
                if base == 0 {
                    return if default > 0 {default} else {0};
                }
                base -= 1;
            }
            let start_of_contig = self.ks.contig[base].twobit;
            if start_of_contig < default {default} else {start_of_contig}
        }
    }
    fn relocate(&mut self, is_template: bool, p: u64) {
        let is_new_seq = false;
        let pos = p & self.pos_mask;
        let plimit = self.get_plimit(is_template, pos);
        let mut eq = ExtQueue::new(pos - ((self.max_no_kmers as u64) << 1), self.max_no_kmers, self.bufsz >> 1);

        // rebuild Extqueue up to where repeat kmer occurs
        while eq.all_kmers_complete() == false || !eq.is_p(p) {
            if let Some(qb) = self.ks.b2.get_mut(eq.loc.p as usize >> 3) {
                let shft = eq.loc.p & 6;
                eq.next(is_template, (*qb >> shft) & 3);
            }
            let _ = self.progress(&mut eq, is_template, is_new_seq);
        }
        assert_eq!(eq.loc.p & (1 << self.in_use_shft), 1 << self.in_use_shft);
        // make sure this position can no longer be used for this extension
        if let Some(e) = self.ks.kmp.get_mut(eq.loc.idx as usize) {
            *e = (eq.loc.p & !self.pos_mask) + (1 << self.pos_ori_bitcount);
        }

        // progress() will extend kmers that were invalidated
        while (eq.loc.p & self.pos_mask) != plimit {
            let _ = self.progress(&mut eq, is_template, is_new_seq);
            if let Some(qb) = self.ks.b2.get_mut(eq.loc.p as usize >> 3) {
                let shft = eq.loc.p & 6;
                eq.next(is_template, (*qb >> shft) & 3);
            }
        }
    }

    fn set_extra(&mut self, ent: &KmerLoc, x: u64) {
        // not necessarily a min / max. written to have some fallback.
        let ext_bits = (x+1) << self.pos_ori_bitcount;
        if let Some(e) = self.ks.kmp.get_mut(ent.idx as usize) {
            // set position, if not yet in use or of lower extension
            if *e <= ext_bits {
                *e = ent.p | ext_bits;
            } else if (*e & !self.pos_mask) == ext_bits {
                // duplicate for this extension, invalidate:
                // can only be overwritten by next+ extension, or any with in_use bit set.
                if (*e & !(1 << self.in_use_shft)) != (ent.p | ext_bits) {
                    *e = (x+2) << self.pos_ori_bitcount;
                }
            }
        }
    }

    fn set_in_use(&mut self, eq: &ExtQueue, x: u64) -> bool {
        if let Some(ent) = eq.outlier[0].front() {
            // for a minimum / maximum, set in_use bit.
            let in_use = 1 << self.in_use_shft;
            let ext_bits = ((x+1) << self.pos_ori_bitcount) | in_use;
            let mut p = 0;
            if let Some(e) = self.ks.kmp.get_mut(ent.idx as usize) {
                // overwritable or already set to this pos: ok.
                if *e <= ext_bits || *e == (ent.p | ext_bits) {
                    *e = ent.p | ext_bits;
                    return true;
                }
                p = *e;
            }
            if (p & !self.pos_mask) == ext_bits {
                // duplicate for this extension + in use, invalidate:
                // only available for next+ extension.
                let is_template = (p & 1) == 1;
                self.relocate(is_template, p & !1);
                // if there is one recurring, there may be more, therefore we
                // store the entire eq. TODO: eq by eq extension: problem: orientation.
                // rev_relocate, afh van orientatie van p & 1?
            }
        }
        false
    }

    // for each added kmer the minimum and extension maxima are updated.
    // returns whether the re was a k-mer not yet in use.
    fn progress(&mut self, eq: &mut ExtQueue, is_template: bool, is_new_sequence: bool) -> bool {
        let mut already_stored = false;

        eq.loc.idx = eq.get_idx(self.bufsz, is_template);

        for x in 0..6 {
            let offs = 1 << x;
            if offs > eq.i {
                return false;
            }
            let loc = eq.get_hash_loc(offs);
            if is_new_sequence {
                // bij relocation zijn posities al geschreven of invalidated.
                self.set_extra(&loc, x as u64);
            }
            if eq.all_kmers_complete() {
                eq.update(&loc, x);
                if already_stored == false && self.set_in_use(eq, x as u64) {
                    already_stored = true;
                }
            }
        }
        eq.i += 1;
        if already_stored == false {
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
    pub fn markcontig(&mut self, seq: &[u8]) {

        // position to b2 offset: Ns and contig offset excluded.
        let mut eq = ExtQueue::new(0, self.max_no_kmers, self.bufsz >> 1);
        let goffs = eq.loc.p & !1;
        let is_template = true;
        let is_new_seq = true;

	    for c in seq {
            // store sequence in twobit
            let mut b2 = 0;
            if let Some(qb) = self.ks.b2.get_mut(eq.loc.p as usize >> 3) {
                b2 = (*c >> 1) & 0x7; // convert ascii to 2bit
                if b2 < 4 {
                    let b2_shft = eq.loc.p & 6;
                    eq.next(is_template, b2);
                    *qb |= b2 << b2_shft;
                }
            }
            if b2 < 4 {
                let _ = self.progress(&mut eq, is_template, is_new_seq);
            } else {
                self.ks.new_contig(&mut eq, goffs);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::KmerIter;
    //use super::KmerMarker;

    const BITLEN: u32 = 8;
    const BUFSZ: u32 = 1 << BITLEN - 1;

    #[test]
    fn test_push_b2() {
        let mut kmi = KmerIter::new(BITLEN);
        for i in 0..6 {
            kmi.add_to_seq(1);
        }
        assert_eq!(kmi.kmer.dna, 0x55);
        assert_eq!(kmi.kmer.rc, 0xff);
        assert_eq!(kmi.p, 0xc);

        kmi.add_to_seq(2);
        assert_eq!(kmi.kmer.dna, 0x95);
        assert_eq!(kmi.kmer.rc, 0xfc);
        assert_eq!(kmi.p, 0xf);
        assert_eq!(kmi.get_idx(BUFSZ), 0x6a);

    }
    #[test]
    fn test_push_b2_rc() {
        let mut kmi = KmerIter::new(BITLEN);
        kmi.add_to_seq(0);
        for _ in 0..3 {
            kmi.add_to_seq(3);
        }
        assert_eq!(kmi.p, 0x8);
        assert_eq!(kmi.kmer.dna, 0xfc);
        assert_eq!(kmi.kmer.rc, 0x95);
        assert_eq!(kmi.get_idx(BUFSZ), 0x6a);
    }
    #[test]
    fn test_ori1() {
        let mut kmi = KmerIter::new(BITLEN);
        kmi.kmer.dna = 0xaa;
        assert_eq!(kmi.kmer.get_ori(), true);
    }

    #[test]
    fn test_ori2() {
        let mut kmi = KmerIter::new(BITLEN);
        kmi.kmer.rc = 0xaa;
        assert_eq!(kmi.kmer.get_ori(), false);
    }

    /*
     #[test]
    fn test_shorten_index_flip() {
        let mut index = 0xffffffffffffffaa;
        index = shorten_index(index, BUFSZ);
        assert_eq!(index, 0x55);
    }

    #[test]
    fn test_shorten_index_no_flip() {
        let mut index = 0xffffffffffffff5a;
        index = shorten_index(index, BUFSZ);
        assert_eq!(index, 0x5a);
    }
    */
}

