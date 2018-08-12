extern crate bio;
extern crate flate2;
extern crate arrayvec;
extern crate rustc_serialize;
extern crate bincode;
extern crate clap;
extern crate bit_reverse;

use flate2::{Compression,write::ZlibEncoder,read::ZlibDecoder};

use std::{io,cmp::{self,Ordering::*},fs::File,collections::VecDeque,collections::HashMap};

use bincode::{SizeLimit,rustc_serialize::{decode_from,encode_into}};

#[derive(RustcEncodable, RustcDecodable)]
pub struct Contig {
    twobit: u64,
    genomic: u64,    // if 0: next contig
}

pub struct StatDeq {
    last: (String, String),
    cap: usize,
    d: VecDeque<String>,
    h: HashMap<String, u32>,
}

impl StatDeq {
    fn new(cap: usize) ->  Self {
        StatDeq {
            last: (String::new(), String::new()),
            cap,
            d: VecDeque::with_capacity(cap),
            h: HashMap::new(),
        }
    }
    fn add(&mut self, fmt: String, msg: String) {
        let last_fmt = self.last.0.to_owned();
        let mut last_msg = self.last.1.to_owned();
        if last_fmt == fmt && last_msg == msg {
            let repeat = self.h.entry(String::from("repeat")).or_insert(0);
            *repeat += 1;
        } else {
            if let Some(repeat) = self.h.get_mut(&String::from("repeat")) {
                if self.d.len() == self.cap {
                    self.d.pop_front();
                }
                if *repeat > 0 {
                    last_msg.push_str(&format!(" (repeated {} times)", repeat));
                }
                *repeat = 0;
            }
            self.d.push_back(last_msg);
            *self.h.entry(last_fmt).or_insert(0) += 1;
            self.last = (fmt, msg);
        }
    }
    fn dump(&mut self) -> &'static str {
        self.add(String::new(), String::new());
        println!("Note stack:");
        for msg in self.d.iter() {
            println!("{}", msg);
        }
        println!("--\nTotal file:line:msg(format) and counts:");
        for (msg, ct) in &self.h {
            if msg != "repeat" {
                println!("{}\t{}", msg, ct);
            }
        }
        "done."
    }
}

static mut STAT_DB: Option<StatDeq> = None;

macro_rules! NB {
    ($fmt:expr) => {{
        unsafe {
            if let Some(ref mut db) = STAT_DB {
                db.add(format!("{}:{}:{}", file!(), line!(), $fmt), format!($fmt));
            }
        }
    }};
    ($fmt:expr, $($arg:tt)*) => {{
        unsafe {
            if let Some(ref mut db) = STAT_DB {
                db.add(format!("{}:{}:{}", file!(), line!(), $fmt), format!($fmt, $($arg)*));
            }
        }
    }};
}

macro_rules! NB_assert {
    ($test:expr) => {{
        unsafe {
            if let Some(ref mut db) = STAT_DB {
                assert!($test, db.dump());
            }
        }
    }};
}

macro_rules! NB_assert_eq {
    ($left:expr, $right:expr) => {{
        unsafe {
            if let Some(ref mut db) = STAT_DB {
                assert_eq!($left, $right, db.dump());
            }
        }
    }};
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
                Equal   => panic!("kmer offset at contig boundary!?"),
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
    /*fn copy(&self) -> KmerLoc {
        KmerLoc { idx: self.idx, p: self.p}
    }*/
    fn next(&mut self, step_and_strand: u64) {
        self.p += step_and_strand ^ (self.p & 1);
    }
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

pub struct ExtQueue {
    loc: KmerLoc,
    is_template: bool,
    max_no_kmers: usize,
    kmer: Kmer,
    d: VecDeque<KmerLoc>,
    outlier: Vec<VecDeque<KmerLoc>>,
}

impl ExtQueue {
    pub fn new(p: u64, readlen: usize, kmerlen: usize, is_template: bool) -> Self {
        let max_no_kmers = readlen - kmerlen;
        let mut outlier = vec![VecDeque::with_capacity(max_no_kmers)];
        for x in 1..6 {
            outlier.push(VecDeque::with_capacity(max_no_kmers - (1 << x)));
        }
        ExtQueue {
            loc: KmerLoc::new(0, p),
            is_template,
            max_no_kmers,
            kmer: Kmer::new((1 << (2 * (kmerlen - 1))) as u64),
            d: VecDeque::with_capacity(max_no_kmers),
            outlier, // outlier; either min or max, dependent on extension (resp. 0 or 1+)
        }
    }
    fn next(&mut self, b2: u8) {
        self.loc.next(if self.is_template {
            self.kmer.add_to_seq(u64::from(b2))
        } else {
            self.kmer.add_to_seq(u64::from(b2 ^ 2)) - 1
        });
    }
    fn update(&mut self, x: usize) -> bool {
        let offs = 1 << x;
        if x != 0 && self.d.len() <= offs {
            NB!("ExtQueue is not yet complete");
            return false;
        }
        let loc = self.get_hash_loc(offs);
        if let Some(outlier) = self.outlier.get_mut(x) {
            let is_same = match outlier.front() {
                Some(front) => front.idx == loc.idx && front.p == loc.p,
                None => false,
            };
            if is_same {
                NB_assert!(false);// XXX: this seems like it should not ever occur (but it does).
                outlier.pop_front();
            }
        }
        if x == 0 {
            assert_eq!(self.max_no_kmers, 48);
            if self.all_kmers_complete() {
                self.d.pop_front();
            }
            self.d.push_back(loc);
        }

        if let Some(outlier) = self.outlier.get_mut(x) {
            while !outlier.is_empty() {
                if let Some(back) = outlier.back() {
                    if (x & 1) == 0 {
                        if back.idx < loc.idx || (back.idx == loc.idx && back.p < loc.p) {
                            break;
                        }
                    } else if back.idx > loc.idx || (back.idx == loc.idx && back.p > loc.p) {
                        break;
                    }
                }
                outlier.pop_back();
            }
            outlier.push_back(loc);
        }
        true
    }

    fn is_p(&self, x: usize, p: u64, pos_mask: u64) -> bool {
        match self.outlier[x].front() {
            Some(outlier) => {
                assert!(outlier.p != 0);
                //println!("ol.p:{}, p:{}", (outlier.p & pos_mask) >> 1, (p & pos_mask) >> 1);
                (outlier.p & pos_mask) == (p & pos_mask)
            },
            None => false,
        }
    }

    fn get_hash_loc(&self, offs: usize) -> KmerLoc {
        if offs == 1 { // bij extensie 0: geen hashing.
            //println!("offs==1 (x==0) => returning self.loc with p:{}", self.loc.p >> 1);
            return self.loc;
        }
        let hash = self.loc.idx ^ self.d[offs].idx;
        let half = offs / 2;
        let ori = if hash != 0 {
            if (hash & hash.wrapping_neg()) != 0 {1} else {0}
        } else {
            self.d[half].p & 1
        };
        KmerLoc::new(hash, (self.loc.p - (half << 1) as u64) ^ ((self.loc.p ^ ori) & 1))
    }

    fn all_kmers_complete(&self) -> bool { self.d.len() == self.max_no_kmers }

    fn clear(&mut self) {
        self.d.clear();
        for outlier in &mut self.outlier {
            outlier.clear();
        }
    }
    fn get_idx(&mut self, bufsz: u64) -> u64 {
        let seq = self.kmer.get_strand(self.is_template);

        // flipped if the top bit is set, to reduce size.
        if (seq & bufsz) == 0 {seq} else {(bufsz - 1) & !seq}
    }
}

pub struct KmerIter {
    pos_mask: u64, // : all these are used only once.
    bufsz: u64,
    pos_ori_bitcount: u32, //
    priority_shft: u32, //
    q: u32,
    kmerlen: usize,
    readlen: usize,
    pub ks: KmerStore,
} //^-^\\

impl KmerIter {
    pub fn new(readlen: usize, bitlen: usize, pos_ori_bitcount: u32) -> Self {
        unsafe {
            STAT_DB = Some(StatDeq::new(50));
        }
        // FIXME: different kmer lengths not supported currently. requires kmer extension
        // windows to be adapted)
        assert_eq!(bitlen, 32);
        // to try to make use of the entire u32 domain.

	KmerIter {
            pos_mask: (1 << pos_ori_bitcount) - 1,
            bufsz: 1 << (bitlen - 1),
            pos_ori_bitcount,
            priority_shft: 63,
            q: 0,
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
        println!("contig[{}].twobit:{}, p:{}, add:{}", i, self.ks.contig[i].twobit, p, add);

        if self.ks.contig[i].twobit < p {
            (cmp::max(p - add, self.ks.contig[i].twobit),
                cmp::min(p + add, self.ks.get_twobit_after(i)))
        } else {
            (cmp::max(p - add, self.ks.get_twobit_before(i)),
                cmp::min(p + add, self.ks.contig[i].twobit))
        }
    }
    fn relocate(&mut self, is_template: bool, p: u64) {
        let in_use = p & (1 << self.priority_shft);
        let plimits = self.get_plimits(self.b2pos(p));
        let plim;
        println!("plimits: ({}, {})", self.b2pos(plimits.0), self.b2pos(plimits.1));
        let mut eq = if is_template {
            plim = plimits.1;
            ExtQueue::new(plimits.0, self.readlen, self.kmerlen, true)
        } else {
            plim = plimits.0;
            ExtQueue::new(plimits.1, self.readlen, self.kmerlen, false)
        };
        self.q += 1; // XXX: for debugging
        if self.q == 4 {
            assert!(false);
        }

        // rebuild Extqueue up to where repeat kmer occurs
        let x = (p & !(1 << self.priority_shft)) >> self.pos_ori_bitcount;
        assert!(x == 0); //remove later: relocated for extended kmer is unlikely in beginning
        println!("x:{}", x);
        // !eq.all_kmers_complete() || 
        loop { // XXX: hot loop.
            if is_template {
                assert!(self.pos(p) > self.pos(eq.loc.p));
            } else {
                assert!(self.pos(p) < self.pos(eq.loc.p));
            }
            //println!("p:{}, eq.loc.p:{}", self.pos(p), self.pos(eq.loc.p));
            if let Some(qb) = self.ks.b2.get_mut(eq.loc.p as usize >> 3) {
                let shft = eq.loc.p & 6;
                eq.next((*qb >> shft) & 3);
            }
            if eq.is_p(x as usize, p, self.pos_mask) {
                break;
            }
            let _ = self.progress(&mut eq); // XXX: endless recursion.
        }
        assert!(eq.all_kmers_complete());
        println!("complete");
        assert_eq!(eq.loc.p & in_use, in_use);
        // make sure this position can no longer be used for this extension
        if let Some(e) = self.ks.kmp.get_mut(eq.loc.idx as usize) {
            *e = (eq.loc.p & !self.pos_mask) + (1 << self.pos_ori_bitcount);
        }

        // progress() will extend kmers that were invalidated
        while (eq.loc.p & self.pos_mask & !1) != plim {
            let _ = self.progress(&mut eq);
            if let Some(qb) = self.ks.b2.get_mut(eq.loc.p as usize >> 3) {
                let shft = eq.loc.p & 6;
                eq.next((*qb >> shft) & 3);
            }
        }
    }

    fn update_kmer_entry(&mut self, new: &KmerLoc, x: u64, priority: u64, is_template: bool) -> bool {

        // for a minimum / maximum only: set priority bit.
        let ext_bits = (x << self.pos_ori_bitcount) | priority;
        loop {
            let lookup = *self.ks.kmp.get_mut(new.idx as usize).expect("ks.kmp[idx] out of bounds?");
            let lookup_had_priority = (lookup & (1 << self.priority_shft)) != 0;
            let pos = self.pos(new.p);
            let olpos = self.pos(lookup);

            if lookup == (new.p | ext_bits) {
                NB!("kmer: already set to {}", pos);
                return true;
            }
            let ori = if (lookup & 1) == (new.p & 1) {is_template} else {!is_template};
            if lookup > ext_bits {
                if (lookup & !self.pos_mask) == ext_bits {
                    NB!("extend previous entry (same ext occurance): p:{} <=> {}", olpos, pos);
                    self.relocate(ori, lookup & !1); // XXX: endless recursion.
                    // progress() -> update_kmer_entry() -> relocate() -> progress()
                }
                // if there is one recurring, there may be more, therefore we
                // store the entire eq. TODO: eq by eq extension: problem: orientation.
                // rev_relocate, afh van orientatie van oldp & 1?
                NB!("kmer: refuse overriding higher priority and or extension entry");
                return false; // irreplaceable
            }
            // lookup == ext_bits: marked as unusable for previous, available for this extension.

            if !lookup_had_priority {
                NB!("kmer: accept override non-priority entry by {}", pos);
                break;
            }
            // in use
            if priority == 0 {
                NB!("kmer: refuse overriding priority entry for non-priorty call");
                return false;
            }
            NB!("extend previous entry (this is retaken). p:{} <=> {}", olpos, pos);
            self.relocate(ori, lookup & !1);
        }
        if let Some(e) = self.ks.kmp.get_mut(new.idx as usize) {
            *e = new.p | ext_bits;
            return true;
        }
        false
    }

    // for each added kmer the minimum and extension maxima are updated.
    // returns whether the re was a k-mer not yet in use.
    fn progress(&mut self, eq: &mut ExtQueue) -> bool {
        let mut already_stored = false;

        eq.loc.idx = eq.get_idx(self.bufsz);
        let in_use = 1 << self.priority_shft;

        for x in 0..6 {
            // TODO: also implement non-priority entries (not min/max).
            //
            //if eq.loc.p >= (self.ks.b2.len() << 3) as u64 { // new sequence
            //    // bij relocation zijn posities al geschreven of invalidated.
            //    let _ = self.update_kmer_entry(eq, x as u64, 0, eq.is_template);
            //}
            if !eq.update(x) {
                break;
            }
            if !already_stored {
                if eq.all_kmers_complete() {
                    if let Some(new) = eq.outlier[0].front() {
                        already_stored = self.update_kmer_entry(new, x as u64, in_use, eq.is_template);
                    }
                }
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

