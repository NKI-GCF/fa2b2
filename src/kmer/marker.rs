extern crate bio;
extern crate flate2;
extern crate arrayvec;
extern crate rustc_serialize;
extern crate bincode;
extern crate clap;

use flate2::Compression;
use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;

use std::io;
use std::fs::File;

use bincode::SizeLimit;
use bincode::rustc_serialize::{decode_from,encode_into};

const DUP: u64 = 1 << 63;

#[derive(RustcEncodable, RustcDecodable)]
struct Contig {
    twobit: u64,
    genome: i64,    // if 0: next contig
}

#[derive(RustcEncodable, RustcDecodable)]
struct Kmer {
    dna: u64,
    rc: u64,
    mask: u64,
    correction: i64, // position to b2 offset: Ns and contig offset excluded.
    p: u64,
    bufsz: usize,
    shift: u32,
    pending: u32
} //^-^\\

impl Kmer {
    pub fn new(bufsz: usize, p: u64) -> Self {
        let shift = bufsz.trailing_zeros() - 1;
        assert_eq!(shift & 1, 0);
	    Kmer {
            dna: 0,
            rc: 0,
            mask: (1 << shift) - 1,
            correction: -(p as i64),
            p: p,
            bufsz: bufsz,
            shift: shift,
            pending: shift,
        }
    }
    fn add_to_seq(&mut self, b2: u64) -> u8 {
        // Add a nucleotide to the sequence. A Nt is added to template(0) in the top two bits.
        self.dna = (self.dna >> 2) | (b2 << self.shift);
        self.rc = ((self.rc & self.mask) << 2) ^ 2 ^ b2;

        let quad_shift = self.p & 6;
        self.p += if self.get_ori() != 0 {3} else {2} - (self.p & 1);
        (b2 as u8) << quad_shift
    }
    // Orientation is 0 if the deviant bit is not set for the template.
    fn get_ori(&self) -> u64 {
        // The deviant bit is the first bit that differs between complements. If there
        // is no deviant bit, use bit 1; the sequence is a palindrome in that case.
        let deviant = self.dna ^ self.rc;
        if deviant != 0 {self.dna & deviant & deviant.wrapping_neg()} else {self.dna & 1}
    }
    fn get_index(&self) ->usize {
        // create an orientation-less k-mer index by separating sequence from orientation.
        let comp = if (self.p & 1) == 0 {self.dna} else {self.rc} as usize;

        // flipped if the top bit is set, to reduce size.
        if (comp & self.bufsz) == 0 {comp} else {(self.bufsz - 1) & !comp}
    }
    fn ambiguous_seq(&mut self) -> u8 {
        self.pending = self.shift;
        self.correction += 2;
        0
    }
    fn is_next_contig(&mut self) -> bool {
        let ret = self.pending == self.shift;
        self.pending -= 2;
        ret
    }
    fn get_contig(&self) -> Contig {
        Contig {twobit: self.p, genome: (self.p as i64 + self.correction) >> 1}
    }
}

#[derive(RustcEncodable, RustcDecodable)]
pub struct KmerMarker {
    contig: Vec<Contig>,
    pub kmp: Vec<u64>,	    // position + strand per k-mer.
    b2: Vec<u8>,
} //^-^\\

impl KmerMarker {
    pub fn new(bitlen: u32) -> Self {
        KmerMarker {
            contig: Vec::new(),
            kmp: vec![0; 1 << bitlen - 1],
            b2: vec![0; 1 << bitlen - 2],
        }
    }
    pub fn markcontig(&mut self, seq: &[u8]) -> u64 {

        let mut kmer = Kmer::new(self.kmp.len(), match self.contig.len() {0 => 0, n => self.contig[n-1].twobit});

	    for c in seq {
            // store sequence in twobit
            if let Some(quadb2) = self.b2.get_mut(kmer.p as usize >> 3) {

                *quadb2 |= match (*c & 0xdf) as char { // mask to include lower case
                    'A' => {kmer.add_to_seq(0)},
                    'C' => {kmer.add_to_seq(2)},
                    'G' => {kmer.add_to_seq(6)},
                    'T'|'U' => {kmer.add_to_seq(4)},
                    _ => {kmer.ambiguous_seq()}
                };
            }
            if kmer.pending == 0 {

                if let Some(kmp) = self.kmp.get_mut(kmer.get_index()) {
                    *kmp |= if *kmp == 0 {kmer.p} else {DUP};
                }

            } else if kmer.is_next_contig() {

                self.contig.push(kmer.get_contig());
            }
        }
        kmer.p
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
}

#[cfg(test)]
mod tests {
    use super::Kmer;
    //use super::KmerMarker;

    const BITLEN: u32 = 4 << 1;
    const BUFSZ: usize = 1 << BITLEN - 1;

    #[test]
    fn test_push_b2() {
        /*let kmermarker = KmerMarker {
            contig: Vec::new(),
            kmp: Vec::new(),
            b2: Vec::new(),
        };*/
        let mut kmer = Kmer::new(BUFSZ, 0);
        kmer.dna = 0x55;
        kmer.rc = 0xff;
        kmer.p = 4;

        let quad_b2 = kmer.add_to_seq(2);
        assert_eq!(quad_b2, 0x20);
        assert_eq!(kmer.dna, 0x95);
        assert_eq!(kmer.rc, 0xfc);
        assert_eq!(kmer.p, 0x7);
        assert_eq!(kmer.get_index(), 0x3);
    }

    #[test]
    fn test_ori1() {
        let mut kmer = Kmer::new(BUFSZ, 0);
        kmer.dna = 0xaa;
        assert_eq!(kmer.get_ori(), 2);
    }

    #[test]
    fn test_ori2() {
        let mut kmer = Kmer::new(BUFSZ, 0);
        kmer.rc = 0xaa;
        assert_eq!(kmer.get_ori(), 0);
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

