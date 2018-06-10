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

const DUP: usize = 1 << 63;

#[derive(RustcEncodable, RustcDecodable)]
struct Kmer {
    seq: Vec<u32>,
} //^-^\\

impl Kmer {
    pub fn new() -> Self {
	    Kmer {seq: vec![0; 2]}
    }
    // Add a nucleotide to the sequence. A Nt is added to template(0) in the top two bits.
    pub fn push(&mut self, b2: u32, shift: u32)
    {
        let mask = (1 << shift) - 1;
        self.seq[0] = (self.seq[0] >> 2) | (b2 << shift);
        self.seq[1] = ((self.seq[1] & mask) << 2) | (b2 ^ 2);
    }
    fn get_ori(&mut self) -> (usize)
    {
        // The deviant bit is the first bit that differs between complements. If there
        // is no deviant bit, use bit 1; the sequence is a palindrome in that case.
        let deviant = self.seq[0] ^ self.seq[1];
        let devbit = if deviant != 0 {deviant & deviant.wrapping_neg()} else {1};

        // Orientation is 0 if the deviant bit is set for the template.
        if (self.seq[0] & devbit) == 0 {1} else {0}
    }
}

#[derive(RustcEncodable, RustcDecodable)]
struct Contig {
    twobit: usize,
    genome: usize,    // if 0: next contig
}

#[derive(RustcEncodable, RustcDecodable)]
pub struct KmerMarker {
    shift: u32,
    contig: Vec<Contig>,
    kmp: Vec<usize>,	    // position + strand per k-mer.
} //^-^\\

impl KmerMarker {
    pub fn new(genomesize: usize) -> Self {

        // bit width, required to store all (cumulative) genomic positions, is used as len
        let bitlen = genomesize.next_power_of_two().leading_zeros();

        println!("Using a kmerlength of {}", bitlen / 2);

        KmerMarker {
            shift: bitlen - 2,
            contig: Vec::new(),
            kmp: vec![0; 1 << (bitlen - 1)],
        }
    }
    pub fn markcontig(&mut self, seq: &[u8], mut p: usize) -> (usize)
    {
        let bufsz = self.kmp.len() as u32;
        let shift = (bufsz.leading_zeros() + 1) - 2;
        let bufmask = bufsz | (bufsz - 1);
        assert_eq!(shift, self.shift); // FIXME: if ok, replace self.shift by shift;

        let mut cleared_bits = self.shift;
        let mut correction: i64 = -(p as i64); // position to b2 offset: Ns and contig offset excluded.

	    let mut kmer = Kmer::new();

	    for c in seq {

            let b2: u32;
            match (*c & 0xdf) as char { // mask to include lower case
                'A' => {b2 = 0},
                'C' => {b2 = 2},
                'G' => {b2 = 6},
                'T'|'U' => {b2 = 4},
                _ => {
                    cleared_bits = self.shift;
                    correction += 1;
                    continue;
                }
            }

            // store sequence in twobit (for later - maybe we can do without?)
            p += 2; // shifted for orientation
            kmer.push(b2, self.shift);

            if cleared_bits == 0 {

                // create an orientation-less k-mer index by separating sequence from orientation.
                let ori = kmer.get_ori();
                let comp = kmer.seq[ori];

                // flipped if the top bit is set, to reduce size.
                let olki = if (comp & bufsz) == 0 {comp} else {bufmask & !comp} as usize;

                if let Some(kmp) = self.kmp.get_mut(olki) {
                    *kmp ^= if *kmp == 0 {ori|p} else {DUP};
                }

            } else {

                if cleared_bits == self.shift {

                    // just after N-stretch or start of contig
                    self.contig.push(Contig {twobit: p, genome: (p as i64 + correction) as usize});
                }
                cleared_bits -= 2;

            }
        }
        p
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
