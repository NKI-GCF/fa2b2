use flate2::{Compression,write::ZlibEncoder,read::ZlibDecoder};
use std::{io,fs::File,cmp::Ordering::{Less,Greater,Equal}};
use bincode::{SizeLimit,rustc_serialize::{decode_from,encode_into}};

use kmerloc::PriExtPosOri;
//use extqueue::ExtQueue;

#[derive(RustcEncodable, RustcDecodable)]
pub struct Contig {
	pub twobit: u64,
	pub genomic: u64,    // if 0: next contig
}

#[derive(RustcEncodable, RustcDecodable)]
pub struct KmerStore<T> {
	pub p_max: u64,
	pub opt: u64,
	pub b2: Vec<u8>,
	pub kmp: Vec<T>,	    // position + strand per k-mer.
	pub contig: Vec<Contig>,
}

impl<T: PriExtPosOri> KmerStore<T>
	where T: rustc_serialize::Encodable + rustc_serialize::Decodable
{
	pub fn new(bitlen: usize) -> Self {
		let shift = bitlen - 2;
		KmerStore {
			p_max: 0,
			opt: 0,
			b2: vec![0; 1 << shift],		// sequence (4 per u8).
			kmp: vec![T::from_u64(0).unwrap(); 1 << (shift + 1)], // kmer positions
			contig: Vec::new(),				// contig info
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
	pub fn push_contig(&mut self, p: u64, goffs: u64) {
		self.contig.push(Contig {
			twobit: p, // TODO: bits could mean something about contig.
			genomic: goffs,
		});
	}
	/// Adjust offset for contig. The
	pub fn offset_contig(&mut self, offset: u64) {
		if let Some(ctg) = self.contig.last_mut() {
			ctg.genomic += offset;
		}
	}
	pub fn get_contig(&self, p: u64) -> usize {
		let mut size = self.contig.len();
		let mut base = 0;
		while size > 0 {
			size /= 2;
			let mid = base + size;
			base = match (self.contig[mid].twobit).cmp(&p) {
				Less	=> mid,
				Greater | Equal => base,
			};
		}
		base
	}
	pub fn get_twobit_after(&self, i: usize) -> Option<u64> {
		if i + 1 < self.contig.len() {Some(self.contig[i+1].twobit)} else {None}
	}
	pub fn get_twobit_before(&self, i: usize) -> Option<u64> {
		if i != 0 {Some(self.contig[i-1].twobit)} else {None}
	}
	pub fn b2_for_p(&self, p: u64) -> Option<u8> {
		self.b2.get(p as usize >> 3).map(|x| (x >> (p & 6)) & 3)
	}
	pub fn is_available(&self, idx: usize, x: usize) -> bool {
		self.kmp[idx].extpos() <= (x as u64) << 48
	}
}


