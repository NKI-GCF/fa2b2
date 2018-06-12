extern crate bio;
extern crate flate2;
extern crate arrayvec;
extern crate rustc_serialize;
extern crate bincode;
extern crate clap;

// target/release/fa2b2 -k 16 /net/NGSanalysis/ref/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa test
//
// target/release/fa2b2 /home/roel/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz test

//use std::fs::File;
//use std::io::BufReader;
//use std::io::prelude::*;
use bio::io::fasta::IndexedReader;
//use flate2::bufread::MultiGzDecoder;
use clap::{App,Arg};

mod kmer;
use kmer::marker::KmerMarker;

fn main() {
    let matches = App::new("fa2b2")
                .version("1.0")
                .about("Read genome fasta file and write 2bit and region data ")
                .arg(Arg::with_name("ref")
                        .short("r")
                        .long("ref")
                        .value_name("FASTA")
                        .help("The faidx'ed reference genome file, optionally bgzipped")
                        .index(1)
                        .required(true)
                        .takes_value(true))
                .arg(Arg::with_name("out")
                        .short("o")
                        .long("out")
                        .value_name("OUTPUT")
                        .help("The output file")
                        .index(2)
                        .required(true)
                        .takes_value(true))
                .get_matches();

        let fa_name =  matches.value_of("ref").unwrap();
        /*let faidx_name = fa_name.to_string() + "fai";

        let fa_file = File::open(fa_name).expect("cannot open fasta");
        let faidx_file = File::open(faidx_name).expect("cannot open fasta index");

        let fa_stream = BufReader::new(fa_file);
        let faidx_stream = BufReader::new(faidx_file);

        let mut gz = MultiGzDecoder::new(fa_stream).expect("cannot open fasta.gz");
        let idxr = IndexedReader::new(gz, faidx_stream).unwrap_or_else(|_| panic!("Error opening reference genome"));

        // Problem: the trait `std::io::Seek` is not implemented for `flate2::bufread::MultiGzDecoder<BufReader<File>>`
        */

        let mut idxr = IndexedReader::from_file(&fa_name).unwrap_or_else(|_| panic!("Error opening reference genome"));
        let chrs = idxr.index.sequences();

        let outfile = matches.value_of("out").unwrap();

        let genomesize: usize = chrs.iter().map(|x| x.len as usize).sum();
        // bit width, required to store all (cumulative) genomic positions, is used as len
        println!("Genome size: {}", genomesize);
        println!("power: {}", genomesize.next_power_of_two());
        let bitlen = genomesize.next_power_of_two().trailing_zeros();
        println!("bitlen: {}", bitlen);
        println!("Using a kmerlength of {}", bitlen / 2);

        let mut kmermarker = KmerMarker::new(bitlen);

        println!("Chromosome\trunning unique count");
        for chr in &chrs {
            let mut seq = Vec::with_capacity(chr.len as usize);
            idxr.fetch_all(&chr.name).expect(&format!("Error fetching {}.", &chr.name));
            idxr.read(&mut seq).expect(&format!("Error reading {}.", &chr.name));

            let p = kmermarker.markcontig(&seq);
            println!("{}\t{}", &chr.name, p);
        }
        let mut dup: u64 = 0;
        let mut uniq: u64 = 0;
        let mut unset: u64 = 0;
        for k in kmermarker.kmp {
            if k == 0 {
                unset += 1;
            } else if (k & (1 << 63)) == 0 {
                uniq += 1;
            } else {
                dup += 1;
            }
        }
        println!("dup: {}\tuniq: {}\tunset: {}", dup, uniq, unset);

        println!("Writing first occurances per kmer to disk");
        //kmermarker.write(outfile).unwrap();
}


