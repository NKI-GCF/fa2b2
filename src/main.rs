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
use kmer::marker::KmerIter;

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


        let mut idxr = IndexedReader::from_file(&fa_name).unwrap_or_else(|_| panic!("Error opening reference genome"));
        let chrs = idxr.index.sequences();

        let outfile = matches.value_of("out").unwrap();

        let genomesize: usize = chrs.iter().map(|x| x.len as usize).sum();
        // bit width, required to store all (cumulative) genomic positions, is used as len
        println!("Genome size: {}", genomesize);
        println!("power: {}", genomesize.next_power_of_two());

        let bitlen = genomesize.next_power_of_two().trailing_zeros() as usize;
        println!("bitlen: {}", bitlen);
        println!("Using a kmerlength of {}", bitlen / 2);

        let readlen = 64;
        let pos_ori_bitcount = 40;
        let mut kmi = KmerIter::new(readlen, bitlen, pos_ori_bitcount);

        println!("Chromosome\trunning unique count");
        for chr in &chrs {
            let mut seq = Vec::with_capacity(chr.len as usize);
            idxr.fetch_all(&chr.name).expect(&format!("Error fetching {}.", &chr.name));
            idxr.read(&mut seq).expect(&format!("Error reading {}.", &chr.name));

            let p = kmi.markcontig(&seq);
            println!("{}\t{}", &chr.name, p);
        }
        let mut set: u64 = 0;
        let mut unset: u64 = 0;
        for k in kmi.ks.kmp {
            if k == 0 {
                unset += 1;
            } else {
                set += 1;
            }
        }
        println!("set: {}\tunset: {}", set, unset);

        println!("Writing first occurances per kmer to disk");
        //kmi.ks.write(outfile).unwrap();
}


