extern crate bio;
extern crate flate2;
extern crate arrayvec;
extern crate rustc_serialize;
extern crate bincode;
extern crate clap;
extern crate kmerconst;
extern crate kmerstore;
extern crate extqueue;

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
use self::kmerstore::KmerStore;
use self::kmerconst::KmerConst;
use self::extqueue::ExtQueue;

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

        let readlen = 64;
        let pos_ori_bitcount = 48;
        let c = KmerConst::new(readlen, chrs.iter().map(|x| x.len as usize).sum(), pos_ori_bitcount);
        let mut ks = KmerStore::new(c.bitlen);
        {
            let mut vq = vec![ExtQueue::new(&mut ks, 0, &c, true)];
            let mut kmi = KmerIter::new(&mut vq);
            println!("Chromosome\trunning unique count");
            for chr in &chrs {
                let mut seq = Vec::with_capacity(chr.len as usize);
                idxr.fetch_all(&chr.name).unwrap_or_else(|_| panic!("Error fetching {}.", &chr.name));
                idxr.read(&mut seq).unwrap_or_else(|_| panic!("Error reading {}.", &chr.name));

                let p = kmi.markcontig(&seq, &c);
                println!("{}\t{}", &chr.name, p);
            }
        }
        let mut set: u64 = 0;
        let mut unset: u64 = 0;
        for k in ks.kmp {
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


