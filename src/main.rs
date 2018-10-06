extern crate bio;
extern crate flate2;
extern crate arrayvec;
extern crate rustc_serialize;
extern crate bincode;
extern crate clap;
extern crate kmerconst;
extern crate kmerstore;
extern crate occurrence;

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
use self::occurrence::Occurrence;

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
        let extent = 48;
        let c = KmerConst::new(readlen, chrs.iter().map(|x| x.len as usize).sum(), extent);
        let mut ks = KmerStore::new(c.bitlen);
        ks.opt |= 1; // 
        //ks.opt |= 2; // if set also non-priority (ext)kmers
        {
            let mut vq = vec![Occurrence::new((0, u64::max_value()), &c, 0, true)];
            let mut kmi = KmerIter::new(&mut ks);
            println!("Chromosome\trunning unique count");
            for chr in &chrs {
                let mut seq = Vec::with_capacity(chr.len as usize);
                idxr.fetch_all(&chr.name).unwrap_or_else(|_| panic!("Error fetching {}.", &chr.name));
                idxr.read(&mut seq).unwrap_or_else(|_| panic!("Error reading {}.", &chr.name));

                let p = kmi.markcontig(&mut vq, &mut seq.iter(), &c);
                println!("{}\t{}", &chr.name, p);
                //break;
            }
        }
        let mut stat = [[0; 8]; 4];
        for k in ks.kmp {
            let priority = k & (1 << c.priority_shft);
            let ext = (k ^ priority) >> c.extent;
            let zeroed = (k & ((1 << c.extent) - 1)) == 0;
            let i = if priority == 0 {
                if zeroed {3} else {1}
            } else {
                if zeroed {2} else {0}
            };
            stat[i][ext as usize] += 1;
        }
        for j in 0..8 {
            println!("Priority position, extension {}: {}", j, stat[0][j]);
        }
        for j in 0..8 {
            println!("Non-priority position, extension {}: {}", j, stat[1][j]);
        }
        for j in 0..8 {
            println!("Priority, zeroed, before extension {}: {}", j, stat[2][j]);
        }
        for j in 0..8 {
            println!("Non-priority, zeroed, extension {}: {}", j, stat[3][j]);
        }


        println!("Writing first occurances per kmer to disk");
        //kmi.ks.write(outfile).unwrap();
}


