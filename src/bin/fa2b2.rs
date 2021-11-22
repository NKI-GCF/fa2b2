extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate flate2;

extern crate fa2b2;

// target/release/fa2b2 -k 16 /net/NGSanalysis/ref/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa test
//
// target/release/fa2b2 /home/roel/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz test

//use std::fs::File;
//use std::io::BufReader;
//use std::io::prelude::*;
use bio::io::fasta::IndexedReader;
//use flate2::bufread::MultiGzDecoder;
use clap::{App, Arg};

use anyhow::Result;
use fa2b2::kmerconst::KmerConst;
use fa2b2::kmerloc::PriExtPosOri;
use fa2b2::kmerstore::KmerStore;
use fa2b2::marker::KmerIter;
use fa2b2::occurrence::Occurrence;

fn main() -> Result<()> {
    let matches = App::new("fa2b2")
        .version("1.0")
        .about("Read genome fasta file and write 2bit and region data ")
        .arg(
            Arg::with_name("ref")
                .short("r")
                .long("ref")
                .value_name("FASTA")
                .help("The faidx'ed reference genome file, optionally bgzipped")
                .index(1)
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("out")
                .short("o")
                .long("out")
                .value_name("OUTPUT")
                .help("The output file")
                .index(2)
                .required(true)
                .takes_value(true),
        )
        .get_matches();

    let fa_name = matches.value_of("ref").unwrap();

    let mut idxr = IndexedReader::from_file(&fa_name)
        .unwrap_or_else(|_| panic!("Error opening reference genome"));
    let chrs = idxr.index.sequences();

    let outfile = matches.value_of("out").unwrap();

    let readlen = 64;
    let kc = KmerConst::new(readlen, chrs.iter().map(|x| x.len as usize).sum());
    let mut ks = KmerStore::<u64>::new(kc.bitlen);
    ks.opt |= 1; //
                 //ks.opt |= 2; // if set also non-priority (ext)kmers
    {
        let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
        let mut kmi = KmerIter::new(&mut ks, &mut occ);
        println!("Chromosome\trunning unique count");
        for chr in &chrs {
            let mut seq = Vec::with_capacity(chr.len as usize);
            idxr.fetch_all(&chr.name)
                .unwrap_or_else(|_| panic!("Error fetching {}.", &chr.name));
            idxr.read(&mut seq)
                .unwrap_or_else(|_| panic!("Error reading {}.", &chr.name));

            let p = kmi.markcontig::<u64>(&mut seq.iter())?;
            println!("{}\t{}", &chr.name, p);
            //break;
        }
    }
    let mut stat = [[0; 8]; 2];
    for k in ks.kmp {
        stat[if k.no_pos() { 0 } else { 1 }][k.x()] += 1;
    }
    println!("Unset: {}", stat[0][0]);
    for j in 1..8 {
        println!("Blacklisted for extension {}: {}", j, stat[0][j]);
    }
    for j in 0..8 {
        println!("Set for extension {}: {}", j, stat[1][j]);
    }

    println!("Writing first occurances per kmer to disk");
    //ks.serialize().unwrap();
    Ok(())
}
