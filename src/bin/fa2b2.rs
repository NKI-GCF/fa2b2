extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate flate2;

extern crate fa2b2;

// target/release/fa2b2 -k 16 /net/NGSanalysis/ref/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa test
//
// target/release/fa2b2 /home/roel/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz test

use std::fs::File;
//use std::io::BufReader;
use bio::io::fasta::IndexedReader;
//use std::io::prelude::*;
//use flate2::bufread::MultiGzDecoder;
use bincode::serialize_into;
use clap::{App, Arg};
use std::collections::HashMap;

use anyhow::Result;
use fa2b2::kmerconst::KmerConst;
use fa2b2::kmerloc::PriExtPosOri;
use fa2b2::kmerstore::KmerStore;
use fa2b2::marker::KmerIter;
use std::io::BufWriter;

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
                //.required(true)
                .takes_value(true),
        )
        .get_matches();

    let fa_name = matches.value_of("ref").unwrap();

    let mut idxr = IndexedReader::from_file(&fa_name)
        .unwrap_or_else(|_| panic!("Error opening reference genome"));
    let chrs = idxr.index.sequences();

    let mut out_file = matches
        .value_of("out")
        .map(|f| BufWriter::new(File::create(f).unwrap()));

    let readlen = 64;
    let kc = KmerConst::new(readlen, chrs.iter().map(|x| x.len as usize).sum());
    let mut ks = KmerStore::new(kc.bitlen);
    ks.opt |= 1; //
                 //ks.opt |= 2; // if set also non-priority (ext)kmers
    {
        let mut kmi = KmerIter::new(&mut ks, &kc);
        for chr in &chrs {
            let mut seq = Vec::with_capacity(chr.len as usize);
            idxr.fetch_all(&chr.name)
                .unwrap_or_else(|_| panic!("Error fetching {}.", &chr.name));
            idxr.read(&mut seq)
                .unwrap_or_else(|_| panic!("Error reading {}.", &chr.name));

            kmi.markcontig::<u64>(&chr.name, &mut seq.iter())?;
            //break;
        }
    }
    dump_stats(&ks);

    if let Some(f) = out_file.as_mut() {
        println!("Writing first occurances per kmer to disk");
        serialize_into(f, &ks).unwrap();
    }
    Ok(())
}

fn dump_stats(ks: &KmerStore<u64>) {
    let mut stat = [[0; 0x1_0000]; 2];
    for k in ks.kmp.iter() {
        stat[if k.is_set() { 1 } else { 0 }][k.x()] += 1;
    }
    println!(
        "Unset: {} of {} \t{:.2}%",
        stat[0][0],
        ks.kmp.len(),
        100.0 * stat[0][0] as f64 / ks.kmp.len() as f64
    );
    for j in 1..=0xFFFF {
        if stat[0][j] != 0 {
            println!(
                "Blacklisted for extension {}: {}\t{:.2}%",
                j,
                stat[0][j],
                100.0 * stat[0][j] as f64 / ks.kmp.len() as f64
            );
        }
    }
    for j in 0..=0xFFFF {
        if stat[1][j] != 0 {
            println!(
                "Set for extension {}: {}\t{:.2}%",
                j,
                stat[1][j],
                100.0 * stat[1][j] as f64 / ks.kmp.len() as f64
            );
        }
    }
    println!("{} sections of non-ambiguous code", ks.contig.len());
    println!("{} sections of repetitive code", ks.repeat.len());

    let mut period_counter = HashMap::new();
    for v in ks.repeat.values() {
        let entry = period_counter.entry(v.0).or_insert(0);
        *entry += 1;
    }
    let mut count_vec: Vec<(&u32, &u32)> = period_counter.iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(a.1));
    println!("repetition period and count (top 20 at most)");
    for cv in count_vec.iter().take(20) {
        println!("{}\t{}", cv.0, cv.1);
    }
    println!("..");
}
