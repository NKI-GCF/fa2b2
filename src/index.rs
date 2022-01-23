use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::marker::KmerIter;
use crate::new_types::extended_position::EXT_MAX;
use anyhow::{anyhow, ensure, Result};
use bincode::serialize_into;
use clap::ArgMatches;
use noodles_fasta as fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;

pub fn index(matches: &ArgMatches) -> Result<()> {
    let fa_name = matches.value_of("ref").unwrap();

    let mut fa = File::open(&fa_name)
        .map(BufReader::new)
        .map(fasta::Reader::new)
        .map_err(|e| anyhow!("Error opening reference genome: {}", e))?;

    let opt_out = if matches.occurrences_of("stats_only") == 1 {
        None
    } else {
        let ks_name = format!("{}.ks", fa_name);
        ensure!(!Path::new(&ks_name).exists(), "{} already exists!", ks_name);
        Some(BufWriter::new(File::create(ks_name)?))
    };

    // I believe transposons are generally 100 to 10_000 bases in length
    let repetition_max_dist = matches
        .value_of("repetition_max_dist")
        .map(|v| v.parse())
        .transpose()?
        .unwrap();

    let read_len = matches
        .value_of("read_length")
        .map(|v| v.parse())
        .transpose()?
        .unwrap();

    let seq_len = matches
        .value_of("sequence_length")
        .map(|v| v.parse())
        .transpose()?
        .unwrap();

    // Ideally the seed should selecting against repetitive k-mers as a median.
    // TODO: find out / theorize what seed may do this.
    if matches.occurrences_of("seed") != 0 {
        eprintln!("Warning, changing the seed for indexing and alignment makes your alignment not portable.");
    }
    let seed = matches
        .value_of("seed")
        .map(|v| v.parse())
        .transpose()?
        .unwrap();

    let kc = KmerConst::new(seq_len, read_len, seed);
    let mut ks = KmerStore::new(kc.bitlen, repetition_max_dist, seed)?;
    let mut kmi = KmerIter::new(&mut ks, &kc);
    for record in fa.records() {
        kmi.markcontig(record?)?;
    }
    make_stats(&ks);

    if let Some(out_file) = opt_out {
        serialize_into(out_file, &ks)?;
    }
    Ok(())
}

fn make_stats(ks: &KmerStore) {
    let mut stat = [[0; EXT_MAX]; 2];
    for k in ks.kmp.iter() {
        stat[if k.is_set() { 1 } else { 0 }][k.x()] += 1;
    }
    println!(
        "Unset: {} of {} (k/x-mers) \t{:.2}%",
        stat[0][0],
        ks.kmp.len(),
        100.0 * stat[0][0] as f64 / ks.kmp.len() as f64
    );
    for j in 1..EXT_MAX {
        if stat[0][j] != 0 {
            println!(
                "Blacklisted for extension {}: {}\t{:.2}%",
                j,
                stat[0][j],
                100.0 * stat[0][j] as f64 / ks.kmp.len() as f64
            );
        }
    }
    let tot_set = ks.kmp.len() - stat[0][0];
    for j in 0..EXT_MAX {
        if stat[1][j] != 0 {
            println!(
                "Set for extension {}: {}\t{:.2}% (of set)",
                j,
                stat[1][j],
                100.0 * stat[1][j] as f64 / tot_set as f64
            );
        }
    }
    println!("{} sections of non-ambiguous code", ks.contig.len());
    println!("{} positions stored for repetitive code", ks.repeat.len());

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
