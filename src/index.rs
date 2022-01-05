use crate::kmerconst::KmerConst;
use crate::kmerloc::PriExtPosOri;
use crate::kmerstore::KmerStore;
use crate::marker::KmerIter;

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
        .unwrap_or(10_000);

    let mut len = 0;
    for record in fa.records() {
        len += record?.sequence().len();
    }
    let kc = KmerConst::new(len, repetition_max_dist);
    let mut ks = KmerStore::new(kc.bitlen);
    let mut kmi = KmerIter::new(&mut ks, &kc);
    for record in fa.records() {
        kmi.markcontig::<u64>(record?)?;
    }
    dump_stats(&ks, kc.extent.len());

    if let Some(out_file) = opt_out {
        serialize_into(out_file, &ks)?;
    }
    Ok(())
}

fn dump_stats(ks: &KmerStore<u64>, extent_len: usize) {
    let mut stat = [[0; 0x100]; 2];
    for k in ks.kmp.iter() {
        stat[if k.is_set() { 1 } else { 0 }][k.x()] += 1;
    }
    println!(
        "Unset: {} of {} (k/x-mers) \t{:.2}%",
        stat[0][0],
        ks.kmp.len(),
        100.0 * stat[0][0] as f64 / ks.kmp.len() as f64
    );
    for j in 1..extent_len {
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
    for j in 0..extent_len {
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
