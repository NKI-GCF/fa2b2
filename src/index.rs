// (c) Roel Kluin, 2023, GPL v3

use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::mark_contig_threads::MarkContigThreads;
use crate::new_types::extended_position::EXT_MAX;
use crate::rdbg::STAT_DB;
use anyhow::{anyhow, ensure, Result};
use bincode::serialize_into;
use clap::Args;
use crossbeam_channel::TryRecvError;
use noodles_fasta as fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::thread::JoinHandle;

#[derive(Args, Debug)]
pub struct IndexCmd {
    /// The faidx'ed reference genome file, optionally bgzipped
    #[arg(short = 'R', long, value_name = "FASTA", required = true)]
    ref_file: PathBuf,

    /// The output file
    #[arg(short = 'N', long)]
    stats_only: bool,

    // I believe transposons are generally 100 to 10_000 bases in length
    /// Maximum distance between kmers to be considered for repetition
    #[arg(short = 'x', long, default_value = "10000")]
    repetition_max_dist: u32,

    /// Length of sequence reads
    #[arg(short = 'l', long, value_name = "read_length", required = true)]
    read_len: u16,

    /// Total length of genome sequence
    #[arg(short, long, default_value = "3099750718")]
    seq_len: u64,

    /// Seed for indexing (affects x-mer choice)
    #[arg(short = 'S', long, value_name = "seed", default_value = "40164")]
    seed: u16,

    /// Total length of genome sequence
    #[arg(short, long, default_value = "8")]
    ct: usize,
}

pub(crate) fn multi_thread<T>(
    ks: &mut KmerStore,
    kc: KmerConst,
    fa: fasta::Reader<T>,
    ct: usize,
) -> Result<()>
where
    T: BufRead,
{
    ensure!(kc.kmerlen <= 16);
    ensure!(ct == 1 || ct.is_power_of_two());
    let mut threads = MarkContigThreads::new(ct, kc.kmerlen, ks.rep_max_dist);
    eprintln!("Using {ct} threads");

    threads.mark_contig(ks, kc, fa)?;

    // receive the data from threads. They should send in order.
    for nr in 0..ct {
        dbg_print!("blocking receive for thread {}..", nr);
        match threads.rx_in_main.try_recv() {
            Err(TryRecvError::Disconnected) => return Err(anyhow!("Main chnnel is diconnected?")),
            Err(TryRecvError::Empty) => {
                // can happen in tests.
                eprintln!("thread {nr}/{ct} yielded nothing.");
            }
            Ok((kmp, max_extended)) => {
                dbg_print!("Received from thread {}.. len {}", nr, kmp.len());

                //XXX actually why not directly write to disk? do we need ks.kmp still?
                ks.kmp.extend(kmp);
                dbg_print!("Thread {} had {} max extended", nr, max_extended.len());
            }
        }
    }
    // drop so the channels close:
    let inner_threads: Vec<_> = threads.threads.drain(..).collect();
    drop(threads);

    for t in inner_threads {
        t.join().unwrap()?;
    }
    Ok(())
}

pub fn parse_fasta_file(fa: PathBuf) -> Result<fasta::Reader<BufReader<File>>> {
    File::open(fa)
        .map(BufReader::new)
        .map(fasta::Reader::new)
        .map_err(|e| anyhow!("Error opening reference genome: {}", e))
}

pub fn index(cmd: IndexCmd) -> Result<()> {
    let mut ks_file = cmd.ref_file.clone();
    ks_file.set_extension("ks");
    eprintln!("Reading {:?}", cmd.ref_file);
    let fa = parse_fasta_file(cmd.ref_file)?;

    let opt_out = if cmd.stats_only {
        None
    } else {
        ensure!(!Path::new(&ks_file).exists(), "{ks_file:?} already exists!");
        eprintln!("Writing {ks_file:?}");
        Some(BufWriter::new(File::create(ks_file)?))
    };

    let ct = if cmd.ct.is_power_of_two() {
        cmd.ct
    } else {
        cmd.ct.next_power_of_two() - 1
    };

    // Ideally the seed should select against repetitive k-mers as a median.
    // TODO: find out / theorize what seed (if any) does this.
    if cmd.seed != 40164 {
        eprintln!("Warning, changing the seed for indexing and alignment makes your alignment not portable.");
    }

    let kc = KmerConst::new(cmd.seq_len, cmd.read_len, cmd.seed);
    let mut ks = KmerStore::new(kc.bitlen, cmd.repetition_max_dist, cmd.seed)?;

    multi_thread(&mut ks, kc, fa, ct)?;
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
