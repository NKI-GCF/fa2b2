use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::marker::KmerIter;
use crate::new_types::extended_position::{ExtPosEtc, EXT_MAX};
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use crate::xmerhasher::XmerHasher;
use anyhow::{anyhow, ensure, Result};
use bincode::serialize_into;
use clap::ArgMatches;
use crossbeam_channel::unbounded;
use itertools::izip;
use noodles_fasta as fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::{BufReader, BufWriter};
use std::iter::repeat;
use std::path::Path;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;
use std::thread::spawn;

pub(crate) fn multi_thread<T>(
    ks: &mut KmerStore,
    kc: KmerConst,
    mut fa: fasta::Reader<T>,
    no_threads: usize,
) -> Result<()>
where
    T: BufRead,
{
    // channels to send XmerLocs to threads
    let (tx_to_thread, rx_from_main): (Vec<_>, Vec<_>) = repeat(())
        .take(no_threads)
        .map(|_| unbounded::<XmerLoc>())
        .unzip();

    // channels for threads to send XmerLocs to one another when extended past their responsibility
    let (mut tx_inter_thread, rx_inter_thread): (Vec<_>, Vec<_>) = repeat(())
        .take(no_threads)
        .map(|_| unbounded::<XmerLoc>())
        .unzip();

    // rotate transmitter back so that its corresponding receiver is passed to the next thread.
    let first = tx_inter_thread.remove(0);
    tx_inter_thread.push(first);

    // to collect the results
    let (tx_to_main, rx_in_main) = unbounded::<(Vec<ExtPosEtc>, Vec<XmerLoc>)>();
    let shutdown_poll = Arc::new(AtomicUsize::new(0));

    let threads: Vec<_> = izip!(
        0..no_threads,
        rx_from_main.into_iter(),
        tx_inter_thread.into_iter(),
        rx_inter_thread.into_iter(),
    )
    .map(|xmer_channels| {
        let tx_to_main = tx_to_main.clone();
        let shutdown_poll = shutdown_poll.clone();
        let rep_max_dist = ks.rep_max_dist.clone();
        spawn(move || {
            XmerHasher::new(
                no_threads,
                kc.kmerlen,
                rep_max_dist,
                xmer_channels,
                tx_to_main,
            )
            .and_then(|mut xh| xh.work(shutdown_poll))
        })
    })
    .collect();
    drop(shutdown_poll);

    // The main thread to read from fasta and prefilter
    let mut kmi = KmerIter::new(ks, &kc, tx_to_thread)?;

    for res in fa.records() {
        let record = res?;
        dbg_print!("Starting with record {}.", record.name());
        kmi.markcontig(&kc, record)?;
        dbg_print!("Finished with record.");
    }

    dbg_print!("Ending transmission to threads.");
    drop(kmi);

    // receive the data from threads. They should send in order.
    for nr in 0..no_threads {
        dbg_print!("blocking receive for thread {}..", nr);
        let (kmp, max_extended) = rx_in_main.recv()?;
        dbg_print!("Received from thread {}.. len {}", nr, kmp.len());

        //XXX actually why not directly write to disk? do we need ks.kmp still?
        ks.kmp.extend(kmp);
        dbg_print!("Thread {} had {} max extended", nr, max_extended.len());
    }

    for t in threads {
        t.join().unwrap()?;
    }
    Ok(())
}

pub fn index(matches: &ArgMatches) -> Result<()> {
    let fa_name = matches.value_of("ref").unwrap();

    let fa = File::open(&fa_name)
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

    let no_threads = matches
        .value_of("no_threads")
        .map(|v| {
            v.parse().map(|n: usize| {
                if n.is_power_of_two() {
                    n
                } else {
                    n.next_power_of_two() - 1
                }
            })
        })
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

    multi_thread(&mut ks, kc, fa, no_threads)?;
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
