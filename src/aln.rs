use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::scope::Scope;
use crate::xmer_location::XmerLoc;
use anyhow::{ensure, Result};
use bincode::deserialize_from;
use clap::ArgMatches;
use noodles_fastq as fastq;
use std::{fs, io, path};

pub fn aln(matches: &ArgMatches) -> Result<()> {
    let fa_name = matches.value_of("ref").unwrap();
    let fq = matches.value_of("fastq1").unwrap();
    let mut reader = fs::File::open(fq)
        .map(io::BufReader::new)
        .map(fastq::Reader::new)?;

    let read_len = matches
        .value_of("read_length")
        .map(|v| v.parse())
        .transpose()?
        .unwrap();

    let ks_name = format!("{}.ks", fa_name);
    ensure!(
        path::Path::new(&ks_name).exists(),
        "{} does not exist!",
        ks_name
    );
    let ks_file = io::BufReader::new(fs::File::create(ks_name)?);
    let ks: KmerStore = deserialize_from(ks_file)?;
    let bitlen = ks.get_bitlen();

    let kc = KmerConst::from_bitlen(bitlen, read_len, ks.seed);

    let scp = Scope::new(&kc);
    let mut aln = Aligner::new(ks, scp);

    aln.align(&mut reader)?;
    Ok(())
}

struct Aligner<'a> {
    ks: KmerStore,
    scp: Scope<'a>,
}

impl<'a> Aligner<'a> {
    fn new(ks: KmerStore, scp: Scope<'a>) -> Self {
        Aligner { ks, scp }
    }
    /// FIXME: pass through the best XmerLocs for mapping - i.e. matching in ks.kmp
    fn filter_median_xmers(&mut self, b: u8) -> Option<XmerLoc> {
        None
    }
    fn align(&mut self, reader: &mut fastq::Reader<io::BufReader<fs::File>>) -> Result<()> {
        let mut n = 0;
        for res in reader.records() {
            let record = res?;
            for b in record
                .sequence()
                .iter()
                .filter_map(|&b| self.filter_median_xmers(b))
            {
                //self.scp.;
            }
            n += 1;
        }
        eprintln!("Processed {} reads", n);
        Ok(())
    }
}
