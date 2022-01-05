use anyhow::{ensure, Result};
use bincode::deserialize_from;
use clap::ArgMatches;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub fn aln(matches: &ArgMatches) -> Result<()> {
    let fa_name = matches.value_of("ref").unwrap();

    let ks_name = format!("{}.ks", fa_name);
    ensure!(Path::new(&ks_name).exists(), "{} does not exist!", ks_name);
    let ks_file = BufReader::new(File::create(ks_name)?);
    let ks = deserialize_from(ks_file)?;

    Ok(())
}
