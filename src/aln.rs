use anyhow::Result;
use bio::io::fasta::IndexedReader;
use clap::ArgMatches;

pub fn aln(matches: &ArgMatches) -> Result<()> {
    let fa_name = matches.value_of("ref").unwrap();

    let mut idxr = IndexedReader::from_file(&fa_name)
        .unwrap_or_else(|_| panic!("Error opening reference genome"));
    let chrs = idxr.index.sequences();
    Ok(())
}
