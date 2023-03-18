#[cfg(test)]
mod tests {
    use super::*;
    use fa2b2::index::parse_fasta_file;
    use noodles_fasta::{
        self as fasta,
        record::{Definition, Sequence},
    };
    use std::io::Write;
    use std::path::PathBuf;

    #[test]
    fn test_parse_fasta() {
        // Create a temporary file with sample FASTA data
        let mut temp_file = tempfile::NamedTempFile::new().unwrap();
        writeln!(temp_file, ">seq1\nATCG\n>seq2\nGCTA").unwrap();
        let file_path = PathBuf::from(temp_file.path());

        // Parse the FASTA file

        let mut reader = parse_fasta_file(file_path).unwrap();
        let mut records = reader.records();

        // Check that the sequence IDs and sequences match the input file
        assert_eq!(
            records.next().transpose().unwrap(),
            Some(fasta::Record::new(
                Definition::new("seq1", None),
                Sequence::from(b"ATCG".to_vec()),
            ))
        );
        assert_eq!(
            records.next().transpose().unwrap(),
            Some(fasta::Record::new(
                Definition::new("seq2", None),
                Sequence::from(b"GCTA".to_vec()),
            ))
        );
        assert_eq!(records.next().transpose().unwrap(), None);
    }
}
