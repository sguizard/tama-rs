//! Bioinformatics format I/O for TAMA.
//!
//! Phase 1 provides transparent (optionally gzip-compressed) file readers/writers
//! and a simple FASTA loader. SAM/BAM parsing via `noodles` is added in Phase 2
//! alongside `tama collapse`.

pub mod fasta;
pub mod sam;

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;

/// Open a file as a buffered reader, transparently decompressing `.gz` inputs.
pub fn open_reader<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn BufRead>> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
    if is_gz {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Create a buffered writer for `path`.
pub fn create_writer<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn Write>> {
    Ok(Box::new(BufWriter::new(File::create(path)?)))
}
