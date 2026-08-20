//! Streaming reader/writer for TAMA BED12 files.

use std::io::{BufRead, Write};

use crate::error::Error;
use crate::model::BedTranscript;

/// Read all transcripts from a TAMA BED source. Blank lines and `track`/`#`
/// header lines are skipped.
pub fn read_bed<R: BufRead>(reader: R) -> Result<Vec<BedTranscript>, Error> {
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("track") {
            continue;
        }
        out.push(BedTranscript::parse_bed_line(&line)?);
    }
    Ok(out)
}

/// Write transcripts as TAMA BED12 lines (each terminated by `\n`).
pub fn write_bed<W: Write>(mut writer: W, transcripts: &[BedTranscript]) -> Result<(), Error> {
    for t in transcripts {
        writer.write_all(t.format_bed_line().as_bytes())?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
