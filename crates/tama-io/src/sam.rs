//! Minimal SAM record access for the collapse driver.
//!
//! TAMA reads a sorted SAM as plain TSV, taking the query name, flag, reference,
//! 1-based position, CIGAR, SEQ, and the optional `XS:A:` strand tag. This mirror
//! keeps that behaviour. (BAM support via noodles is a planned addition.)

use std::io::BufRead;
use std::path::Path;

use crate::open_reader;

/// The fields of a SAM alignment record used by collapse.
#[derive(Debug, Clone)]
pub struct SamRecord {
    pub read_id: String,
    pub flag: u32,
    pub scaffold: String,
    pub start_pos: i64,
    pub cigar: String,
    pub seq: String,
    /// Strand from the `XS:A:` tag (`+`/`-`), or `None` if absent.
    pub xs_strand: Option<char>,
}

impl SamRecord {
    /// Parse a single non-header SAM line. Returns `None` for header/blank lines.
    pub fn parse(line: &str) -> Option<SamRecord> {
        if line.is_empty() || line.starts_with('@') {
            return None;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            return None;
        }
        let xs_strand = f.iter().find_map(|field| {
            field
                .strip_prefix("XS:A:")
                .and_then(|s| s.chars().next())
        });
        Some(SamRecord {
            read_id: f[0].to_string(),
            flag: f[1].parse().ok()?,
            scaffold: f[2].to_string(),
            start_pos: f[3].parse().ok()?,
            cigar: f[5].to_string(),
            seq: f[9].to_string(),
            xs_strand,
        })
    }
}

/// Read all alignment records from a SAM file in file order.
pub fn read_sam<P: AsRef<Path>>(path: P) -> std::io::Result<Vec<SamRecord>> {
    let reader = open_reader(path)?;
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rec) = SamRecord::parse(&line) {
            out.push(rec);
        }
    }
    Ok(out)
}

/// TAMA's SAM-flag categories.
pub fn mapped_flag(flag: u32) -> &'static str {
    match flag {
        0 => "forward_strand",
        4 => "unmapped",
        16 => "reverse_strand",
        256 | 272 => "not_primary",
        2048 | 2064 => "chimeric",
        _ => "unknown",
    }
}
