//! Minimal FASTA loading.
//!
//! TAMA's genome inputs are small scaffolds, so an in-memory map keyed by the
//! first whitespace-delimited token of each header is sufficient and matches how
//! the Python code indexes sequences by scaffold name.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use crate::open_reader;

/// Load a FASTA file into `name -> sequence` (uppercased), keyed by the header
/// ID (text up to the first whitespace).
pub fn load_fasta<P: AsRef<Path>>(path: P) -> std::io::Result<HashMap<String, String>> {
    let reader = open_reader(path)?;
    let mut map: HashMap<String, String> = HashMap::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = cur_name.take() {
                map.insert(name, std::mem::take(&mut cur_seq));
            }
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            cur_name = Some(id);
        } else {
            cur_seq.push_str(line.trim_end());
        }
    }
    if let Some(name) = cur_name.take() {
        map.insert(name, cur_seq);
    }
    for seq in map.values_mut() {
        seq.make_ascii_uppercase();
    }
    Ok(map)
}
