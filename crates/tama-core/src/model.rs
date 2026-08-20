//! Core transcript model and the TAMA BED12 dialect.
//!
//! TAMA uses a standard 12-column BED where column 4 (`name`) packs the gene and
//! transcript IDs as `gene_id;trans_id`, and columns 7/8 (`thickStart`/`thickEnd`)
//! carry the CDS region (equal to the transcript bounds when there is no CDS).
//!
//! Coordinate convention follows `tama_merge.py::Transcript.add_bed_info`:
//! `trans_start`/`trans_end` are the 0-based BED bounds, and each absolute exon
//! coordinate is `trans_start + relative_block_offset`. Writing inverts this
//! exactly, so parse→format round-trips byte-for-byte for any TAMA BED line.

use std::fmt::Write as _;

use crate::error::Error;

/// A transcript model as represented in a TAMA BED12 line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BedTranscript {
    pub scaffold: String,
    /// 0-based BED start (`trans_start`).
    pub trans_start: i64,
    /// BED end (`trans_end`).
    pub trans_end: i64,
    /// Full `name` field, e.g. `G1;G1.1`.
    pub id_line: String,
    pub score: String,
    pub strand: char,
    /// CDS start (BED `thickStart`).
    pub cds_start: i64,
    /// CDS end (BED `thickEnd`).
    pub cds_end: i64,
    pub rgb: String,
    /// Absolute exon start coordinates (`trans_start + relative_offset`).
    pub exon_start_list: Vec<i64>,
    /// Absolute exon end coordinates.
    pub exon_end_list: Vec<i64>,
}

impl BedTranscript {
    /// Gene ID: the part of `id_line` before the first `;`.
    pub fn gene_id(&self) -> &str {
        self.id_line.split(';').next().unwrap_or(&self.id_line)
    }

    /// Transcript ID: the part of `id_line` after the first `;`, or the whole
    /// field if there is no `;`.
    pub fn trans_id(&self) -> &str {
        match self.id_line.split_once(';') {
            Some((_, t)) => t,
            None => &self.id_line,
        }
    }

    pub fn num_exons(&self) -> usize {
        self.exon_start_list.len()
    }

    /// Parse one TAMA BED12 line. Trailing commas in the block fields are
    /// tolerated, matching the Python parser.
    pub fn parse_bed_line(line: &str) -> Result<Self, Error> {
        let cols: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
        if cols.len() < 12 {
            return Err(Error::Bed(format!(
                "expected >=12 tab-separated columns, found {}: {line:?}",
                cols.len()
            )));
        }
        let parse_i = |s: &str, what: &str| -> Result<i64, Error> {
            s.trim()
                .parse::<i64>()
                .map_err(|_| Error::Bed(format!("invalid {what}: {s:?}")))
        };

        let trans_start = parse_i(cols[1], "trans_start")?;
        let trans_end = parse_i(cols[2], "trans_end")?;
        let strand = cols[5].chars().next().unwrap_or('+');
        let num_exons = parse_i(cols[9], "block count")? as usize;

        fn split_blocks(s: &str) -> Vec<&str> {
            s.split(',').filter(|f| !f.is_empty()).collect()
        }
        let block_sizes = split_blocks(cols[10]);
        let block_starts = split_blocks(cols[11]);
        if block_sizes.len() != block_starts.len() {
            return Err(Error::Bed(format!(
                "block size/start count mismatch ({} vs {})",
                block_sizes.len(),
                block_starts.len()
            )));
        }

        let mut exon_start_list = Vec::with_capacity(block_sizes.len());
        let mut exon_end_list = Vec::with_capacity(block_sizes.len());
        for (bs, bstart) in block_sizes.iter().zip(&block_starts) {
            let rel_start = parse_i(bstart, "block start")?;
            let size = parse_i(bs, "block size")?;
            exon_start_list.push(trans_start + rel_start);
            exon_end_list.push(trans_start + rel_start + size);
        }
        if num_exons != exon_start_list.len() {
            return Err(Error::Bed(format!(
                "block count {num_exons} does not match number of blocks {}",
                exon_start_list.len()
            )));
        }

        Ok(BedTranscript {
            scaffold: cols[0].to_string(),
            trans_start,
            trans_end,
            id_line: cols[3].to_string(),
            score: cols[4].to_string(),
            strand,
            cds_start: parse_i(cols[6], "cds_start")?,
            cds_end: parse_i(cols[7], "cds_end")?,
            rgb: cols[8].to_string(),
            exon_start_list,
            exon_end_list,
        })
    }

    /// Serialize back to a TAMA BED12 line (no trailing newline).
    pub fn format_bed_line(&self) -> String {
        let mut sizes = String::new();
        let mut starts = String::new();
        for (i, (&s, &e)) in self
            .exon_start_list
            .iter()
            .zip(&self.exon_end_list)
            .enumerate()
        {
            if i > 0 {
                sizes.push(',');
                starts.push(',');
            }
            let _ = write!(sizes, "{}", e - s);
            let _ = write!(starts, "{}", s - self.trans_start);
        }

        [
            self.scaffold.clone(),
            self.trans_start.to_string(),
            self.trans_end.to_string(),
            self.id_line.clone(),
            self.score.clone(),
            self.strand.to_string(),
            self.cds_start.to_string(),
            self.cds_end.to_string(),
            self.rgb.clone(),
            self.num_exons().to_string(),
            sizes,
            starts,
        ]
        .join("\t")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LINE: &str = "AADN04001032.1\t866\t2296\tG1;G1.1\t40\t+\t866\t2296\t255,0,0\t3\t100,50,80\t0,500,1350";

    #[test]
    fn round_trips() {
        let t = BedTranscript::parse_bed_line(LINE).unwrap();
        assert_eq!(t.num_exons(), 3);
        assert_eq!(t.gene_id(), "G1");
        assert_eq!(t.trans_id(), "G1.1");
        assert_eq!(t.exon_start_list, vec![866, 1366, 2216]);
        assert_eq!(t.exon_end_list, vec![966, 1416, 2296]);
        assert_eq!(t.format_bed_line(), LINE);
    }

    #[test]
    fn tolerates_trailing_commas() {
        let line = LINE.to_string() + ",";
        let with_trailing = line.replace("0,500,1350", "0,500,1350,");
        let t = BedTranscript::parse_bed_line(&with_trailing).unwrap();
        assert_eq!(t.num_exons(), 3);
    }
}
