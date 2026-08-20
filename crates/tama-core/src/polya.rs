//! Poly-A run-on detection, ported from `detect_polya` in `tama_collapse.py`.
//!
//! Looks at the genomic window just 3' of a transcript (reverse-complemented on
//! the minus strand) and reports the fraction of A's — a signal that the read
//! ran into a genomic poly-A stretch rather than a real transcript end.

use crate::seq::reverse_complement;

/// Result of poly-A detection for one transcript.
#[derive(Debug, Clone, PartialEq)]
pub struct PolyA {
    pub downstream_seq: String,
    pub dseq_length: i64,
    pub a_count: i64,
    pub n_count: i64,
    pub a_percent: f64,
    pub n_percent: f64,
}

/// Detect a downstream poly-A window. `genome` is the full (uppercased) scaffold;
/// `start_pos`/`end_pos` are the transcript's 1-based bounds as produced by
/// [`crate::cigar::trans_coordinates`]. `a_window` defaults to 20 in the original.
pub fn detect_polya(
    genome: &[u8],
    strand: char,
    start_pos: i64,
    end_pos: i64,
    a_window: i64,
) -> PolyA {
    let downstream_seq: String = if strand == '+' {
        let s = end_pos.max(0) as usize;
        let e = ((end_pos + a_window).max(0) as usize).min(genome.len());
        let s = s.min(genome.len());
        String::from_utf8_lossy(&genome[s..e]).into_owned()
    } else {
        let trans_end = start_pos;
        let a_window_start = (trans_end - a_window).max(0);
        let s = (a_window_start.max(0) as usize).min(genome.len());
        let e = (trans_end.max(0) as usize).min(genome.len());
        let s = s.min(e);
        let raw = String::from_utf8_lossy(&genome[s..e]).into_owned();
        reverse_complement(&raw)
    };

    let mut dseq_length = downstream_seq.len() as i64;
    if dseq_length == 0 {
        dseq_length = 1;
    }
    let a_count = downstream_seq.bytes().filter(|&b| b == b'A').count() as i64;
    let n_count = downstream_seq.bytes().filter(|&b| b == b'N').count() as i64;
    let a_percent = a_count as f64 / dseq_length as f64;
    let n_percent = n_count as f64 / dseq_length as f64;

    PolyA {
        downstream_seq,
        dseq_length,
        a_count,
        n_count,
        a_percent,
        n_percent,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plus_strand_window() {
        // genome positions 0..: last aligned base at 1-based end_pos-1; window is
        // genome[end_pos..end_pos+window].
        let genome = b"CCCCCAAAAAAAAAA"; // A run starts at index 5
        let p = detect_polya(genome, '+', 1, 5, 10);
        assert_eq!(p.downstream_seq, "AAAAAAAAAA");
        assert_eq!(p.a_count, 10);
        assert!((p.a_percent - 1.0).abs() < 1e-9);
    }

    #[test]
    fn minus_strand_revcomp_window() {
        // On the minus strand we take the window upstream and reverse-complement.
        let genome = b"TTTTTGGGGG"; // upstream of pos 5 is TTTTT -> revcomp AAAAA
        let p = detect_polya(genome, '-', 5, 10, 5);
        assert_eq!(p.downstream_seq, "AAAAA");
        assert_eq!(p.a_count, 5);
    }
}
