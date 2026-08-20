//! CIGAR parsing and coordinate derivation, ported from `tama_collapse.py`.
//!
//! Mirrors `cigar_list`, `mapped_seq_length`, and `trans_coordinates`. As in the
//! original, `=` and `X` ops are normalised to `M` (minimap2 emits these).

use crate::error::Error;

/// One CIGAR operation: a length and a single-character op code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    pub len: i64,
    pub op: u8,
}

/// Parse a CIGAR string into `(len, op)` pairs. `=`/`X` become `M`.
pub fn cigar_list(cigar: &str) -> Result<Vec<CigarOp>, Error> {
    let mut ops = Vec::new();
    let mut num = String::new();
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            let len = num
                .parse::<i64>()
                .map_err(|_| Error::Cigar(format!("bad length before {c:?} in {cigar:?}")))?;
            num.clear();
            let op = match c {
                '=' | 'X' => b'M',
                other => other as u8,
            };
            ops.push(CigarOp { len, op });
        }
    }
    if !num.is_empty() {
        return Err(Error::Cigar(format!("trailing digits in {cigar:?}")));
    }
    Ok(ops)
}

/// Length of the mapped reference span: sum of `M` and `D` op lengths.
pub fn mapped_seq_length(cigar: &str) -> Result<i64, Error> {
    let mut len = 0;
    for c in cigar_list(cigar)? {
        match c.op {
            b'M' | b'D' => len += c.len,
            _ => {}
        }
    }
    Ok(len)
}

/// Result of walking a CIGAR from a 1-based genomic `start_pos`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TransCoordinates {
    /// 1-based coordinate just past the last aligned base.
    pub end_pos: i64,
    /// Absolute 1-based exon start coordinates.
    pub exon_start_list: Vec<i64>,
    /// Absolute exon end coordinates.
    pub exon_end_list: Vec<i64>,
}

/// Derive transcript/exon coordinates from a CIGAR. `N` ops split exons.
pub fn trans_coordinates(start_pos: i64, cigar: &str) -> Result<TransCoordinates, Error> {
    let mut end_pos = start_pos;
    let mut exon_start_list = vec![start_pos];
    let mut exon_end_list = Vec::new();

    for c in cigar_list(cigar)? {
        match c.op {
            b'M' | b'D' => end_pos += c.len,
            b'N' => {
                exon_end_list.push(end_pos);
                end_pos += c.len;
                exon_start_list.push(end_pos);
            }
            _ => {} // H, S, I consume no reference span
        }
    }
    exon_end_list.push(end_pos);

    Ok(TransCoordinates {
        end_pos,
        exon_start_list,
        exon_end_list,
    })
}

/// IUPAC ambiguity codes: does reference `code` match query base `base`?
/// Ported from `nuc_char_dict` in `tama_collapse.py`. Both are compared as
/// uppercase ASCII.
pub fn iupac_match(code: u8, base: u8) -> bool {
    let code = code.to_ascii_uppercase();
    let base = base.to_ascii_uppercase();
    let set: &[u8] = match code {
        b'A' => b"A",
        b'T' => b"T",
        b'C' => b"C",
        b'G' => b"G",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"AC",
        b'M' => b"GT",
        b'Y' => b"AG",
        b'R' => b"CT",
        b'V' => b"CGT",
        b'H' => b"AGT",
        b'D' => b"ACT",
        b'B' => b"ACG",
        b'N' => b"ATCG",
        _ => return code == base,
    };
    set.contains(&base)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_normalises() {
        let ops = cigar_list("10M5=3X2N4M").unwrap();
        assert_eq!(ops[1], CigarOp { len: 5, op: b'M' });
        assert_eq!(ops[2], CigarOp { len: 3, op: b'M' });
    }

    #[test]
    fn mapped_length_counts_m_and_d() {
        // 10M + 1D + 8M = 19; N and S do not count.
        assert_eq!(mapped_seq_length("5S10M1D8M3N4M").unwrap(), 23);
    }

    #[test]
    fn splits_exons_on_n() {
        let tc = trans_coordinates(100, "50M20N30M").unwrap();
        assert_eq!(tc.exon_start_list, vec![100, 170]);
        assert_eq!(tc.exon_end_list, vec![150, 200]);
        assert_eq!(tc.end_pos, 200);
    }

    #[test]
    fn iupac_codes() {
        assert!(iupac_match(b'N', b'a'));
        assert!(iupac_match(b'S', b'G'));
        assert!(!iupac_match(b'S', b'A'));
        assert!(iupac_match(b'A', b'A'));
    }
}
