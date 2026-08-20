//! Sequence utilities: reverse complement and related helpers.

/// IUPAC-aware complement of a single base, preserving case for ACGT/N and
/// passing through anything unexpected unchanged.
pub fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'N' => b'N',
        b'n' => b'n',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    }
}

/// Reverse complement of a nucleotide sequence.
pub fn reverse_complement(seq: &str) -> String {
    let bytes: Vec<u8> = seq.bytes().rev().map(complement_base).collect();
    String::from_utf8(bytes).expect("complement of ASCII is ASCII")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_basic() {
        assert_eq!(reverse_complement("ACGTN"), "NACGT");
        assert_eq!(reverse_complement("acgt"), "acgt");
    }
}
