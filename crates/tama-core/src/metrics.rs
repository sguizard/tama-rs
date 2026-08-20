//! Per-read coverage/identity metrics, ported from the `Transcript` methods
//! `calc_coverage`, `calc_identity`, and `make_error_line` in `tama_collapse.py`.

use crate::error_calc::ErrorRate;

/// Identity calculation method (`-icm`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdentMethod {
    /// `ident_cov`: includes hard/soft clipping in the denominator (default).
    IdentCov,
    /// `ident_map`: excludes hard/soft clipping.
    IdentMap,
}

/// Coverage/identity and the derived read report fields.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadMetrics {
    /// Full read length including hard-clipped bases (`seq_length + h_count`).
    pub length: i64,
    pub percent_coverage: f64,
    pub percent_identity: f64,
    /// `h;s;i;d;mis`.
    pub error_line: String,
}

/// Compute read metrics. `read_seq_len` is the length of the SAM SEQ field
/// (hard-clipped bases excluded); the hard-clip count is added back exactly as
/// the original does in `add_mismatch`.
pub fn read_metrics(read_seq_len: i64, e: &ErrorRate, method: IdentMethod) -> ReadMetrics {
    let seq_length = read_seq_len + e.h_count;
    let seq_length_f = seq_length as f64;

    let percent_coverage = (seq_length - e.h_count - e.s_count) as f64 / seq_length_f * 100.0;

    let percent_identity = match method {
        IdentMethod::IdentCov => {
            let nonmatch = e.h_count + e.s_count + e.i_count + e.d_count + e.mis_count;
            (seq_length - nonmatch) as f64 / seq_length_f * 100.0
        }
        IdentMethod::IdentMap => {
            let map_len = seq_length - e.h_count - e.s_count;
            let nonmatch = e.i_count + e.d_count + e.mis_count;
            (map_len - nonmatch) as f64 / map_len as f64 * 100.0
        }
    };

    ReadMetrics {
        length: seq_length,
        percent_coverage,
        percent_identity,
        error_line: format!(
            "{};{};{};{};{}",
            e.h_count, e.s_count, e.i_count, e.d_count, e.mis_count
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_read_is_full_coverage_identity() {
        let e = ErrorRate::default();
        let m = read_metrics(100, &e, IdentMethod::IdentCov);
        assert_eq!(m.length, 100);
        assert!((m.percent_coverage - 100.0).abs() < 1e-9);
        assert!((m.percent_identity - 100.0).abs() < 1e-9);
        assert_eq!(m.error_line, "0;0;0;0;0");
    }
}
