//! Grouping transcripts into genes by exon overlap.
//!
//! Ported from `gene_group` in `tama_collapse.py`. The original iterates over all
//! transcript pairs and merges their gene groups whenever any exon overlaps —
//! i.e. it computes the connected components of the "shares an overlapping exon"
//! graph. This implementation does the same via union-find, which yields an
//! identical partition. Strand is not considered here (callers split by strand
//! first). Genes are returned ordered by their minimum first-exon start, matching
//! the original's `start_gene_list.sort()`.

/// A transcript participating in gene grouping.
pub struct GeneMember<'a> {
    pub id: &'a str,
    pub exon_starts: &'a [i64],
    pub exon_ends: &'a [i64],
}

/// One gene group: the minimum first-exon start and its member transcript IDs
/// (in input order for determinism).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneGroup {
    pub gene_start: i64,
    pub trans_ids: Vec<String>,
}

/// Do any exons of `a` and `b` overlap? Half-open-inclusive test matching the
/// original `exon_start <= o_exon_end and exon_end >= o_exon_start`.
fn any_exon_overlap(a: &GeneMember, b: &GeneMember) -> bool {
    a.exon_starts.iter().zip(a.exon_ends).any(|(&as_, &ae)| {
        b.exon_starts
            .iter()
            .zip(b.exon_ends)
            .any(|(&bs, &be)| as_ <= be && ae >= bs)
    })
}

fn find(parent: &mut [usize], x: usize) -> usize {
    let mut root = x;
    while parent[root] != root {
        root = parent[root];
    }
    let mut cur = x;
    while parent[cur] != root {
        let next = parent[cur];
        parent[cur] = root;
        cur = next;
    }
    root
}

/// Group transcripts into genes by exon overlap.
#[allow(clippy::needless_range_loop)]
pub fn gene_group(members: &[GeneMember]) -> Vec<GeneGroup> {
    let n = members.len();
    let mut parent: Vec<usize> = (0..n).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            if any_exon_overlap(&members[i], &members[j]) {
                let ri = find(&mut parent, i);
                let rj = find(&mut parent, j);
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }

    // Collect members per root, preserving input order.
    let mut groups: Vec<(usize, GeneGroup)> = Vec::new();
    let mut root_to_idx: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        let first_start = members[i].exon_starts.first().copied().unwrap_or(i64::MAX);
        match root_to_idx.get(&r) {
            Some(&gi) => {
                let g = &mut groups[gi].1;
                g.trans_ids.push(members[i].id.to_string());
                if first_start < g.gene_start {
                    g.gene_start = first_start;
                }
            }
            None => {
                root_to_idx.insert(r, groups.len());
                groups.push((
                    r,
                    GeneGroup {
                        gene_start: first_start,
                        trans_ids: vec![members[i].id.to_string()],
                    },
                ));
            }
        }
    }

    let mut out: Vec<GeneGroup> = groups.into_iter().map(|(_, g)| g).collect();
    out.sort_by_key(|g| g.gene_start);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m<'a>(id: &'a str, s: &'a [i64], e: &'a [i64]) -> GeneMember<'a> {
        GeneMember { id, exon_starts: s, exon_ends: e }
    }

    #[test]
    fn separates_non_overlapping_and_orders_by_start() {
        let a_s = [100, 300];
        let a_e = [200, 400];
        let b_s = [1000];
        let b_e = [1100];
        let members = vec![m("b", &b_s, &b_e), m("a", &a_s, &a_e)];
        let groups = gene_group(&members);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].gene_start, 100); // sorted ascending
        assert_eq!(groups[0].trans_ids, vec!["a"]);
        assert_eq!(groups[1].trans_ids, vec!["b"]);
    }

    #[test]
    fn merges_transitively() {
        // a-b overlap, b-c overlap, a-c do not: all end in one gene.
        let a_s = [100];
        let a_e = [200];
        let b_s = [180];
        let b_e = [320];
        let c_s = [300];
        let c_e = [400];
        let members = vec![m("a", &a_s, &a_e), m("b", &b_s, &b_e), m("c", &c_s, &c_e)];
        let groups = gene_group(&members);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].trans_ids, vec!["a", "b", "c"]);
        assert_eq!(groups[0].gene_start, 100);
    }
}
