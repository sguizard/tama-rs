//! No-cap transcript grouping, ported from `simplify_gene_nocap` and the nocap
//! parts of the `TransGroup` class in `tama_collapse.py`.
//!
//! In no-cap mode `compare_transcripts` never returns `same_transcript` (that is
//! capped-only), so `merge_a_b_groups_nocap` is never reached. The grouping
//! therefore reduces to `add_a_to_b_group`:
//!
//! * same exon level, `same_three_prime_same_exons`: a 5'-shorter read is added
//!   to the hunter's group (isoforms sharing a 3' backbone collapse together);
//! * lower exon level, `same_three_prime_diff_exons`: a 5'-degraded shorter read
//!   is absorbed into the longer read's group(s).
//!
//! A read may belong to multiple groups (a degraded read can support several
//! full-length models), so groups are returned as (possibly overlapping) index
//! lists.

use std::collections::HashMap;

use indexmap::{IndexMap, IndexSet};

use crate::collapse::{compare_transcripts, CollapseParams, CompFlag, Transcript};

/// Group-membership bookkeeping (nocap subset of the original `TransGroup`).
#[derive(Default)]
struct TransGroup {
    trans_group: IndexMap<String, IndexSet<i64>>,
    group_trans: IndexMap<i64, IndexSet<String>>,
    group_count: i64,
}

impl TransGroup {
    fn has(&self, t: &str) -> bool {
        self.trans_group.contains_key(t)
    }

    fn same_group(&self, a: &str, b: &str) -> bool {
        match (self.trans_group.get(a), self.trans_group.get(b)) {
            (Some(ga), Some(gb)) => ga.iter().any(|g| gb.contains(g)),
            _ => false,
        }
    }

    fn new_group_a(&mut self, t: &str) {
        self.group_count += 1;
        let g = self.group_count;
        self.trans_group.entry(t.to_string()).or_default().insert(g);
        let mut set = IndexSet::new();
        set.insert(t.to_string());
        self.group_trans.insert(g, set);
    }

    /// Add short read `a` to every group of long read `b`.
    fn add_a_to_b_group(&mut self, a: &str, b: &str) {
        // Drop `a`'s initial lone self-identity group before absorbing it.
        if let Some(groups) = self.trans_group.get(a) {
            if groups.len() == 1 {
                let g = *groups.iter().next().unwrap();
                if self
                    .group_trans
                    .get(&g)
                    .map(|s| s.len() == 1 && s.contains(a))
                    .unwrap_or(false)
                {
                    self.group_trans.shift_remove(&g);
                    self.trans_group.shift_remove(a);
                }
            }
        }
        self.trans_group.entry(a.to_string()).or_default();

        let b_groups: Vec<i64> = self
            .trans_group
            .get(b)
            .map(|s| s.iter().copied().collect())
            .unwrap_or_default();
        for g in b_groups {
            self.trans_group.get_mut(a).unwrap().insert(g);
            self.group_trans.entry(g).or_default().insert(a.to_string());
        }
    }
}

/// Group a gene's transcripts for no-cap collapsing. Returns member index lists
/// (a transcript index may appear in more than one group).
pub fn simplify_gene_nocap(trans: &[Transcript], p: &CollapseParams) -> Vec<Vec<usize>> {
    let idx: HashMap<&str, usize> = trans
        .iter()
        .enumerate()
        .map(|(i, t)| (t.cluster_id.as_str(), i))
        .collect();
    let strand = trans[0].strand;

    // Transcripts grouped by exon count, in input order.
    let mut exon_trans: IndexMap<usize, Vec<String>> = IndexMap::new();
    for t in trans {
        exon_trans
            .entry(t.num_exons())
            .or_default()
            .push(t.cluster_id.clone());
    }
    let mut level_list: Vec<usize> = exon_trans.keys().copied().collect();
    level_list.sort_unstable_by(|a, b| b.cmp(a)); // descending

    let mut tg = TransGroup::default();
    let mut sub_length: IndexSet<String> = IndexSet::new();

    let five_prime = |id: &str| -> i64 {
        let t = &trans[idx[id]];
        if strand == '+' {
            t.start_pos
        } else {
            t.end_pos
        }
    };

    for exon_num_index in 0..level_list.len() {
        let level = level_list[exon_num_index];
        // ungrouped clusters at this level (self-mutating across outer iterations)
        let mut ungrouped: IndexSet<String> = exon_trans[&level].iter().cloned().collect();

        while !ungrouped.is_empty() {
            // hunter: smallest 5' coord (largest for minus strand), then input order
            let hunter = {
                let mut ids: Vec<String> = ungrouped.iter().cloned().collect();
                ids.sort_by(|a, b| {
                    let (fa, fb) = (five_prime(a), five_prime(b));
                    if strand == '+' {
                        fa.cmp(&fb)
                    } else {
                        fb.cmp(&fa)
                    }
                });
                ids[0].clone()
            };
            ungrouped.shift_remove(&hunter);

            let mut unsearched: IndexSet<String> = IndexSet::new();
            let mut next_hunter = Some(hunter);

            loop {
                let h = match next_hunter.take() {
                    Some(h) => h,
                    None => {
                        if unsearched.is_empty() {
                            break;
                        }
                        let mut ks: Vec<String> = unsearched.iter().cloned().collect();
                        ks.sort();
                        let h = ks[0].clone();
                        unsearched.shift_remove(&h);
                        h
                    }
                };
                if !tg.has(&h) {
                    tg.new_group_a(&h);
                }

                // same exon level search (skip if this read was absorbed upward)
                if !sub_length.contains(&h) {
                    let preys: Vec<String> = ungrouped.iter().cloned().collect();
                    for prey in preys {
                        if prey == h {
                            continue;
                        }
                        if !tg.has(&prey) {
                            tg.new_group_a(&prey);
                        }
                        if tg.same_group(&h, &prey) {
                            continue;
                        }
                        let flag = compare_transcripts(
                            &trans[idx[h.as_str()]],
                            &trans[idx[prey.as_str()]],
                            p,
                            strand,
                        );
                        if flag == CompFlag::SameThreePrimeSameExons {
                            tg.add_a_to_b_group(&prey, &h);
                            unsearched.insert(prey);
                        }
                    }
                }

                // lower exon level search: absorb shorter, 3'-compatible reads
                for &prey_level in &level_list[exon_num_index + 1..] {
                    let preys: Vec<String> = exon_trans[&prey_level].clone();
                    for prey in preys {
                        if prey == h {
                            continue;
                        }
                        if !tg.has(&prey) {
                            tg.new_group_a(&prey);
                        }
                        if tg.same_group(&h, &prey) {
                            continue;
                        }
                        let flag = compare_transcripts(
                            &trans[idx[h.as_str()]],
                            &trans[idx[prey.as_str()]],
                            p,
                            strand,
                        );
                        if flag == CompFlag::SameThreePrimeDiffExons {
                            tg.add_a_to_b_group(&prey, &h);
                            sub_length.insert(prey);
                        }
                    }
                }

                if unsearched.is_empty() {
                    break;
                }
            }
        }
    }

    tg.group_trans
        .values()
        .map(|members| members.iter().map(|c| idx[c.as_str()]).collect())
        .collect()
}
