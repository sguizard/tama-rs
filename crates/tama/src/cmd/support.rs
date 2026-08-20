//! `tama support` — read-support tracking tools.

use std::io::{BufRead, Write};

use anyhow::bail;
use clap::{Args as ClapArgs, Subcommand};
use indexmap::{IndexMap, IndexSet};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Read support from a collapse trans_read.bed + cluster file. (tama_read_support_collapse_cluster)
    CollapseCluster {
        /// Collapse `<prefix>_trans_read.bed`.
        collapse: std::path::PathBuf,
        /// Cluster file (Iso-Seq1 CSV, Iso-Seq3, or a BED for no-cluster mode).
        cluster: std::path::PathBuf,
        /// Output file.
        output: std::path::PathBuf,
    },
    /// Read support levels. (tama_read_support_levels)
    Levels,
    /// Read support after merging collapse outputs. (tama_read_support_merge_collapse)
    MergeCollapse {
        /// tama merge `<prefix>_merge.txt`.
        merge: std::path::PathBuf,
        /// Filelist: `support_filename<TAB>prefix<TAB>dir/` per source.
        filelist: std::path::PathBuf,
        /// Output file.
        output: std::path::PathBuf,
    },
}

pub fn run(args: Args) -> anyhow::Result<()> {
    match args.cmd {
        Cmd::CollapseCluster { collapse, cluster, output } => {
            collapse_cluster(&collapse, &cluster, &output)
        }
        Cmd::Levels => Err(super::not_implemented("support levels")),
        Cmd::MergeCollapse { merge, filelist, output } => {
            merge_collapse(&merge, &filelist, &output)
        }
    }
}

/// Read support per merged transcript from collapse support files.
/// Ports `tama_read_support_merge_collapse`.
fn merge_collapse(
    merge: &std::path::Path,
    filelist: &std::path::Path,
    output: &std::path::Path,
) -> anyhow::Result<()> {
    // prefix -> trans_id -> trans_num_reads
    let mut trans_read: IndexMap<String, IndexMap<String, i64>> = IndexMap::new();
    // prefix -> gene -> trans -> set(read)
    type ReadTree = IndexMap<String, IndexMap<String, IndexMap<String, IndexSet<String>>>>;
    let mut gene_trans_read: ReadTree = IndexMap::new();
    let mut prefix_list: Vec<String> = Vec::new();

    for fline in read_lines(filelist)? {
        if fline.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = fline.split('\t').collect();
        let (filename, prefix, fpath) = (f[0], f[1].to_string(), f[2]);
        prefix_list.push(prefix.clone());
        let support_path = std::path::PathBuf::from(format!("{fpath}{filename}"));
        if trans_read.contains_key(&prefix) {
            bail!("duplicate prefix {prefix}");
        }
        trans_read.insert(prefix.clone(), IndexMap::new());
        gene_trans_read.insert(prefix.clone(), IndexMap::new());

        for line in read_lines(&support_path)? {
            if line.starts_with("gene_id") {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            let (gene_id, trans_id, trans_num_reads, cluster_line) =
                (c[0].to_string(), c[1].to_string(), c[3].parse::<i64>()?, c[4]);
            let mut reads: Vec<String> = Vec::new();
            for cg in cluster_line.split(';') {
                if let Some(rl) = cg.split(':').nth(1) {
                    reads.extend(rl.split(',').map(String::from));
                }
            }
            trans_read.get_mut(&prefix).unwrap().insert(trans_id.clone(), trans_num_reads);
            let g = gene_trans_read
                .get_mut(&prefix)
                .unwrap()
                .entry(gene_id)
                .or_default()
                .entry(trans_id)
                .or_default();
            for r in reads {
                g.insert(r);
            }
        }
    }

    // parse merge file
    let mut merge_gene_list: Vec<String> = Vec::new();
    // merge_gene -> merge_trans -> prefix -> collapse_id -> set(read)
    type MergeTree = IndexMap<String, IndexMap<String, IndexMap<String, IndexMap<String, IndexSet<String>>>>>;
    let mut merge_tree: MergeTree = IndexMap::new();
    let mut merge_gene_trans: IndexMap<String, Vec<String>> = IndexMap::new();

    for line in read_lines(merge)? {
        let cols: Vec<&str> = line.split('\t').collect();
        let id_split: Vec<&str> = cols[3].split(';').collect();
        let trans_id = id_split[0].to_string();
        let support_id = id_split[1];
        let gene_id = trans_id.split('.').next().unwrap_or(&trans_id).to_string();
        let prefix = support_id.split('_').next().unwrap_or("").to_string();
        let collapse_id = support_id.split_once('_').map(|x| x.1).unwrap_or("").to_string();

        if !merge_gene_trans.contains_key(&gene_id) {
            merge_gene_list.push(gene_id.clone());
            merge_gene_trans.insert(gene_id.clone(), Vec::new());
        }
        if !merge_gene_trans[&gene_id].contains(&trans_id) {
            merge_gene_trans.get_mut(&gene_id).unwrap().push(trans_id.clone());
        }

        let source_gene_id = collapse_id.split('.').next().unwrap_or(&collapse_id).to_string();
        let reads: Vec<String> = gene_trans_read
            .get(&prefix)
            .and_then(|g| g.get(&source_gene_id))
            .and_then(|t| t.get(&collapse_id))
            .map(|s| s.iter().cloned().collect())
            .unwrap_or_default();

        let slot = merge_tree
            .entry(gene_id)
            .or_default()
            .entry(trans_id)
            .or_default()
            .entry(prefix)
            .or_default()
            .entry(collapse_id)
            .or_default();
        for r in reads {
            slot.insert(r);
        }
    }

    let mut out = tama_io::create_writer(output)?;
    writeln!(out, "merge_gene_id\tmerge_trans_id\tgene_read_support\ttrans_read_support\tsource_prefix\tsource_trans_line\tsource_read_line")?;

    for gene in &merge_gene_list {
        let mut gene_reads: IndexSet<String> = IndexSet::new();
        for tmap in merge_tree[gene].values() {
            for pmap in tmap.values() {
                for reads in pmap.values() {
                    gene_reads.extend(reads.iter().cloned());
                }
            }
        }
        let gene_read_support = gene_reads.len();

        for trans_id in &merge_gene_trans[gene] {
            let mut trans_reads: IndexSet<String> = IndexSet::new();
            let mut prefixes: Vec<String> = Vec::new();
            let mut support_trans: Vec<String> = Vec::new();
            let mut read_lines_out: Vec<String> = Vec::new();
            for (prefix, cmap) in &merge_tree[gene][trans_id] {
                prefixes.push(prefix.clone());
                for (collapse_id, reads) in cmap {
                    support_trans.push(format!("{prefix}_{collapse_id}"));
                    let rl: Vec<String> = reads.iter().cloned().collect();
                    read_lines_out.push(rl.join(","));
                    trans_reads.extend(reads.iter().cloned());
                }
            }
            writeln!(
                out,
                "{gene}\t{trans_id}\t{gene_read_support}\t{}\t{}\t{}\t{}",
                trans_reads.len(),
                prefixes.join(","),
                support_trans.join(","),
                read_lines_out.join(";")
            )?;
        }
    }
    Ok(())
}

#[derive(PartialEq)]
enum ClusterType {
    V1,
    NoCluster,
}

/// Read support per collapsed transcript from a clustering file.
/// Ports `tama_read_support_collapse_cluster`.
fn collapse_cluster(
    collapse: &std::path::Path,
    cluster: &std::path::Path,
    output: &std::path::Path,
) -> anyhow::Result<()> {
    let cluster_lines: Vec<String> = read_lines(cluster)?;
    if cluster_lines.is_empty() {
        bail!("empty cluster file");
    }
    // Detect cluster file type from the first line.
    let first: Vec<&str> = cluster_lines[0].split('\t').collect();
    let ctype = if first.len() != 12 {
        ClusterType::V1
    } else if first[3].starts_with('G') && first[3].split('.').count() == 2 {
        ClusterType::NoCluster
    } else {
        bail!("cannot understand cluster file type");
    };

    // cluster_id -> set(read_id)
    let mut cluster_read: IndexMap<String, IndexSet<String>> = IndexMap::new();
    for line in &cluster_lines {
        if line.starts_with("cluster_id") || line.starts_with("from") {
            continue;
        }
        let (cluster_id, read_id) = match ctype {
            ClusterType::V1 => {
                let f: Vec<&str> = line.split(',').collect();
                (f[0].to_string(), f[1].to_string())
            }
            ClusterType::NoCluster => {
                let f: Vec<&str> = line.split('\t').collect();
                let read = f[3].split(';').nth(1).unwrap_or("").to_string();
                (read.clone(), read)
            }
        };
        cluster_read.entry(cluster_id).or_default().insert(read_id);
    }

    // Walk the collapse trans_read.bed, mapping transcripts -> clusters.
    let mut gene_order: Vec<String> = Vec::new();
    let mut gene_trans: IndexMap<String, Vec<String>> = IndexMap::new();
    let mut trans_cluster: IndexMap<String, IndexSet<String>> = IndexMap::new();

    for line in read_lines(collapse)? {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            continue;
        }
        let id_split: Vec<&str> = cols[3].split(';').collect();
        let trans_id = id_split[0].to_string();
        let cluster_id = match ctype {
            ClusterType::V1 => id_split[1].split('/').next().unwrap_or("").to_string(),
            ClusterType::NoCluster => id_split[1].to_string(),
        };
        let gene_id = trans_id.split('.').next().unwrap_or(&trans_id).to_string();

        if !gene_trans.contains_key(&gene_id) {
            gene_order.push(gene_id.clone());
        }
        if !trans_cluster.contains_key(&trans_id) {
            gene_trans.entry(gene_id).or_default().push(trans_id.clone());
        }
        trans_cluster.entry(trans_id).or_default().insert(cluster_id);
    }

    let mut out = tama_io::create_writer(output)?;
    writeln!(out, "gene_id\ttrans_id\tgene_num_reads\ttrans_num_reads\tcluster_line")?;

    for gene_id in &gene_order {
        let gene_trans_list = &gene_trans[gene_id];
        // gene read count
        let mut gene_num_reads = 0usize;
        for trans_id in gene_trans_list {
            for cluster_id in &trans_cluster[trans_id] {
                gene_num_reads += cluster_read.get(cluster_id).map(|s| s.len()).unwrap_or(0);
            }
        }
        // per transcript
        for trans_id in gene_trans_list {
            let mut trans_num_reads = 0usize;
            let mut cluster_parts: Vec<String> = Vec::new();
            for cluster_id in &trans_cluster[trans_id] {
                let reads: Vec<String> = cluster_read
                    .get(cluster_id)
                    .map(|s| s.iter().cloned().collect())
                    .unwrap_or_default();
                trans_num_reads += reads.len();
                cluster_parts.push(format!("{}:{}", cluster_id, reads.join(",")));
            }
            writeln!(
                out,
                "{gene_id}\t{trans_id}\t{gene_num_reads}\t{trans_num_reads}\t{}",
                cluster_parts.join(";")
            )?;
        }
    }
    Ok(())
}

fn read_lines(path: &std::path::Path) -> anyhow::Result<Vec<String>> {
    let reader = tama_io::open_reader(path)?;
    let mut out = Vec::new();
    for line in reader.lines() {
        out.push(line?);
    }
    // mimic rstrip("\n").split("\n"): drop a single trailing empty line
    while out.last().map(|l| l.is_empty()).unwrap_or(false) {
        out.pop();
    }
    Ok(out)
}
