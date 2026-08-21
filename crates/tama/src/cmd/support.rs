//! `tama support` — read-support tracking tools.

use std::io::{BufRead, Write};

use anyhow::bail;
use clap::{Args as ClapArgs, Subcommand, ValueEnum};
use indexmap::{IndexMap, IndexSet};

/// Layout of the merge file passed to `support levels`. (`-mt`)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "snake_case")]
pub enum MergeType {
    /// A `tama merge` `_merge.txt` / `trans_read.bed` file.
    Tama,
    /// A cupcake collapse group file.
    Cupcake,
    /// A `tama filter` report.
    Filter,
}

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
    /// Read support levels file. (tama_read_support_levels)
    Levels {
        /// Filelist: `source_name<TAB>transread_file<TAB>file_type`. (`-f`)
        #[arg(short = 'f', long)]
        filelist: std::path::PathBuf,
        /// Merge `_merge.txt`, or `no_merge`. (`-m`)
        #[arg(short = 'm', long)]
        merge: String,
        /// Output prefix (writes `<prefix>_read_support.txt`). (`-o`)
        #[arg(short = 'o', long)]
        output: String,
        /// Merge file layout. (`-mt`)
        #[arg(long = "mt", value_enum, default_value = "tama")]
        merge_type: MergeType,
    },
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
        Cmd::CollapseCluster {
            collapse,
            cluster,
            output,
        } => collapse_cluster(&collapse, &cluster, &output),
        Cmd::Levels {
            filelist,
            merge,
            output,
            merge_type,
        } => levels(&filelist, &merge, &output, merge_type),
        Cmd::MergeCollapse {
            merge,
            filelist,
            output,
        } => merge_collapse(&merge, &filelist, &output),
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
            let (gene_id, trans_id, trans_num_reads, cluster_line) = (
                c[0].to_string(),
                c[1].to_string(),
                c[3].parse::<i64>()?,
                c[4],
            );
            let mut reads: Vec<String> = Vec::new();
            for cg in cluster_line.split(';') {
                if let Some(rl) = cg.split(':').nth(1) {
                    reads.extend(rl.split(',').map(String::from));
                }
            }
            trans_read
                .get_mut(&prefix)
                .unwrap()
                .insert(trans_id.clone(), trans_num_reads);
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
    type MergeTree =
        IndexMap<String, IndexMap<String, IndexMap<String, IndexMap<String, IndexSet<String>>>>>;
    let mut merge_tree: MergeTree = IndexMap::new();
    let mut merge_gene_trans: IndexMap<String, Vec<String>> = IndexMap::new();

    for line in read_lines(merge)? {
        let cols: Vec<&str> = line.split('\t').collect();
        let id_split: Vec<&str> = cols[3].split(';').collect();
        let trans_id = id_split[0].to_string();
        let support_id = id_split[1];
        let gene_id = trans_id.split('.').next().unwrap_or(&trans_id).to_string();
        let prefix = support_id.split('_').next().unwrap_or("").to_string();
        let collapse_id = support_id
            .split_once('_')
            .map(|x| x.1)
            .unwrap_or("")
            .to_string();

        if !merge_gene_trans.contains_key(&gene_id) {
            merge_gene_list.push(gene_id.clone());
            merge_gene_trans.insert(gene_id.clone(), Vec::new());
        }
        if !merge_gene_trans[&gene_id].contains(&trans_id) {
            merge_gene_trans
                .get_mut(&gene_id)
                .unwrap()
                .push(trans_id.clone());
        }

        let source_gene_id = collapse_id
            .split('.')
            .next()
            .unwrap_or(&collapse_id)
            .to_string();
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
            gene_trans
                .entry(gene_id)
                .or_default()
                .push(trans_id.clone());
        }
        trans_cluster
            .entry(trans_id)
            .or_default()
            .insert(cluster_id);
    }

    let mut out = tama_io::create_writer(output)?;
    writeln!(
        out,
        "gene_id\ttrans_id\tgene_num_reads\ttrans_num_reads\tcluster_line"
    )?;

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

/// Build a read-support levels file. Ports `tama_read_support_levels` (trans_read/
/// ref_anno/read_support/cluster inputs; tama/cupcake/filter merge; or no_merge).
fn levels(
    filelist: &std::path::Path,
    merge: &str,
    output_prefix: &str,
    merge_type: MergeType,
) -> anyhow::Result<()> {
    // source -> trans_id -> set(read)
    let mut src_trans_read: IndexMap<String, IndexMap<String, IndexSet<String>>> = IndexMap::new();
    let mut source_list: Vec<String> = Vec::new();
    let mut source_trans_list: Vec<String> = Vec::new();

    for line in read_lines(filelist)? {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let (source_name, transread_file, file_type) = (f[0].to_string(), f[1], f[2]);
        src_trans_read.insert(source_name.clone(), IndexMap::new());
        source_list.push(source_name.clone());
        let map = src_trans_read.get_mut(&source_name).unwrap();

        for tl in read_lines(std::path::Path::new(transread_file))? {
            match file_type {
                "trans_read" | "ref_anno" => {
                    let cols: Vec<&str> = tl.split('\t').collect();
                    let id_split: Vec<&str> = cols[3].split(';').collect();
                    let trans_id = id_split[0].to_string();
                    let read_id = id_split[1].to_string();
                    if !map.contains_key(&trans_id) {
                        source_trans_list.push(trans_id.clone());
                    }
                    map.entry(trans_id).or_default().insert(read_id);
                }
                "read_support" => {
                    if tl.starts_with("merge_gene_id") {
                        continue;
                    }
                    let cols: Vec<&str> = tl.split('\t').collect();
                    let trans_id = cols[1].to_string();
                    if !map.contains_key(&trans_id) {
                        source_trans_list.push(trans_id.clone());
                    }
                    let slot = map.entry(trans_id).or_default();
                    for src_read in cols[5].split(';') {
                        let parts: Vec<&str> = src_read.split(':').collect();
                        // reads are after the first (source) field; supports ':' in read ids
                        let read_line = parts[1..].join(":");
                        for r in read_line.split(',') {
                            slot.insert(r.to_string());
                        }
                    }
                }
                "cluster" => {
                    if tl.starts_with("cluster_id") {
                        continue;
                    }
                    let cols: Vec<&str> = tl.split(',').collect();
                    let trans_id = cols[0].to_string();
                    let read_id = cols[1].to_string();
                    if !map.contains_key(&trans_id) {
                        source_trans_list.push(trans_id.clone());
                    }
                    map.entry(trans_id).or_default().insert(read_id);
                }
                other => bail!("unknown file type {other:?}"),
            }
        }
    }

    let mut out = tama_io::create_writer(format!("{output_prefix}_read_support.txt"))?;
    writeln!(out, "merge_gene_id\tmerge_trans_id\tgene_read_count\ttrans_read_count\tsource_line\tsupport_line")?;

    if merge != "no_merge" {
        // merge_trans -> source -> source_trans -> read_list
        let mut merge_trans: IndexMap<String, IndexMap<String, IndexMap<String, Vec<String>>>> =
            IndexMap::new();
        let mut merge_trans_list: Vec<String> = Vec::new();
        let mut trans_gene: IndexMap<String, String> = IndexMap::new();
        let mut gene_read: IndexMap<String, IndexMap<String, IndexSet<String>>> = IndexMap::new();
        let has_trans_read = merge.contains("trans_read.bed");

        for line in read_lines(std::path::Path::new(merge))? {
            if merge_type == MergeType::Filter && line.starts_with("old_gene_id") {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            let (merge_trans_id, source_name, source_trans_id, merge_gene_id) = match merge_type {
                MergeType::Tama => {
                    let id_split: Vec<&str> = cols[3].split(';').collect();
                    let mtid = id_split[0].to_string();
                    let stid_line = id_split[1];
                    let (sname, stid) = if has_trans_read {
                        (source_list[0].clone(), stid_line.to_string())
                    } else {
                        let parts: Vec<&str> = stid_line.split('_').collect();
                        if parts.len() == 2 {
                            (parts[0].to_string(), parts[1].to_string())
                        } else {
                            let stid = parts[parts.len() - 1].to_string();
                            let sname = parts[..parts.len() - 1].join("_");
                            (sname, stid)
                        }
                    };
                    let mg = mtid.split('.').next().unwrap_or(&mtid).to_string();
                    (mtid, sname, stid, mg)
                }
                MergeType::Cupcake => {
                    // handled per source_trans below; loop expansion done here
                    let mtid = cols[0].to_string();
                    let parts: Vec<&str> = mtid.split('.').collect();
                    let mg = format!("{}.{}", parts[0], parts.get(1).copied().unwrap_or(""));
                    // expand comma list
                    for stid in cols[1].split(',') {
                        add_merge_entry(
                            &mut merge_trans,
                            &mut merge_trans_list,
                            &mut trans_gene,
                            &mut gene_read,
                            &src_trans_read,
                            &mtid,
                            &source_list[0],
                            stid,
                            &mg,
                        )?;
                    }
                    continue;
                }
                MergeType::Filter => {
                    let mtid = cols[5].to_string();
                    if mtid == "removed_transcript" {
                        continue;
                    }
                    (
                        mtid,
                        source_list[0].clone(),
                        cols[1].to_string(),
                        cols[4].to_string(),
                    )
                }
            };
            add_merge_entry(
                &mut merge_trans,
                &mut merge_trans_list,
                &mut trans_gene,
                &mut gene_read,
                &src_trans_read,
                &merge_trans_id,
                &source_name,
                &source_trans_id,
                &merge_gene_id,
            )?;
        }

        for mtid in &merge_trans_list {
            let mut this_reads: IndexSet<String> = IndexSet::new();
            let mut this_source_read: IndexMap<String, IndexSet<String>> = IndexMap::new();
            let mut this_source_list: Vec<String> = Vec::new();
            for source_name in &source_list {
                if let Some(stmap) = merge_trans[mtid].get(source_name) {
                    this_source_list.push(source_name.clone());
                    let slot = this_source_read.entry(source_name.clone()).or_default();
                    for reads in stmap.values() {
                        for r in reads {
                            this_reads.insert(r.clone());
                            slot.insert(r.clone());
                        }
                    }
                }
            }
            let merge_gene_id = &trans_gene[mtid];
            let mut gene_read_num = 0usize;
            let mut support_list: Vec<String> = Vec::new();
            for source_name in &source_list {
                if let Some(reads) = this_source_read.get(source_name) {
                    let line: Vec<String> = reads.iter().cloned().collect();
                    support_list.push(format!("{source_name}:{}", line.join(",")));
                }
                if let Some(gr) = gene_read
                    .get(merge_gene_id)
                    .and_then(|g| g.get(source_name))
                {
                    gene_read_num += gr.len();
                }
            }
            writeln!(
                out,
                "{merge_gene_id}\t{mtid}\t{gene_read_num}\t{}\t{}\t{}",
                this_reads.len(),
                this_source_list.join(","),
                support_list.join(";")
            )?;
        }
    } else {
        // no_merge: single source
        let source_name = &source_list[0];
        let mut gene_read: IndexMap<String, IndexSet<String>> = IndexMap::new();
        let is_cluster = source_name == "cluster";
        if !is_cluster {
            for mtid in &source_trans_list {
                let mg = mtid.split('.').next().unwrap_or(mtid).to_string();
                let slot = gene_read.entry(mg).or_default();
                for r in &src_trans_read[source_name][mtid] {
                    slot.insert(r.clone());
                }
            }
        }
        for mtid in &source_trans_list {
            let reads = &src_trans_read[source_name][mtid];
            let read_line: Vec<String> = reads.iter().cloned().collect();
            let (merge_gene_id, gene_read_num) = if is_cluster {
                ("NA".to_string(), "NA".to_string())
            } else {
                let mg = mtid.split('.').next().unwrap_or(mtid).to_string();
                let n = gene_read[&mg].len();
                (mg, n.to_string())
            };
            writeln!(
                out,
                "{merge_gene_id}\t{mtid}\t{gene_read_num}\t{}\t{source_name}\t{source_name}:{}",
                reads.len(),
                read_line.join(",")
            )?;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn add_merge_entry(
    merge_trans: &mut IndexMap<String, IndexMap<String, IndexMap<String, Vec<String>>>>,
    merge_trans_list: &mut Vec<String>,
    trans_gene: &mut IndexMap<String, String>,
    gene_read: &mut IndexMap<String, IndexMap<String, IndexSet<String>>>,
    src_trans_read: &IndexMap<String, IndexMap<String, IndexSet<String>>>,
    merge_trans_id: &str,
    source_name: &str,
    source_trans_id: &str,
    merge_gene_id: &str,
) -> anyhow::Result<()> {
    trans_gene
        .entry(merge_trans_id.to_string())
        .or_insert_with(|| merge_gene_id.to_string());
    if !merge_trans.contains_key(merge_trans_id) {
        merge_trans_list.push(merge_trans_id.to_string());
    }
    let read_list: Vec<String> = src_trans_read
        .get(source_name)
        .and_then(|t| t.get(source_trans_id))
        .map(|s| s.iter().cloned().collect())
        .unwrap_or_default();
    merge_trans
        .entry(merge_trans_id.to_string())
        .or_default()
        .entry(source_name.to_string())
        .or_default()
        .insert(source_trans_id.to_string(), read_list.clone());
    let gr = gene_read
        .entry(merge_gene_id.to_string())
        .or_default()
        .entry(source_name.to_string())
        .or_default();
    for r in read_list {
        gr.insert(r);
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
