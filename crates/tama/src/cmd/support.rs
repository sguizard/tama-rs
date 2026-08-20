//! `tama support` — read-support tracking tools.

use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs)]
pub struct Args {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Read support from collapse cluster files. (tama_read_support_collapse_cluster)
    CollapseCluster,
    /// Read support levels. (tama_read_support_levels)
    Levels,
    /// Read support after merging collapse outputs. (tama_read_support_merge_collapse)
    MergeCollapse,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let name = match args.cmd {
        Cmd::CollapseCluster => "support collapse-cluster",
        Cmd::Levels => "support levels",
        Cmd::MergeCollapse => "support merge-collapse",
    };
    Err(super::not_implemented(name))
}
