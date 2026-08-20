//! Subcommand implementations, one module per original TAMA tool group.

pub mod cleanup;
pub mod collapse;
pub mod filter;
pub mod format;
pub mod merge;
pub mod orf;
pub mod split;
pub mod stats;
pub mod support;
pub mod variants;

/// Uniform error for tools that are scaffolded but not yet implemented.
pub(crate) fn not_implemented(tool: &str) -> anyhow::Error {
    anyhow::anyhow!("`{tool}` is not implemented yet in this build")
}
