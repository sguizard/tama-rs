#!/usr/bin/env bash
# Regenerate the committed shell completion scripts from the clap definitions.
#
# Run this after changing anything about the CLI surface (subcommands, flags, or
# the value enums in crates/tama/src/cmd/opts.rs). CI checks that the committed
# copy matches what the current binary emits.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"

cargo run --quiet --manifest-path "$root/Cargo.toml" -p tama -- completions fish \
  > "$here/tama.fish"

echo "Regenerated $here/tama.fish"
