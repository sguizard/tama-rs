# How this project was built

This document is a transparent account of how `tama-rs` was developed, for anyone
who wants to understand its provenance and judge its trustworthiness.

## AI-assisted development

`tama-rs` is an **AI-assisted port**. The Rust code was written with
[Claude Code](https://claude.com/claude-code) (Anthropic's Claude Opus 4.8),
directed and reviewed by a human maintainer ([@sguizard](https://github.com/sguizard))
who set the scope, made all design and release decisions, and validated the
results.

This is recorded in the git history: commits produced with AI assistance carry a
trailer, e.g.

```
Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
```

## What it is a port of

It is a from-scratch Rust reimplementation of
[GenomeRIK/tama](https://github.com/GenomeRIK/tama) (TAMA), a Python 2 toolkit
(GPL-3.0). The original source was read to reproduce each tool's behaviour; it is
**not** copied. A copy of the upstream Python is kept locally under `reference/`
purely as a porting aid and is **git-ignored** — it is never redistributed here.
`tama-rs` is licensed GPL-3.0-or-later to match upstream, and users are asked to
cite the original TAMA paper (see [`CITATION.cff`](CITATION.cff) and the README).

## Method: behavioural equivalence, verified

The guiding principle was **functional equivalence with the original**, checked
empirically rather than assumed:

1. Each tool was reimplemented in Rust (workspace crates `tama-core`,
   `tama-io`, `tama`).
2. The **original Python 2 tool** and the Rust tool were run on the **same
   input**, and their outputs were diffed ("golden" tests).
3. The Rust output is **byte-identical** to the original wherever the original is
   deterministic. The handful of exceptions are columns whose ordering derives
   from Python-2 dictionary iteration (hash order) — these are not reproducible
   and not semantically meaningful, so they are compared as sets and documented
   in code, the README, and the `CHANGELOG`.

The reference environment for regenerating golden outputs is Python 2.7 +
BioPython 1.76 (a conda/mamba env); see `tests/golden/regenerate.sh`. Golden
fixtures live under `tests/golden*`, and the Rust integration tests under
`crates/tama/tests/` diff against them.

## Reproducing the validation

```sh
# build + run the full test suite (unit + golden integration tests)
cargo test
cargo clippy --all-targets   # lints (CI runs with -D warnings)
cargo fmt --all --check      # formatting
```

To regenerate the golden reference outputs from the original Python 2 tools you
need the Python 2.7 + BioPython environment described in
`tests/golden/regenerate.sh`.

## Scope and honesty about limitations

- Not every secondary code path is ported yet; current gaps are listed in the
  README ("Status") and the `CHANGELOG` — e.g. some `collapse` options (BAM
  input, multimap handling, `low_mem` mode) and the documented Python-2
  dict-order artifacts in certain report columns.
- Where behaviour intentionally differs (for example, the toolkit is quiet by
  default rather than printing per-record progress), it is noted in the docs.

## Why this matters

Publishing the method, keeping the AI co-authorship in the commit trail, and
gating every tool behind byte-level golden tests against the original are what
make an AI-assisted port auditable. If you find a case where `tama-rs` diverges
from the original TAMA in a way that isn't a documented artifact, please open an
issue with the input and both outputs.
