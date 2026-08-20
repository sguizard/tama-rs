# Benchmarks

A comparison of `tama collapse` — the compute-heavy core tool — between the
original Python 2 implementation and this Rust port. Both produce
**byte-identical** output at every input size tested, so this measures the cost
of the same work, not a different result.

## Results

`tama collapse -x capped`, median wall time and peak resident memory
(`/usr/bin/time`):

| Reads  | Python 2 time | Python 2 mem | Rust time | Rust mem | Speedup | Mem ratio |
|-------:|--------------:|-------------:|----------:|---------:|:-------:|:---------:|
| 104    | 0.26 s        | 26 MB        | 0.01 s    | 5.5 MB   | **26×** | 4.8×      |
| 1,040  | 1.87 s        | 81 MB        | 0.05 s    | 13 MB    | **37×** | 6.1×      |
| 10,400 | 88.3 s        | 802 MB       | 1.46 s    | 84 MB    | **60×** | 9.5×      |
| 52,000 | aborted¹      | —            | 26.3 s    | 424 MB   | —       | —         |

¹ The Python 500× run was stopped after exceeding an hour; extrapolating its
superlinear curve it would take well over that. Rust finished the same input in
26 seconds.

### Observations

- **~26–60× faster**, and the gap widens with input size.
- **~5–10× less memory**, also widening.
- Both implementations scale **superlinearly** on this input, but Python much more
  steeply: going from 1,040 → 10,400 reads (10× the data) costs Python ~47× more
  time versus ~29× for Rust. The superlinearity comes from the collapse step's
  pairwise per-locus comparison (see the caveat below).
- **Correctness is preserved under load:** the `.bed` output was diffed against the
  Python output at every size and is identical.

## Methodology

- **Command:** `tama collapse -s <input.sam> -f test_genome.fa -p <out> -x capped`
  for both implementations.
- **Input:** the repository's real `test_data/gmap_test.sam` (104 mapped reads)
  duplicated N× — each copy gets a unique read ID and the file stays
  coordinate-sorted — against the single-scaffold `test_data/test_genome.fa`.
- **Timing:** `/usr/bin/time -f "%e %M"`; the reported value is the **median** of
  repeated runs (7 runs at 1×/10×, 3 at 100×, 2 at 500×).
- **Equivalence check:** the Rust and Python `.bed` outputs were compared with
  `diff` at each size.

### Caveat: synthetic worst case

Duplicating the same 104 reads concentrates **all** reads onto the same handful of
genomic loci. That is a deliberate stress test of collapse's per-locus pairwise
comparison — the part that scales superlinearly — rather than a natural read
distribution. Real long-read data spreads reads across thousands of genes, so
absolute times will be lower and scaling flatter for both tools. The **relative**
speedup and memory advantage shown here are representative; the absolute
superlinear curve is specific to this concentrated input.

The lighter TAMA GO tools (format converters, filters, splitters) are I/O-bound;
the Rust versions additionally avoid the Python interpreter and BioPython import
on every invocation, so they start and finish in milliseconds.

## Environment

| | |
|---|---|
| CPU | Intel Xeon Silver 4514Y |
| OS | Rocky Linux 9.7 |
| Rust | rustc 1.98.0, `--release` build |
| Python | 2.7.15 + BioPython 1.76 (conda env) |

## Reproducing

The original Python 2 tools and a Python 2.7 + BioPython environment are required
(see `tests/golden/regenerate.sh`). Scale the bundled test SAM to the desired
size, then time both implementations with the same arguments. Because output is
byte-identical, a `diff` of the two `.bed` files is a good correctness gate before
trusting any timing.
