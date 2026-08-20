#!/usr/bin/env bash
# Regenerate golden reference outputs from the original Python 2 TAMA.
#
# Requires conda/mamba. Creates a Python 2.7 + BioPython env on first run.
# The original Python sources are expected under ../../reference (gitignored;
# fetch from https://github.com/GenomeRIK/tama if missing).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
ref="$root/reference"
data="$root/test_data"

if [ ! -f "$ref/tama_collapse.py" ]; then
  echo "Missing $ref/tama_collapse.py — clone GenomeRIK/tama into $ref first." >&2
  exit 1
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | grep -q '^tama_py2 '; then
  mamba create -y -n tama_py2 -c conda-forge python=2.7 'biopython=1.76'
fi

conda run -n tama_py2 python "$ref/tama_collapse.py" \
  -s "$data/gmap_test.sam" -f "$data/test_genome.fa" -p "$here/collapse" -x capped

nocap_dir="$root/tests/golden_nocap"
mkdir -p "$nocap_dir"
conda run -n tama_py2 python "$ref/tama_collapse.py" \
  -s "$data/gmap_test.sam" -f "$data/test_genome.fa" -p "$nocap_dir/collapse" -x no_cap

# merge golden: single capped source (the capped collapse output)
merge_in="$root/test_data/merge"
merge_gold="$root/tests/golden_merge"
mkdir -p "$merge_in" "$merge_gold"
cp "$here/collapse.bed" "$merge_in/tama_collapse_test.bed"
printf 'tama_collapse_test.bed\tcapped\t1,1,1\ttestmerge\n' > "$merge_in/filelist_merge.txt"
( cd "$merge_in" && conda run -n tama_py2 python "$ref/tama_merge.py" \
    -f filelist_merge.txt -p "$merge_gold/merged" )

echo "Golden outputs regenerated in $here, $nocap_dir and $merge_gold"
