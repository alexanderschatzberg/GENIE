#!/bin/bash
# Run GENIE with profiling on simulation data

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${1:-$ROOT_DIR/sim_data/n500_N500_L1}"
OUTPUT_DIR="$ROOT_DIR/profile_results"

mkdir -p "$OUTPUT_DIR"

"$ROOT_DIR/build/GENIE" \
  -g "$DATA_DIR/genotypes" \
  -p "$DATA_DIR/phenotype.csv" \
  -c "$DATA_DIR/covariates.csv" \
  -e "$DATA_DIR/environment.csv" \
  -m G+GxE \
  -k 10 \
  -jn 10 \
  -o "$OUTPUT_DIR/results.out" \
  --profile \
  --profile_out "$OUTPUT_DIR/timing.csv"

echo "Results written to $OUTPUT_DIR/"
