#!/bin/bash
# Run GENIE with profiling on all simulation datasets

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="$ROOT_DIR/profile_results"

mkdir -p "$OUTPUT_DIR"

for DATA_DIR in "$ROOT_DIR"/sim_data/*/; do
  NAME=$(basename "$DATA_DIR")
  echo "=== Running profiler on $NAME ==="

  "$ROOT_DIR/build/GENIE" \
    -g "$DATA_DIR/genotypes" \
    -p "$DATA_DIR/phenotype.csv" \
    -c "$DATA_DIR/covariates.csv" \
    -e "$DATA_DIR/environment.csv" \
    -m G+GxE \
    -k 10 \
    -jn 10 \
    -o "$OUTPUT_DIR/${NAME}_results.out" \
    --profile \
    --profile_out "$OUTPUT_DIR/${NAME}_timing.csv"

  echo ""
done

echo "Results written to $OUTPUT_DIR/"
