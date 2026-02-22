#!/bin/bash
# Profile GENIE on a single dataset at various thread counts

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <data_dir> [output_dir]"
    echo "Example: $0 sim_data/n500_N500_L1"
    exit 1
fi

DATA_DIR="$1"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${2:-$ROOT_DIR/profile_threads_results}"
THREADS=(1 2 4 8 16 32 64 128 256 512)

# Extract dataset name from path
NAME=$(basename "$DATA_DIR")

mkdir -p "$OUTPUT_DIR"

echo "Profiling $NAME at thread counts: ${THREADS[*]}"
echo "Results will be written to $OUTPUT_DIR"
echo ""

for T in "${THREADS[@]}"; do
    echo "=== Running with $T thread(s) ==="

    "$ROOT_DIR/build/GENIE" \
        -g "$DATA_DIR/genotypes" \
        -p "$DATA_DIR/phenotype.csv" \
        -c "$DATA_DIR/covariates.csv" \
        -e "$DATA_DIR/environment.csv" \
        -m G+GxE \
        -k 10 \
        -jn 100 \
        -t "$T" \
        -o "$OUTPUT_DIR/${NAME}_t${T}_results.out" \
        --profile \
        --profile_out "$OUTPUT_DIR/${NAME}_t${T}_timing.csv"

    echo ""
done

echo "Done! Results written to $OUTPUT_DIR/"
echo ""
echo "To plot results:"
echo "  python scripts/plot_thread_scaling.py -i $OUTPUT_DIR"
