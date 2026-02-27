#!/usr/bin/env python3
"""Plot profiling results from GENIE runs.

Reads timing CSV files from a results directory and creates a grouped bar
plot comparing operations across datasets.

Usage:
    python plot_profiling.py                            # defaults
    python plot_profiling.py -i my_results/ -o plot.png
    python plot_profiling.py --filter matvec_Xv,matvec_Xt_v
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import argparse

# Colormap for plots (change this to customize colors)
COLORMAP = plt.cm.tab20
OPACITY = 0.95


def load_profiling_data(results_dir: Path) -> pd.DataFrame:
    """Load all timing CSV files from the results directory."""
    records = []
    for csv_file in sorted(results_dir.glob("*_timing.csv")):
        dataset = csv_file.stem.replace("_timing", "")
        df = pd.read_csv(csv_file)
        df["dataset"] = dataset
        records.append(df)
    if not records:
        return pd.DataFrame()
    return pd.concat(records, ignore_index=True)


def plot_profiling(df: pd.DataFrame, output_path: Path = None,
                   filter_ops: list = None):
    """Create a grouped bar plot of profiling results."""
    # Filter out empty rows
    df = df[df["name"].notna() & (df["name"] != "")]

    # Apply operation filter if specified
    if filter_ops:
        df = df[df["name"].isin(filter_ops)]
        if df.empty:
            print(f"Error: no matching operations for filter {filter_ops}")
            return

    datasets = sorted(df["dataset"].unique())
    operations = sorted(df["name"].unique())

    x = np.arange(len(datasets))
    width = 0.8 / max(len(operations), 1)

    fig, ax = plt.subplots(figsize=(max(10, len(datasets) * 1.5), 6))
    colors = [(*c[:3], OPACITY)
              for c in COLORMAP(np.linspace(0, 1, len(operations)))]

    for i, op in enumerate(operations):
        times = []
        for dataset in datasets:
            mask = (df["dataset"] == dataset) & (df["name"] == op)
            if mask.any():
                times.append(df.loc[mask, "self_seconds"].values[0])
            else:
                times.append(0)

        offset = (i - len(operations) / 2 + 0.5) * width
        ax.bar(x + offset, times, width, label=op, color=colors[i])

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Self-Time (seconds)")
    ax.set_title("GENIE Profiling Results by Dataset")
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, rotation=30, ha="right")
    ax.legend(loc="upper left", fontsize=7, ncol=max(1, len(operations) // 8))
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved plot to {output_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot GENIE profiling results")
    parser.add_argument(
        "--input", "-i",
        type=Path,
        default=Path(__file__).parent.parent / "profile_results",
        help="Directory containing timing CSV files"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Output file path (e.g., plot.png). If not specified, displays plot."
    )
    parser.add_argument(
        "--filter", "-f",
        type=str,
        default=None,
        help="Comma-separated list of operation names to include (e.g., "
             "Xv,Xt_v,jackknife)"
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: {args.input} does not exist")
        return 1

    df = load_profiling_data(args.input)
    if df.empty:
        print(f"Error: No timing CSV files found in {args.input}")
        return 1

    filter_ops = None
    if args.filter:
        filter_ops = [s.strip() for s in args.filter.split(",")]

    plot_profiling(df, args.output, filter_ops)
    return 0


if __name__ == "__main__":
    exit(main())
