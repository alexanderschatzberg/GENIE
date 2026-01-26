#!/usr/bin/env python3
"""Plot thread scaling results from GENIE profiling."""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import argparse
import re


def load_thread_data(results_dir: Path) -> pd.DataFrame:
    """Load timing CSV files and extract thread counts from filenames."""
    records = []
    for csv_file in sorted(results_dir.glob("*_timing.csv")):
        # Extract thread count from filename (e.g., n500_N500_L1_t4_timing.csv)
        match = re.search(r"_t(\d+)_timing\.csv$", csv_file.name)
        if not match:
            continue
        threads = int(match.group(1))

        df = pd.read_csv(csv_file)
        df["threads"] = threads
        records.append(df)

    if not records:
        return pd.DataFrame()
    return pd.concat(records, ignore_index=True)


def plot_thread_scaling(df: pd.DataFrame, output_path: Path = None):
    """Create thread scaling plots."""
    df = df[df["name"].notna() & (df["name"] != "")]

    threads = sorted(df["threads"].unique())
    operations = sorted(df["name"].unique())

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    colors = plt.cm.tab10(np.linspace(0, 1, len(operations)))

    # Plot 1: Absolute time vs threads
    ax1 = axes[0]
    for i, op in enumerate(operations):
        times = []
        for t in threads:
            mask = (df["threads"] == t) & (df["name"] == op)
            if mask.any():
                times.append(df.loc[mask, "total_seconds"].values[0])
            else:
                times.append(np.nan)
        ax1.plot(threads, times, "o-", label=op, color=colors[i])

    ax1.set_xlabel("Threads")
    ax1.set_ylabel("Time (seconds)")
    ax1.set_title("Absolute Time vs Thread Count")
    ax1.set_xscale("log", base=2)
    ax1.set_xticks(threads)
    ax1.set_xticklabels(threads)
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Speedup vs threads (relative to single-threaded)
    ax2 = axes[1]
    for i, op in enumerate(operations):
        times = []
        for t in threads:
            mask = (df["threads"] == t) & (df["name"] == op)
            if mask.any():
                times.append(df.loc[mask, "total_seconds"].values[0])
            else:
                times.append(np.nan)

        if times and times[0] > 0:
            speedup = [times[0] / t if t > 0 else np.nan for t in times]
            ax2.plot(threads, speedup, "o-", label=op, color=colors[i])

    # Plot ideal scaling line
    ax2.plot(threads, threads, "k--", label="Ideal", alpha=0.5)

    ax2.set_xlabel("Threads")
    ax2.set_ylabel("Speedup (T1 / Tn)")
    ax2.set_title("Speedup vs Thread Count")
    ax2.set_xscale("log", base=2)
    ax2.set_xticks(threads)
    ax2.set_xticklabels(threads)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved plot to {output_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot GENIE thread scaling results")
    parser.add_argument(
        "--input", "-i",
        type=Path,
        default=Path(__file__).parent.parent / "profile_threads_results",
        help="Directory containing timing CSV files"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Output file path (e.g., plot.png). If not specified, displays plot."
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: {args.input} does not exist")
        return 1

    df = load_thread_data(args.input)
    if df.empty:
        print(f"Error: No thread timing CSV files found in {args.input}")
        return 1

    plot_thread_scaling(df, args.output)
    return 0


if __name__ == "__main__":
    exit(main())
