# GENIE: Getting Started and Profiling Guide

GENIE (**G**ene-**EN**vironment **I**nteraction **E**stimator) is a C++ tool for
estimating heritability components including additive genetic (G),
gene-environment interaction (GxE), and heterogeneous noise (NxE) effects.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Building GENIE](#building-genie)
3. [Quick Start](#quick-start)
4. [Data Format Reference](#data-format-reference)
5. [Running with Your Own Data](#running-with-your-own-data)
6. [Profiling](#profiling)
   - [Building with Profiling](#building-with-profiling)
   - [Single Dataset Profiling](#single-dataset-profiling)
   - [Batch Profiling Across Datasets](#batch-profiling-across-datasets)
   - [Thread Scaling Analysis](#thread-scaling-analysis)
   - [Profiled Regions](#profiled-regions)
   - [Interpreting Results](#interpreting-results)
   - [Plotting Results](#plotting-results)
7. [Model Specifications](#model-specifications)
8. [Troubleshooting](#troubleshooting)

---

## Prerequisites

- C++ compiler with C++11 support (GCC 5+, Clang 5+)
- CMake 3.10+
- Eigen3 (header-only, included or installed via package manager)
- Python 3.7+ with `pandas`, `matplotlib`, `numpy` (for plotting scripts)

## Building GENIE

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The `GENIE` executable will be created in the `build/` directory.

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `-DCMAKE_BUILD_TYPE=Release` | Release | Build type (Debug, Release, RelWithDebInfo, MinSizeRel) |
| `-DENABLE_SSE=ON` | OFF | SSE intrinsics for x86_64 platforms |
| `-DUSE_DOUBLE=ON` | OFF | Double precision instead of single precision |
| `-DGENIE_PROFILE=ON` | OFF | Enable internal profiling instrumentation |

## Quick Start

Run the included toy example to verify your build:

```bash
cd example
chmod +x test.sh
./test.sh
```

This runs GENIE on a small bundled dataset with the `G+GxE+NxE` model.

## Data Format Reference

GENIE requires input files in the following formats:

**Genotype files** — PLINK binary format (BED/BIM/FAM). Pass the prefix
(without the `.bed`/`.bim`/`.fam` extension) to the `-g` flag.

**Phenotype / Covariate / Environment files** — Space or tab-delimited text
with a header row:

```
FID IID var1 var2 ... varN
1001 1001 0.52 1.3
1002 1002 -0.71 0.8
```

**Annotation files** — Space-delimited M x K matrix with **no header**. Each
row corresponds to a SNP (in the same order as the BIM file) and each column
is a binary annotation indicator (1 = SNP belongs to annotation, 0 = not).

### Important ordering rules

- Phenotype, covariate, environment, and genotype files must list individuals
  in the **same order**.
- Annotation rows must match the BIM file SNP order.
- SNPs with MAF = 0 must be removed before running GENIE.
- Individuals with `NA` or `-9` in phenotype/environment are automatically
  excluded.

## Running with Your Own Data

### Basic run (additive genetic effects only)

```bash
./build/GENIE \
    -g path/to/genotypes \
    -p path/to/phenotype.txt \
    -c path/to/covariates.txt \
    -m G \
    -k 10 -jn 100 \
    -o results.out
```

### Gene-environment interaction model

```bash
./build/GENIE \
    -g path/to/genotypes \
    -p path/to/phenotype.txt \
    -c path/to/covariates.txt \
    -e path/to/environment.txt \
    -m G+GxE \
    -k 10 -jn 100 \
    -o results.out
```

### Full model with heterogeneous noise

```bash
./build/GENIE \
    -g path/to/genotypes \
    -p path/to/phenotype.txt \
    -c path/to/covariates.txt \
    -e path/to/environment.txt \
    --annot path/to/annotations.txt \
    -m G+GxE+NxE \
    -k 10 -jn 100 \
    -o results.out
```

### Key parameters

| Flag | Long form | Description |
|------|-----------|-------------|
| `-g` | `--genotype` | PLINK BED prefix |
| `-p` | `--phenotype` | Phenotype file |
| `-c` | `--covariate` | Covariate file |
| `-e` | `--environment` | Environment variable file |
| `--annot` | | Annotation file |
| `-m` | `--model` | Model: `G`, `G+GxE`, `G+GxE+NxE` |
| `-k` | | Number of random vectors for trace estimation (default: 10) |
| `-jn` | `--num-jack` | Number of jackknife blocks (default: 100) |
| `-t` | `--threads` | Number of threads |
| `-o` | `--output` | Output file path |
| `--memeff` | | Memory-efficient mode (process SNPs in blocks) |

---

## Profiling

GENIE includes a built-in profiler that records wall-clock timing and call
counts for key computational regions. When disabled at compile time, all
profiling code compiles to zero-overhead no-ops.

### Building with Profiling

```bash
mkdir build && cd build
cmake -DGENIE_PROFILE=ON ..
make -j$(nproc)
```

### Single Dataset Profiling

Run GENIE with `--profile` and `--profile_out` to capture timing data:

```bash
./build/GENIE \
    -g path/to/genotypes \
    -p path/to/phenotype.txt \
    -c path/to/covariates.txt \
    -e path/to/environment.txt \
    -m G+GxE \
    -k 10 -jn 10 \
    --profile \
    --profile_out timing_output.csv
```

The output format is determined by the file extension: `.csv` for CSV,
`.json` for JSON.

**CSV output example:**

```csv
name,calls,total_seconds,avg_seconds,total_bytes
jackknife_block,10,28.98,2.898,0
matvec_Xv,600,15.23,0.0254,0
matvec_Xt_v,400,12.11,0.0303,0
solve_normal_equations,11,0.003,0.0003,0
io_read_genotype,1,0.45,0.45,0
```

### Batch Profiling Across Datasets

Use the provided script to profile GENIE across all datasets in a directory.
Each subdirectory should contain genotype, phenotype, covariate, and
environment files.

```bash
# Profile all datasets under sim_data/
scripts/run_profiler.sh

# Results are written to profile_results/
ls profile_results/
#   n100_N100_L1_timing.csv
#   n200_N200_L1_timing.csv
#   ...
```

**Expected directory layout for each dataset:**

```
sim_data/
  my_dataset/
    genotypes.bed
    genotypes.bim
    genotypes.fam
    phenotype.csv
    covariates.csv
    environment.csv
```

To profile datasets in a custom location, edit `run_profiler.sh` and change
the `DATA_DIR` glob or pass arguments directly.

### Thread Scaling Analysis

Measure how GENIE scales across different thread counts on a single dataset:

```bash
# Profile at 1, 2, 4, 8, 16, and 32 threads
scripts/profile_threads.sh sim_data/n500_N500_L1

# Or specify a custom output directory
scripts/profile_threads.sh sim_data/n500_N500_L1 my_thread_results/
```

Output files follow the naming convention `<dataset>_t<threads>_timing.csv`.

### Profiled Regions

The profiler instruments the following regions across all GENIE executables:

| Region | Description |
|--------|-------------|
| `matvec_Xv` | Genotype matrix x vector multiplication (X * v) |
| `matvec_Xt_v` | Transposed genotype matrix x vector (X^T * v) |
| `jackknife_pass1` | First pass over jackknife blocks |
| `jackknife_pass2` | Second pass over jackknife blocks |
| `jackknife_block` | Processing of individual jackknife blocks |
| `trace_assembly` | Trace matrix assembly for variance component estimation |
| `solve_normal_equations` | Solving the normal equations for heritability |
| `regress_covariates` | Regressing covariates from phenotypes |
| `io_read_genotype` | Reading PLINK BED/BIM/FAM files |
| `io_read_bed_block` | Reading a block of BED data (memory-efficient and multi-pheno modes) |
| `io_read_phenotype` | Reading the phenotype file |
| `io_read_covariate` | Reading the covariate file |
| `io_read_environment` | Reading the environment file |
| `io_read_annotation` | Reading the annotation file |
| `compute_yXXy` | Computing y'X X'y quadratic forms |
| `compute_XXz` | Computing X X'z products |
| `random_vector_init` | Initializing random vectors for trace estimation |
| `matmult_init` | Initializing matrix multiplication structures |

Not all regions appear in every run — it depends on the model and executable
used.

### Interpreting Results

Results are sorted by total time descending, so hotspots appear first.

- **Matrix operations** (`matvec_Xv`, `matvec_Xt_v`) typically dominate for
  large datasets and should scale with N x M (individuals x SNPs).
- **Jackknife loops** (`jackknife_block`, `jackknife_pass1/2`) wrap the
  matrix operations and should account for most wall-clock time.
- **Trace assembly** (`trace_assembly`) and **solver**
  (`solve_normal_equations`) are usually fast relative to matrix operations.
- **I/O** regions should be a small fraction of total time. If I/O dominates,
  consider using faster storage or the memory-efficient mode (`--memeff`).
- **Thread scaling**: compare `total_seconds` for `matvec_Xv` across
  different thread counts. Diminishing returns above the number of physical
  cores is expected.

### Plotting Results

Two plotting scripts are provided in `scripts/`:

**1. Compare operations across datasets:**

```bash
# Plot from the default profile_results/ directory
python scripts/plot_profiling.py

# Specify a custom directory and save to file
python scripts/plot_profiling.py -i my_results/ -o profiling_comparison.png

# Filter to specific operations
python scripts/plot_profiling.py -i my_results/ --filter matvec_Xv,matvec_Xt_v
```

**2. Thread scaling analysis:**

```bash
# Plot from the default profile_threads_results/ directory
python scripts/plot_thread_scaling.py

# Specify a custom directory and dataset filter
python scripts/plot_thread_scaling.py -i my_thread_results/ -o scaling.png

# Filter to specific operations
python scripts/plot_thread_scaling.py -i my_thread_results/ --filter matvec_Xv,jackknife_block
```

Both scripts support `--filter` to focus on specific operations of interest.

---

## Model Specifications

| Model | Flag | Components estimated |
|-------|------|---------------------|
| Additive genetic only | `-m G` | Additive genetic variance (h^2) per annotation |
| GxE interaction | `-m G+GxE` | Additive + gene-environment interaction |
| Full heterogeneous noise | `-m G+GxE+NxE` | Additive + GxE + heterogeneous noise |

## Troubleshooting

**"profiling not compiled in" warning** — Rebuild with
`cmake -DGENIE_PROFILE=ON ..`.

**Empty profiling output** — Ensure GENIE completed successfully. Check
stderr for errors.

**Some regions show 0 calls** — Those code paths were not exercised by your
particular model/input combination.

**Plotting scripts fail with ImportError** — Install dependencies:
`pip install pandas matplotlib numpy`.

**Thread scaling shows no improvement** — Ensure your dataset is large
enough. Small datasets (< 1000 individuals) have insufficient work to
benefit from parallelism. Also check that OpenMP / threading is enabled in
your build.
