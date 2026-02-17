# Profiling GENIE

Step-by-step instructions for profiling GENIE on your datasets.

## 1. Build with Profiling Enabled

```bash
mkdir build && cd build
cmake -DGENIE_PROFILE=ON ..
make -j$(nproc)
```

When built without `-DGENIE_PROFILE=ON`, all profiling code compiles to
zero-overhead no-ops. Using `--profile` at runtime without the build flag
will print a warning and disable profiling.

## 2. Profile a Single Dataset

Add `--profile` and optionally `--profile_out <file>` to your normal GENIE
command:

```bash
./build/GENIE \
    -g path/to/genotypes \
    -p path/to/phenotype.txt \
    -c path/to/covariates.txt \
    -e path/to/environment.txt \
    -m G+GxE \
    -k 10 -jn 10 \
    --profile \
    --profile_out timing.csv
```

The output format is determined by the file extension:
- `.csv` — comma-separated (default: `genie_profile.csv`)
- `.json` — JSON array

### CSV output columns

| Column | Description |
|--------|-------------|
| `name` | Profiled region name |
| `calls` | Number of times the region was entered |
| `total_seconds` | Cumulative wall-clock time (6 decimal places) |
| `avg_seconds` | `total_seconds / calls` |
| `total_bytes` | Bytes transferred (0 if not tracked for that region) |

Rows are sorted by `total_seconds` descending so hotspots appear first.

### Example CSV output

```csv
name,calls,total_seconds,avg_seconds,total_bytes
jackknife_block,10,28.980000,2.898000,0
jackknife_pass1,1,18.450000,18.450000,0
matvec_Xv,600,15.230000,0.025383,0
matvec_Xt_v,400,12.110000,0.030275,0
jackknife_pass2,1,10.520000,10.520000,0
compute_yXXy,200,8.340000,0.041700,0
compute_XXz,400,7.920000,0.019800,0
trace_assembly,11,1.230000,0.111818,0
regress_covariates,1,0.520000,0.520000,0
io_read_genotype,1,0.450000,0.450000,0
io_read_bed_block,10,0.380000,0.038000,0
random_vector_init,1,0.120000,0.120000,0
matmult_init,1,0.085000,0.085000,0
io_read_phenotype,1,0.023000,0.023000,0
io_read_covariate,1,0.018000,0.018000,0
io_read_environment,1,0.015000,0.015000,0
io_read_annotation,1,0.012000,0.012000,0
solve_normal_equations,11,0.003000,0.000273,0
```

Not all regions appear in every run — which regions are active depends on
the model (`-m`) and executable used.

## 3. Profile Multiple Datasets

Use the batch script to profile every dataset in a directory:

```bash
scripts/run_profiler.sh
```

Each subdirectory under `sim_data/` should contain genotype (BED/BIM/FAM),
phenotype, covariate, and environment files. Results are written to
`profile_results/`.

To use a different data directory, edit `DATA_DIR` in `run_profiler.sh`.

## 4. Thread Scaling Analysis

Measure how GENIE scales across thread counts on a single dataset:

```bash
scripts/profile_threads.sh sim_data/n500_N500_L1

# Or with a custom output directory
scripts/profile_threads.sh sim_data/n500_N500_L1 my_thread_results/
```

This runs GENIE at 1, 2, 4, 8, 16, and 32 threads. Output files are named
`<dataset>_t<threads>_timing.csv`.

## 5. Plotting

**Compare operations across datasets:**

```bash
python scripts/plot_profiling.py
python scripts/plot_profiling.py -i my_results/ -o comparison.png
python scripts/plot_profiling.py -i my_results/ --filter matvec_Xv,matvec_Xt_v
```

**Thread scaling plot:**

```bash
python scripts/plot_thread_scaling.py
python scripts/plot_thread_scaling.py -i my_thread_results/ -o scaling.png
python scripts/plot_thread_scaling.py --filter matvec_Xv,jackknife_block
```

Both scripts require `pandas`, `matplotlib`, and `numpy`.

## 6. Interpreting Results

- **Matrix operations** (`matvec_Xv`, `matvec_Xt_v`, `compute_yXXy`,
  `compute_XXz`) typically dominate and scale with N x M.
- **Jackknife loops** (`jackknife_block`, `jackknife_pass1/2`) wrap the
  matrix operations and account for most wall-clock time.
- **Trace assembly** and **solve** are usually fast relative to matrix ops.
- **I/O regions** should be a small fraction of total time. If I/O
  dominates, try faster storage or `--memeff` mode.
- **Thread scaling**: compare `total_seconds` for `matvec_Xv` across thread
  counts. Diminishing returns above the number of physical cores is expected.

## Profiled Regions Reference

| Region | Description |
|--------|-------------|
| `matvec_Xv` | Genotype matrix x vector (X * v) |
| `matvec_Xt_v` | Transposed genotype matrix x vector (X^T * v) |
| `compute_yXXy` | y'X X'y quadratic forms |
| `compute_XXz` | X X'z products |
| `jackknife_pass1` | First pass over jackknife blocks |
| `jackknife_pass2` | Second pass over jackknife blocks |
| `jackknife_block` | Individual jackknife block processing |
| `trace_assembly` | Trace matrix assembly for variance components |
| `solve_normal_equations` | Solving normal equations for heritability |
| `regress_covariates` | Regressing covariates from phenotypes |
| `random_vector_init` | Initializing random vectors for trace estimation |
| `matmult_init` | Initializing matrix multiplication structures |
| `io_read_genotype` | Reading PLINK BED/BIM/FAM files |
| `io_read_bed_block` | Reading a BED block (memory-efficient / multi-pheno modes) |
| `io_read_phenotype` | Reading phenotype file |
| `io_read_covariate` | Reading covariate file |
| `io_read_environment` | Reading environment file |
| `io_read_annotation` | Reading annotation file |

## Troubleshooting

- **"profiling not compiled in" warning** — Rebuild with `cmake -DGENIE_PROFILE=ON ..`.
- **Empty profiling output** — Check stderr; GENIE may have exited early.
- **Missing regions** — That code path wasn't exercised by your model/input.
- **No thread scaling improvement** — Dataset may be too small (< 1000 individuals), or OpenMP may not be enabled.
- **Plotting ImportError** — `pip install pandas matplotlib numpy`.
