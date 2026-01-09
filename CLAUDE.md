# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

GENIE (**G**ene-**EN**viroment **I**nteraction **E**stimator) is a C++ tool for estimating heritability components including additive genetic (G), gene-environment (GxE), and heterogeneous noise (NxE) effects. The repository also contains FLEX, a related tool for LD-aware estimation of gene-level heritability.

## Build and Development Commands

### Building GENIE

```bash
mkdir build
cd build
cmake ..
make
```

The main executable `GENIE` will be in the `build/` directory.

### Build Options

CMake supports several configuration options:

- **Build type**: `-DCMAKE_BUILD_TYPE=<Debug|Release|MinSizeRel|RelWithDebInfo>` (default: Release)
- **SSE support**: `-DENABLE_SSE=ON` (only for amd64 platforms, default: OFF)
- **Precision**: `-DUSE_DOUBLE=ON` (use double precision instead of float, default: OFF)

### Debug/Sanitizer Flags

For development and debugging:

```bash
cmake .. -DENABLE_ASAN=ON       # AddressSanitizer
cmake .. -DENABLE_UBSAN=ON      # UndefinedBehaviorSanitizer
cmake .. -DENABLE_TSAN=ON       # ThreadSanitizer (mutually exclusive with ASAN/UBSAN)
cmake .. -DENABLE_FP_TRAPS=ON   # Floating-point traps
cmake .. -DENABLE_ASSERTIONS=ON # libstdc++ debug assertions
cmake .. -DENABLE_EIGEN_NAN_INIT=ON # Initialize Eigen matrices to NaN
```

### Building FLEX

```bash
cd FLEX
mkdir build && cd build
cmake ..
make
```

The FLEX executable `FLEX_h2` will be in `FLEX/build/`.

### Testing

Run the toy example to verify everything works:

```bash
cd example
chmod +x test.sh
./test.sh
```

For Python tests:

```bash
python test/test.py
```

### Python Package

GENIE can be built as a Python package (`rhe`) using scikit-build-core. The package exposes Python wrappers for running GENIE from Python scripts.

## Architecture

### Core Components

**Genotype Representation (Mailman)**
- `include/mailman.h`, `include/storage.h`, `src/storage.cpp`
- GENIE uses a memory-efficient genotype encoding called "Mailman" that compresses genotype data
- SNPs are stored in a packed 2-bit format in `storage.h`
- Fast matrix multiplication routines in `mailman.h` operate directly on compressed data
- SSE-optimized variants available for x86_64 when `-DENABLE_SSE=ON`

**Matrix Operations**
- `include/matmult.h`, `src/matmult.cpp`
- `include/genomult.h`, `src/genomult.cpp`
- High-performance genotype-matrix multiplication
- Handles both normal and SSE-optimized code paths based on `SSE_SUPPORT` macro

**Genotype I/O**
- `include/genotype.h`, `src/genotype.cpp`
- Reads PLINK BED/BIM/FAM format files
- Manages SNP metadata and sample information

**Argument Parsing and Configuration**
- `include/arguments.h`
- Defines `options` struct containing all command-line parameters
- Supports both command-line flags and config file input
- Config files use format: `key=value` (one per line)

**Utility Functions**
- `include/auxillary.h`, `src/auxillary.cpp`: Matrix operations, jackknife utilities
- `include/functions.h`, `src/functions.cpp`: General utility functions
- `include/vectorfn.h`, `include/statsfn.h`: Vector and statistical functions
- `include/helper.h`: Helper macros and utilities
- `include/io.h`: I/O utilities and debugging support

### Main Executables

**GENIE Main** (`src/ge_flexible.cpp`)
- Primary entry point for standard GENIE analysis
- Supports models: `G`, `G+GxE`, `G+GxE+NxE`
- Uses stochastic trace estimation with random vectors (`-k` parameter)
- Jackknife-based variance estimation (`-jn` parameter)

**GENIE Variants**
- `ge_flexible_multi_pheno.cpp`: Multiple phenotype support
- `ge_hetro_flexible.cpp`, `ge_hetro_flexible_multi_pheno.cpp`: Heterogeneous noise models
- `ge_mem_flexible.cpp`: Memory-efficient version for large annotation files
- `ge_hetro_flexible_multi_pheno_trace_step1.cpp`, `ge_hetro_flexible_multi_pheno_trace_step2.cpp`: Two-stage trace estimation

### FLEX Subproject

FLEX has a similar architecture but focused on gene-level heritability:

**C++ Core**
- `FLEX/src/flex_h2.cpp`: Main binary for individual-level heritability
- Shares genotype infrastructure with GENIE (`genotype.h`, `storage.h`, `matmult.h`)

**Python Modules**
- `FLEX/src/flex_cond_test.py`: Conditional testing module
- `FLEX/src/flex_summ_h2/`: Summary statistics mode (two-stage LD score computation)

## Data Flow and Key Concepts

### Memory-Efficient Mode
- Enable with `-mem` or `--memeff` flag
- Processes SNPs in blocks to reduce memory footprint
- Block size controlled by `-mem_Nsnp` (default: 10 SNPs)
- Alternative: `--opt2 0` reads entire jackknife blocks at once

### Jackknife Schemes
Controlled by `jack_scheme` option (internal):
- 0: No jackknife
- 1: Constant number of SNPs per block (default)
- 2: Adjusting for chromosome boundaries
- 3: Constant physical length with chromosome boundaries

### Annotation System
- Annotation files are M×K matrices (M=SNPs, K=annotations)
- Binary indicators: 1 if SNP belongs to annotation, 0 otherwise
- By default, GENIE fits a single GxE component; use `-eXannot` to partition GxE by annotations
- In FLEX, the first annotation serves as the flanking/LD region for conditioning

### Environment Variables
- Environment files specify variables for GxE interactions
- By default, environment variables are added to covariates as fixed effects
- Must match individual order in genotype/phenotype files

## Precision Control

GENIE supports both single and double precision via the `USE_DOUBLE` compile-time flag:
- When `USE_DOUBLE=0` (default): Uses `float` and `MatrixXf`
- When `USE_DOUBLE=1`: Uses `double` and `MatrixXd`
- The `MatrixXdr` typedef in code adapts based on this setting

## Common Gotchas

- **SNP ordering**: Annotation file must have rows in same order as BIM file
- **Individual ordering**: Phenotype, genotype, covariate, and environment files must have individuals in identical order
- **MAF filtering**: SNPs with MAF=0 must be removed before running GENIE
- **Missing values**: Individuals with NA or -9 in phenotype/environment are automatically excluded
- **No headers in annotation files**: Unlike phenotype/covariate files, annotation files have no header row
- **Config files use full keys**: When using `--config`, use full parameter names (e.g., `genotype=`) not shortcuts (e.g., `-g`)

## File Format Requirements

**Phenotype/Covariate/Environment**: Space/tab-delimited with header:
```
FID IID var1 var2 ... varN
```

**Annotation**: Space-delimited M×K matrix with NO header
- M rows (one per SNP in BIM file order)
- K columns (one per annotation)

**Genotype**: PLINK BED/BIM/FAM format (use prefix for `-g` flag)

## Model Specifications

Use `-m` or `--model` flag:
- `G`: Additive genetic effects only (equivalent to RHE-mc)
- `G+GxE`: Additive genetic + gene-environment interactions (homogeneous noise)
- `G+GxE+NxE`: Full model with heterogeneous noise effects

## Python Interface

The `rhe` Python package provides wrappers:
```python
from rhe import run_genie, run_genie_mem, run_genie_multi_pheno

# Arguments can be list or string
run_genie(['-g', 'genotype', '-p', 'pheno', ...])
run_genie('-g genotype -p pheno ...')
```
