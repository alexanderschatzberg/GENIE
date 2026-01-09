# GENIE Profiling Guide

GENIE includes an internal lightweight profiler to help analyze performance bottlenecks and understand where computation time is spent during execution.

## Building with Profiling Support

To enable profiling instrumentation, build GENIE with the `GENIE_PROFILE` CMake option:

```bash
mkdir build && cd build
cmake -DGENIE_PROFILE=ON ..
make -j
```

**Note**: If you build without `-DGENIE_PROFILE=ON`, the profiling code will compile to no-ops with zero runtime overhead.

## Running GENIE with Profiling

Once built with profiling support, enable profiling output using the `-prof` or `--profile` flag:

```bash
./GENIE --genotype example/test \
        --phenotype example/pheno.txt \
        --environment example/env.txt \
        --model G+GxE+NxE \
        --profile \
        --profile_out genie_profile.csv
```

### Profiling Options

- `-prof`, `--profile`: Enable profiling and output timing data (requires `GENIE_PROFILE` build flag)
- `--profile_out <path>`: Path to profiling output file (default: `genie_profile.csv`)

The output format is automatically determined by the file extension:
- `.csv` → CSV format (default)
- `.json` → JSON format

## Understanding the Profiling Output

### CSV Format

The CSV output contains the following columns:

```csv
name,calls,total_seconds,avg_seconds,total_bytes
matvec_Xv,1500,45.23,0.030153,0
jackknife_pass1,1,120.45,120.45,0
solve_normal_equations,100,2.34,0.0234,0
io_read_genotype,3,0.89,0.296667,0
```

- `name`: Name of the profiled region
- `calls`: Number of times the region was called
- `total_seconds`: Total wall-clock time spent in this region (sum across all calls)
- `avg_seconds`: Average time per call (`total_seconds / calls`)
- `total_bytes`: Total bytes transferred (currently used for I/O tracking)

### JSON Format

The JSON output provides the same information in a structured format:

```json
[
  {
    "name": "matvec_Xv",
    "calls": 1500,
    "total_seconds": 45.234567
  },
  {
    "name": "jackknife_pass1",
    "calls": 1,
    "total_seconds": 120.456789,
    "total_bytes": 0
  }
]
```

## Profiled Regions

GENIE instruments the following key computational regions:

### Matrix-Vector Operations
- `matvec_Xv`: Genotype matrix × vector multiplication (Y = X * v)
- `matvec_Xt_v`: Transposed genotype matrix × vector (Z = X^T * v)

### Main Computation Loops
- `jackknife_pass1`: First pass over jackknife blocks
- `jackknife_pass2`: Second pass over jackknife blocks (if applicable)
- `jackknife_block`: Processing of individual jackknife blocks

### Linear Algebra
- `solve_normal_equations`: Solving the normal equations for variance components

### I/O Operations
- `io_read_genotype`: Reading PLINK BED/BIM/FAM files
- `io_read_phenotype`: Reading phenotype file
- `io_read_covariate`: Reading covariate file
- `io_read_environment`: Reading environment file
- `io_read_annotation`: Reading annotation file

## Interpreting Results

The profiler output is sorted by total time (descending), making it easy to identify performance bottlenecks:

1. **Identify hotspots**: Regions with the highest `total_seconds` are where most time is spent
2. **Check call counts**: High `calls` with low `avg_seconds` might indicate loop overhead
3. **Validate expected behavior**:
   - Matrix operations should dominate for large datasets
   - I/O should be fast relative to computation
   - Jackknife loops should scale with the number of jackknife blocks

## Example Analysis

```csv
name,calls,total_seconds,avg_seconds,total_bytes
matvec_Xv,2400,180.5,0.0752,0
jackknife_pass1,1,195.2,195.2,0
matvec_Xt_v,2400,165.3,0.0689,0
solve_normal_equations,100,5.6,0.056,0
io_read_genotype,3,2.1,0.7,0
```

**Interpretation**:
- Matrix-vector operations consume ~345s total (180.5 + 165.3)
- This represents ~89% of the jackknife pass1 time (345/195.2)
- I/O is only 1% of total time (2.1/195.2), indicating I/O is not a bottleneck
- Linear solver is fast (~5.6s total), suggesting convergence is good

## Profiling Overhead

When compiled with `-DGENIE_PROFILE=ON`, the profiler adds minimal overhead:
- Thread-local start/stop timing (~10-20 nanoseconds per call)
- Atomic aggregation when timers stop
- No heap allocations in hot paths

For production runs where profiling is not needed, simply build without `-DGENIE_PROFILE=ON` and all profiling code compiles to zero-overhead no-ops.

## Troubleshooting

**Q: I get "profiling not compiled in" warning**
A: Rebuild GENIE with `cmake -DGENIE_PROFILE=ON ..`

**Q: Profiling output file is empty**
A: Check that the program completed successfully and didn't exit early due to errors

**Q: Some regions show 0 calls**
A: Those code paths may not have been executed for your particular input/model configuration

## Advanced Usage

### Combining with External Profilers

The `-fno-omit-frame-pointer` flag is automatically added when profiling is enabled, making it easier to use external profilers like `perf` or `Instruments.app`:

```bash
# Build with profiling
cmake -DGENIE_PROFILE=ON ..
make -j

# Run with both internal and external profiling
perf record -g ./GENIE <args> --profile

# Analyze
perf report
```

### Custom Profile Regions

To add custom profiling to your own code modifications:

```cpp
#include "profiler.h"

void my_function() {
    ScopedTimer timer("my_custom_region");
    // Your code here
    // Timer automatically stops when it goes out of scope
}
```

The profiler is thread-safe and can be used in multithreaded code.
