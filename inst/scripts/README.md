# DuckDBGRanges benchmark scripts

These scripts produce the numbers rendered in the *Benchmarking DuckDBGRanges*
vignette. They follow the
[HDF5Array performance vignette](https://bioconductor.org/packages/release/bioc/vignettes/HDF5Array/inst/doc/HDF5Array_performance.html)
precompute pattern: the full-scale benchmark runs **offline** and writes a results
file that the vignette renders, so the vignette itself builds quickly.

## Files

- **`run_vignette_benchmarks.R`** — generates synthetic ranges at scale (a
  scATAC-seq peak scenario and a variant scenario), times six operations on both
  in-memory `GRanges` and on-disk `DuckDBGRanges`, records each object's memory
  footprint, and writes `benchmark_results.rds`.
- **`make_timings_table.R`** — helpers the vignette `source()`s to render the
  results (`load_vignette_timings`, `make_timings_table`, `make_memory_table`,
  `timings_config_note`). If `benchmark_results.rds` is absent, the vignette shows
  a short note instead of tables.

## Regenerating the results

```sh
# full run (1M peaks, 10M variants) — do this on a machine with ample RAM
Rscript inst/scripts/run_vignette_benchmarks.R
# then copy benchmark_results.rds into inst/scripts/ and commit it
```

Environment variables (all optional):

| Variable | Default | Meaning |
|----------|---------|---------|
| `BENCH_CORES` | available − 1 | core budget for DuckDB |
| `BENCH_SCATAC_PEAKS` | 1000000 | peak count for the scATAC-seq scenario |
| `BENCH_VARIANTS` | 10000000 | variant count for the variant scenario |
| `BENCH_ANNOTATION` | 50000 | annotation features used as the overlap subject |
| `BENCH_SCENARIOS` | `scATAC,variant` | subset of scenarios to run |

To smoke-test the harness quickly before a full run, set small counts:

```sh
BENCH_SCATAC_PEAKS=5000 BENCH_VARIANTS=5000 Rscript inst/scripts/run_vignette_benchmarks.R
```

`GRanges` runs single-threaded C code; `DuckDBGRanges` pushes filters into the
Parquet reader and autotunes DuckDB's internal threads up to the core budget. The
per-run configuration is recorded in `attr(results, "config")` and shown beneath
the vignette tables.
