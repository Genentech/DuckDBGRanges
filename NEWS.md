# DuckDBGRanges 0.9.3

## Documentation

- Restructured the vignettes into a user-first set, replacing the single
  internals-heavy *Architecture of the DuckDBGRanges Package*:
  - *Introduction to DuckDBGRanges* --- motivation, construction, and the common
    operations (accessors, subsetting, range operations, the filter-then-materialize
    workflow, `DuckDBGRangesList`).
  - *Benchmarking DuckDBGRanges* --- a best-effort comparison against in-memory
    `GRanges` on scATAC-seq (1M peaks) and variant (10M variants) scenarios,
    rendered from precomputed results so the vignette builds quickly.
  - *Design and extension of DuckDBGRanges* --- the five coordinate columns, class
    structure, SQL translation, and the `LIST[]` representation of grouped ranges,
    for developers.
- Added `inst/scripts/` with the offline benchmark generator
  (`run_vignette_benchmarks.R`) and the vignette table helpers
  (`make_timings_table.R`), following the `HDF5Array` performance-vignette
  precompute pattern.
- Rewrote the README.
