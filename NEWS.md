# DuckDBGRanges 0.9.5

## Bug fixes

- `nearest(x)` and `distanceToNearest(x)` (no `subject`) now exclude self-hits,
  matching base GenomicRanges' `drop.self = TRUE` for the missing-subject form:
  a range is never its own nearest neighbour, so its next-nearest is returned
  instead of itself (distance 0). Previously both delegated to the two-argument
  form against `x` itself, so every range matched itself. The explicit
  `nearest(x, x)` / `distanceToNearest(x, x)` forms still keep self-hits, as in
  base. (`precede`/`follow` were already correct — their strict inequalities
  never admit a self-match.) Added a self-query test comparing to the base oracle.
- `nearest`/`precede`/`follow`/`distanceToNearest` now treat strand `"*"` as
  compatible with any strand, matching base GenomicRanges: `"+"` pairs with
  `{+,*}`, `"-"` with `{-,*}`, and `"*"` with all. The neighbour join used a
  strict equi-join on strand, which silently dropped every `"*"` pair (a `"*"`
  query against `"+"`/`"-"` subjects returned `NA`). The join is now on
  `seqnames` only, followed by a strand-compatibility filter.
- `precede()` and `follow()` are now strand-directional. Base defines them in the
  transcription direction — for a `"-"` strand query the roles of upstream and
  downstream are reversed, and for a `"*"` query the direction is chosen per
  subject strand. The previous implementation used a fixed genomic-coordinate
  direction regardless of strand, so results were inverted on the `"-"` strand.
  The convention is now selected per row (`subj_start > end` vs `subj_end <
  start`) via `use_minus = strand == "-" | (strand == "*" & subj_strand == "-")`,
  ranking by the transcription-direction gap; `ignore.strand = TRUE` collapses to
  the `"+"` convention as before. `follow()` keeps its base `select = "last"`
  (largest index) tie-break, `precede()` its `select = "first"` (smallest index).
  Added oracle tests against base GenomicRanges across `+`/`-`/`*` queries with
  mixed-strand subjects, plus tie-break tests.

# DuckDBGRanges 0.9.4

## Bug fixes

- `resize(x, width, fix = "center")` now matches base `IRanges::resize()`:
  `new_start = start + (width(x) - width) %/% 2`, anchored at the (strand-independent)
  center with an exact integer width. The previous implementation used the midpoint
  `(start + end) / 2` with true division --- yielding fractional coordinates for
  even-width ranges --- and `center + width %/% 2 - 1`, which returned `width - 1`
  for odd target widths. Added `fix = "center"` cases to the `resize` tests.

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
