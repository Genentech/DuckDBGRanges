# DuckDBGRanges 0.99.1

## Documentation

- Added `\value` sections (roxygen `@return`) to the `DuckDBGRanges` and
  `DuckDBGRangesList` class and utility man pages, documenting the values
  returned by the constructors, accessors, replacement, coercion, and
  range-operation methods. Resolves the `R CMD BiocCheck` "missing \value"
  WARNING.

# DuckDBGRanges 0.9.8

## Testing and diagnostics

- Added a query-plan **regression guard** for interval overlaps (no change to
  overlap results or their current speed). A new internal `.explainQuery()`
  returns DuckDB's plan for a lazy query, and a test uses it to assert the
  overlap join stays an IEJoin / range join rather than degrading to a
  `NESTED_LOOP_JOIN` over the full cross product (which OOMs on skewed inputs).
  Interval overlap is a range join, not an ASOF (nearest-match) join. The
  overlap join was factored into an internal `.overlap_join_tbl()` helper so its
  plan can be inspected without collecting.

# DuckDBGRanges 0.9.7

## Changes

- Relicensed under the MIT License.

# DuckDBGRanges 0.9.6

## Bug fixes

- `precede()`/`follow()` with `select = "all"` now return only the subject(s) at
  the **nearest** distance (the ties), matching base GenomicRanges — previously
  they returned every directional subject. The `select = "first"`/`"last"` paths
  were already correct; the `"all"` branch was missing the min-distance filter.
- `distance()` now returns `NA` for a paired range on a **different seqname** or
  (unless `ignore.strand = TRUE`) an **incompatible strand** (`+` vs `-`; `*`
  matches any strand), matching base `GenomicRanges::distance`. Previously it
  ignored seqnames and strand and returned a bare coordinate gap. Added an
  `ignore.strand` argument.

# DuckDBGRanges 0.9.5

## Bug fixes

- `narrow()` now treats a **negative `start`** as a position counting back from
  the range end (`-1` is the last base), matching base `IRanges::narrow`
  (`solveUserSEW`) and the existing negative-`end` behaviour. Previously a
  negative `start` was applied as `start + (start - 1)`, walking left of the
  range start and producing an invalid interval. The start/end/width expressions
  are rebuilt so every sign-and-`width` combination is correct; because the
  `datacols` expressions flatten to SQL without parentheses, each is kept as a
  `column +/- scalar` (never a subtraction of a compound). Oracle-tested against
  base `narrow` across positive/negative `start`/`end` and their `width`
  combinations.
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
