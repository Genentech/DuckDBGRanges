# DuckDBGRanges

*A DuckDB-backed GenomicRanges — out-of-core genomic intervals, stored as Parquet.*

## Overview

`DuckDBGRanges` provides a [GenomicRanges](https://bioconductor.org/packages/GenomicRanges)
`GRanges` (and `GRangesList`) backed by [DuckDB](https://duckdb.org) over columnar
**Parquet**. It exposes the full `GenomicRanges` API --- `seqnames()`, `start()`,
`restrict()`, `findOverlaps()` --- but the coordinates stay **on disk** and
operations are recorded as lazy SQL queries, so you can filter by region,
chromosome, or metadata **without loading the ranges into memory**. The R object
is a small, constant ~175 KB whether it fronts a thousand ranges or a hundred
million.

It builds on the `DuckDBDataFrame` foundation of the **BiocDuckDB** suite and is
designed for the *filter-first* shape of most genomic work: use DuckDB to cut a
large interval set down to the ranges of interest, then materialize that subset to
an ordinary `GRanges` for the interval-tree operations R does best.

## Installation

```r
# once available from Bioconductor:
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("DuckDBGRanges")
```

## Quick start

```r
library(DuckDBGRanges)
library(GenomicRanges)
library(arrow)

gr <- GRanges(c("chr1:100-150:+", "chr1:200-275:-", "chr2:300-360:+"),
              score = c(10, 20, 15))
gr_df <- as.data.frame(gr)
path <- tempfile(fileext = ".parquet")
write_parquet(gr_df, path)

gr_ddb <- DuckDBGRanges(path, seqnames = "seqnames", start = "start",
                        end = "end", strand = "strand", mcols = "score")

gr_ddb[seqnames(gr_ddb) == "chr1" & gr_ddb$score > 12]   # lazy, on-disk filter
```

Standard `GRanges` accessors and range operations (`restrict()`, `reduce()`,
`shift()`, `coverage()`, ...) work as usual; `as(x, "GRanges")` materializes a
(typically filtered) subset into memory.

## Performance

Full numbers are in the *Benchmarking DuckDBGRanges* vignette; the shape of the
result, comparing in-memory `GRanges` with `DuckDBGRanges`:

| | GRanges | DuckDBGRanges |
|---|---|---|
| Object footprint | proportional to N | **constant (~175 KB; 150–1500x smaller)** |
| Region / metadata filtering | scans all ranges | **predicate pushdown (2–13x faster)** |
| `distanceToNearest` | interval trees | **SQL aggregation (10–70x faster)** |
| `reduce` / `disjoin` / `coverage` | optimized C | SQL (comparable; faster at scale) |
| `findOverlaps` / `subsetByOverlaps` | **interval trees (fast)** | SQL join (slower) |

The idiomatic **hybrid workflow** plays to both strengths:

```r
# 1. filter with DuckDB (memory-efficient, predicate pushdown)
sub <- restrict(gr_ddb, start = 1, end = 50e6)
sub <- sub[sub$score > 12]

# 2. materialize the small result, then use GRanges interval trees for overlaps
findOverlaps(as(sub, "GRanges"), annotation)
```

## When to use DuckDBGRanges

A good fit when the interval set is **larger than memory**, when the workload is
**filter-heavy** (restricting by region/chromosome/metadata before the expensive
step), when you want to keep a huge annotation *referenceable* at constant memory
cost, or when the data already lives on disk as **Parquet** other tools can read.
An in-memory `GRanges` remains preferable when the data fits in RAM and for
overlap-heavy work.

## Documentation

- **Introduction to DuckDBGRanges** — motivation, construction, and the common
  operations (`vignettes/DuckDBGRanges.Rmd`).
- **Benchmarking DuckDBGRanges** — a best-effort comparison against in-memory
  `GRanges` on scATAC-seq and variant-scale data (`vignettes/DuckDBGRanges-benchmark.Rmd`).
- **Design and extension of DuckDBGRanges** — class structure, SQL translation, and
  the `LIST[]` representation of grouped ranges, for developers
  (`vignettes/DuckDBGRanges-design.Rmd`).

## License

MIT License. Copyright Genentech, Inc., 2026.
