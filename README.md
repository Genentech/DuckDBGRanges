# DuckDBGRanges

## Memory-Efficient DuckDB Backend for GenomicRanges

DuckDBGRanges provides DuckDB-backed implementations of `GRanges` and `GRangesList` for efficient, out-of-memory operations on large genomic interval datasets. With a constant ~175 KB memory footprint regardless of data size, it enables analysis of millions of genomic features on memory-constrained systems.

### Why DuckDB for Genomic Ranges?

Genomic workflows often involve filtering millions of variants, peaks, or features by region, chromosome, or metadata before performing overlap operations. DuckDBGRanges excels at this filtering phase:

- **Constant memory**: ~175 KB footprint for 1M or 100M ranges
- **Fast filtering**: SQL predicate pushdown makes region/metadata filters 5-17x faster
- **Lazy evaluation**: Build complex filter pipelines without loading data
- **Hybrid workflows**: Filter with DuckDB, then materialize to GRanges for overlaps

## Performance Highlights

### Variant Analysis (10M variants)

| Operation | GRanges | DuckDBGRanges | Winner |
|-----------|---------|---------------|--------|
| Memory footprint | 76 MB | 174 KB | **DuckDB (450x smaller)** |
| `restrict()` (50Mb region) | 2.7 sec | 0.16 sec | **DuckDB (17.5x faster)** |
| Filter pipeline | 2.1 sec | 0.60 sec | **DuckDB (3.5x faster)** |
| `findOverlaps()` | 0.018 sec | 6.4 sec | GRanges (interval trees) |
| `subsetByOverlaps()` | 0.016 sec | 7.0 sec | GRanges (interval trees) |

### scATAC-seq Peak Analysis (1M peaks)

| Operation | GRanges | DuckDBGRanges | Winner |
|-----------|---------|---------------|--------|
| Memory footprint | 30 MB | 175 KB | **DuckDB (176x smaller)** |
| `restrict()` (filter) | 0.42 sec | 0.08 sec | **DuckDB (5.3x faster)** |
| `distanceToNearest()` | 0.49 sec | 0.24 sec | **DuckDB (2x faster)** |
| `reduce()` | 0.19 sec | 0.70 sec | GRanges |
| `findOverlaps()` | 0.09 sec | 34.3 sec | GRanges (interval trees) |

## Core Classes

### DuckDBGRanges

Extends GenomicRanges `GRanges`:

```r
library(DuckDBGRanges)
library(GenomicRanges)

# Create genomic ranges data
gr_df <- data.frame(
    seqnames = rep(paste0("chr", 1:22), each = 50000),
    start = sample(1:1e8, 1.1e6),
    end = start + sample(1:1000, 1.1e6, replace = TRUE),
    strand = sample(c("+", "-", "*"), 1.1e6, replace = TRUE),
    score = runif(1.1e6),
    gene_id = paste0("GENE", 1:1.1e6)
)

# Write to Parquet
library(arrow)
write_parquet(gr_df, "ranges.parquet")

# Create DuckDBGRanges (only ~175 KB in memory!)
gr_ddb <- DuckDBGRanges(
    "ranges.parquet",
    seqnames = "seqnames",
    start = "start",
    end = "end",
    strand = "strand",
    mcols = c("score", "gene_id")
)

# Standard GRanges operations
length(gr_ddb)
seqnames(gr_ddb)
start(gr_ddb)
mcols(gr_ddb)
```

### DuckDBGRangesList

Extends GenomicRanges `GRangesList`:

```r
# Grouped ranges (e.g., exons by transcript)
grl_ddb <- DuckDBGRangesList(
    "transcripts.parquet",
    seqnames = "seqnames",
    start = "start",
    end = "end",
    strand = "strand",
    group = "transcript_id"
)

# Access by group
grl_ddb[[1]]
lengths(grl_ddb)
```

## Key Features

### Memory-Efficient Storage

```r
# In-memory GRanges: 76 MB
gr <- GRanges(
    seqnames = rep(paste0("chr", 1:22), each = 50000),
    ranges = IRanges(start = sample(1:1e8, 1.1e6),
                     width = sample(1:1000, 1.1e6)),
    score = runif(1.1e6)
)
object.size(gr)  # 79,621,984 bytes

# DuckDBGRanges: 175 KB metadata + disk storage
object.size(gr_ddb)  # 179,176 bytes (450x smaller!)
```

### Fast Region Filtering

```r
# Restrict to chromosome and region (17.5x faster than GRanges)
chr1_region <- restrict(gr_ddb, 
                        seqnames = "chr1",
                        start = 1e6, 
                        end = 10e6)

# Filter by metadata
high_score <- gr_ddb[gr_ddb$score > 0.8]

# Complex SQL queries
filtered <- gr_ddb[gr_ddb$score > 0.5 & gr_ddb$gene_id %in% genes_of_interest]
```

### Hybrid Workflow Pattern

**Recommended**: Filter with DuckDB, then materialize to GRanges for overlaps:

```r
# Step 1: Filter with DuckDB (fast, memory-efficient)
chr1_variants <- restrict(gr_ddb, seqnames = "chr1", start = 1e6, end = 50e6)
high_qual <- chr1_variants[chr1_variants$qual > 30]

# Step 2: Materialize to GRanges (only filtered subset!)
high_qual_gr <- as(high_qual, "GRanges")

# Step 3: Use GRanges for overlaps (fast interval trees)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
overlaps <- findOverlaps(high_qual_gr, genes)
```

This hybrid approach gives you:
- ✅ Memory efficiency during filtering (DuckDB)
- ✅ Fast overlap operations (GRanges interval trees)
- ✅ Best of both worlds!

## GenomicRanges Operations

### When to Use DuckDB vs GRanges

**DuckDBGRanges excels at:**
- `restrict()` - Region filtering
- `seqnames()` filtering - Chromosome selection  
- `mcols()` filtering - Metadata queries
- `distanceToNearest()` - Distance calculations
- Pipeline operations with lazy evaluation

**GRanges excels at:**
- `findOverlaps()` - Overlap detection
- `subsetByOverlaps()` - Overlap-based subsetting
- `reduce()` / `disjoin()` - Range reduction
- `union()` / `intersect()` / `setdiff()` - Set operations
- `coverage()` - Coverage calculation

### SQL-Optimized Operations

| R Function | SQL Translation | Performance Gain |
|------------|-----------------|------------------|
| `restrict()` | `WHERE start >= x AND end <= y` | **17.5x faster** |
| `seqnames()` filter | `WHERE seqnames = 'chr1'` | Lazy evaluation |
| `distanceToNearest()` | `MIN(ABS(pos - subject_pos))` | **2x faster** |
| `mcols()` filter | `WHERE score > x` | Predicate pushdown |

## Quick Start

### Variant Analysis Workflow

```r
library(DuckDBGRanges)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load variants as DuckDBGRanges
variants <- DuckDBGRanges(
    "variants.parquet",
    seqnames = "chr",
    start = "pos",
    end = "pos",
    strand = "*",
    mcols = c("ref", "alt", "qual", "filter")
)

# Filter by quality (fast SQL)
high_qual <- variants[variants$qual > 30 & variants$filter == "PASS"]

# Restrict to gene region (fast SQL)
brca1_region <- restrict(high_qual, 
                         seqnames = "chr17",
                         start = 43044295,
                         end = 43125483)

# Materialize for overlap annotation
brca1_gr <- as(brca1_region, "GRanges")

# Annotate with GRanges (fast interval trees)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
exons <- exons(txdb)
coding_variants <- subsetByOverlaps(brca1_gr, exons)
```

### scATAC-seq Peak Analysis

```r
# Load 1M peaks as DuckDBGRanges
peaks <- DuckDBGRanges(
    "peaks.parquet",
    seqnames = "seqnames",
    start = "start",
    end = "end",
    strand = "*",
    mcols = c("score", "peak_id", "sample")
)

# Filter by chromosome and score (5.3x faster)
chr1_peaks <- restrict(peaks, seqnames = "chr1")
high_peaks <- chr1_peaks[chr1_peaks$score > 10]

# Find nearest TSS (2x faster than GRanges)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tss <- promoters(genes(txdb), upstream = 0, downstream = 1)
nearest_tss <- distanceToNearest(high_peaks, as(tss, "DuckDBGRanges"))

# For complex overlaps, materialize first
high_peaks_gr <- as(high_peaks, "GRanges")
promoter_peaks <- subsetByOverlaps(high_peaks_gr, promoters(genes(txdb)))
```

## Data Persistence

### Writing GRanges to Parquet

```r
# Simple pattern - convert to data.frame and write
gr <- GRanges("chr1", IRanges(1:1000, width = 100), score = rnorm(1000))
gr_df <- as.data.frame(gr)
arrow::write_parquet(gr_df, "ranges.parquet")

# Load back
gr_ddb <- DuckDBGRanges("ranges.parquet",
                        seqnames = "seqnames",
                        start = "start",
                        end = "end",
                        strand = "strand",
                        mcols = "score")
```

For complex objects with nested metadata, use `BiocDuckDB::writeParquet()` which handles schema and metadata automatically.

## When to Use DuckDBGRanges

**Recommended for:**
- Variant analysis with millions of variants
- scATAC-seq peak calling and analysis
- Memory-constrained environments
- Cloud-based workflows (Parquet on S3/GCS)
- Initial filtering before overlap operations
- Cross-language data sharing (Parquet is universal)

**Consider GRanges when:**
- Data fits in memory
- Primary operations are overlap-based
- You need interval tree performance
- Complex inter-range operations (`reduce`, `disjoin`)

**Hybrid workflow** (recommended):
1. Filter/query with DuckDBGRanges (memory-efficient, fast SQL)
2. Materialize to GRanges for overlaps (fast interval trees)
3. Best of both worlds!

## Documentation

- **[DuckDBGRanges Classes](vignettes/DuckDBGRanges-classes.Rmd)**: Architecture and design
- **[Variant Analysis](vignettes/DuckDBGRanges-variant-analysis.Rmd)**: Comprehensive variant filtering and annotation workflow
- **[scATAC-seq Analysis](vignettes/DuckDBGRanges-scATAC-seq.Rmd)**: Peak analysis and chromatin accessibility workflows

## Installation

```r
# Requires DuckDBDataFrame
# install.packages("remotes")
remotes::install_github("your-org/DuckDBDataFrame")
remotes::install_github("your-org/DuckDBGRanges")
```

## Dependencies

DuckDBGRanges depends on:
- **DuckDBDataFrame**: Foundation for DuckDB-backed structures
- **Bioconductor**: GenomicRanges, IRanges, S4Vectors, Seqinfo, BiocGenerics
- **Database**: DBI, dplyr, dbplyr

## Contributing

Contributions are welcome! Please:
- Report issues through GitHub
- Include benchmarks for performance claims
- Test with real genomic datasets
- Follow Bioconductor standards

## License

DuckDBGRanges is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgements

Special thanks to:
- The GenomicRanges team for the excellent API
- The DuckDB team for query optimization
- The Bioconductor project for infrastructure
- The Apache Arrow project for Parquet format
