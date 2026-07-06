#!/usr/bin/env Rscript
# Best-effort comparison for the "Benchmarking DuckDBGRanges" vignette.
#
# Compares in-memory GRanges against on-disk DuckDBGRanges on two synthetic
# genomic-range scenarios at scale:
#
#   * scATAC : peak-like intervals (a few hundred bp wide) clustered across the
#              autosomes -- the shape of multi-sample scATAC-seq consensus peaks.
#   * variant: SNV-like point positions across the autosomes -- the shape of a
#              large variant catalogue.
#
# For each scenario it times the range operations that matter for these workflows
# on BOTH backends, and records the in-memory object footprint of each:
#
#   restrict, reduce, disjoin           -- intra/inter-range transforms
#   findOverlaps, subsetByOverlaps      -- overlap queries against an annotation set
#   distanceToNearest                   -- nearest-feature distance
#
# GRanges is single-threaded C code; DuckDBGRanges pushes filters into the Parquet
# reader and autotunes DuckDB's internal threads up to the core budget. Each timing
# is wrapped so an unsupported/failed op records NA rather than aborting the run.
#
# Writes benchmark_results.rds: a tidy data.frame (Scenario, Operation, Backend,
# Seconds) with attr(., "memory") (object footprint) and attr(., "config") (run
# parameters). The vignette renders it via inst/scripts/make_timings_table.R.
#
# Usage:  Rscript run_vignette_benchmarks.R
# Env:    BENCH_CORES         core budget for DuckDB (default = available - 1)
#         BENCH_SCATAC_PEAKS  scATAC peak count      (default 1000000)
#         BENCH_VARIANTS      variant count          (default 10000000)
#         BENCH_ANNOTATION    annotation feature count for overlaps (default 50000)
#         BENCH_SCENARIOS     comma-separated subset  (default "scATAC,variant")
#         # Tip: set small counts (e.g. BENCH_SCATAC_PEAKS=5000 BENCH_VARIANTS=5000)
#         # to smoke-test the harness quickly before the full run.

suppressPackageStartupMessages({
    library(DuckDBGRanges)
    library(GenomicRanges)
    library(DuckDBDataFrame)
    library(arrow)
    library(DBI)
})

## ---- configuration -----------------------------------------------------------
# Cores AVAILABLE to this process. detectCores() reports the whole machine even
# when a scheduler grants fewer; prefer the SLURM allocation, then nproc, then
# detectCores() as a last resort.
available_cores <- function() {
    for (v in c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE")) {
        x <- suppressWarnings(as.integer(Sys.getenv(v, "")))
        if (length(x) == 1L && !is.na(x) && x > 0L) return(x)
    }
    np <- tryCatch(as.integer(system("nproc", intern = TRUE, ignore.stderr = TRUE)),
                   error = function(e) NA_integer_, warning = function(w) NA_integer_)
    if (length(np) == 1L && !is.na(np) && np > 0L) return(np)
    max(1L, parallel::detectCores())
}
cores <- as.integer(Sys.getenv("BENCH_CORES",
                                as.character(max(1L, available_cores() - 1L))))
scatac_peaks <- as.numeric(Sys.getenv("BENCH_SCATAC_PEAKS", "1000000"))
variants     <- as.numeric(Sys.getenv("BENCH_VARIANTS", "10000000"))
n_annot      <- as.numeric(Sys.getenv("BENCH_ANNOTATION", "50000"))
scenarios <- trimws(strsplit(Sys.getenv("BENCH_SCENARIOS", "scATAC,variant"), ",")[[1]])
scenarios <- scenarios[nzchar(scenarios)]

# hg38 autosome lengths -- positions are drawn proportionally to length.
chrom_len <- c(
    chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
    chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
    chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
    chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
    chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
    chr21 = 46709983, chr22 = 50818468)

cat(sprintf("Config: scATAC = %s peaks | variant = %s variants | annotation = %s | cores = %d\n",
            format(scatac_peaks, big.mark = ","), format(variants, big.mark = ","),
            format(n_annot, big.mark = ","), cores))

## ---- helpers -----------------------------------------------------------------
elapsed <- function(expr)
    tryCatch(unname(system.time(force(expr))["elapsed"]),
             error = function(e) { message("    (failed: ", conditionMessage(e), ")"); NA_real_ })

mb <- function(x) tryCatch(round(as.numeric(object.size(x)) / 2^20, 4),
                           error = function(e) NA_real_)

# Draw `n` intervals of the given width distribution across the autosomes.
make_ranges_df <- function(n, width_fun) {
    n <- as.integer(n)
    p <- chrom_len / sum(chrom_len)
    seqn <- sample(names(chrom_len), n, replace = TRUE, prob = p)
    w <- width_fun(n)
    ## vectorized start draw: uniform in [1, chrom_len - width] per row
    maxstart <- pmax(1L, chrom_len[seqn] - w)
    start <- as.integer(floor(runif(n) * maxstart)) + 1L
    data.frame(seqnames = seqn, start = start, end = start + w - 1L,
               strand = sample(c("+", "-", "*"), n, replace = TRUE),
               score = runif(n), stringsAsFactors = FALSE)
}

# Build both backends (in-memory GRanges + on-disk DuckDBGRanges) from a df.
# DuckDBGRanges is constructed WITHOUT an explicit keycol so it keeps its small,
# constant footprint (a row_number identifier is generated internally).
make_backends <- function(df) {
    gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
    path <- tempfile(fileext = ".parquet")
    write_parquet(df, path)
    ddb <- DuckDBGRanges(path, seqnames = "seqnames", start = "start",
                         end = "end", strand = "strand", mcols = "score")
    list(gr = gr, ddb = ddb)
}

# The operation set, as per-backend thunks. `subj` is the annotation set.
op_thunks <- function(x, subj) {
    list(
        restrict         = function() restrict(x, start = 1L, end = 50000000L),
        reduce           = function() reduce(x),
        disjoin          = function() disjoin(x),
        findOverlaps     = function() findOverlaps(x, subj),
        subsetByOverlaps = function() subsetByOverlaps(x, subj),
        distanceToNearest = function() distanceToNearest(x, subj))
}

## ---- run ---------------------------------------------------------------------
# Cap DuckDB's internal threads at the core budget on the shared connection.
conn <- DuckDBDataFrame::acquireDuckDBConn()
try(dbExecute(conn, sprintf("SET threads = %d;", cores)), silent = TRUE)

rows <- list(); mem_rows <- list()
for (scn in scenarios) {
    cat(sprintf("\n=== scenario: %s ===\n", scn))
    if (scn == "scATAC") {
        df <- make_ranges_df(scatac_peaks, function(n) sample(200:1000, n, replace = TRUE))
    } else if (scn == "variant") {
        df <- make_ranges_df(variants, function(n) rep(1L, n))  # SNVs
    } else {
        message("  unknown scenario '", scn, "', skipping"); next
    }
    annot_df <- make_ranges_df(n_annot, function(n) sample(500:5000, n, replace = TRUE))

    b <- make_backends(df)
    a <- make_backends(annot_df)
    mem_rows[[length(mem_rows) + 1L]] <-
        data.frame(Scenario = scn, Backend = "GRanges", MB = mb(b$gr))
    mem_rows[[length(mem_rows) + 1L]] <-
        data.frame(Scenario = scn, Backend = "DuckDBGRanges", MB = mb(b$ddb))
    cat(sprintf("  footprint: GRanges %.1f MB | DuckDBGRanges %.4f MB\n",
                mb(b$gr), mb(b$ddb)))

    gr_ops  <- op_thunks(b$gr,  a$gr)
    ddb_ops <- op_thunks(b$ddb, a$ddb)
    for (op in names(gr_ops)) {
        cat(sprintf("  --- %s ---\n", op))
        s_gr  <- elapsed(gr_ops[[op]]())
        s_ddb <- elapsed(ddb_ops[[op]]())
        rows[[length(rows) + 1L]] <-
            data.frame(Scenario = scn, Operation = op, Backend = "GRanges", Seconds = s_gr)
        rows[[length(rows) + 1L]] <-
            data.frame(Scenario = scn, Operation = op, Backend = "DuckDBGRanges", Seconds = s_ddb)
        cat(sprintf("    GRanges %s | DuckDBGRanges %s\n",
                    if (is.na(s_gr)) "NA" else sprintf("%.3f s", s_gr),
                    if (is.na(s_ddb)) "NA" else sprintf("%.3f s", s_ddb)))
    }
    rm(b, a, df, annot_df); gc(FALSE)
}

results <- do.call(rbind, rows)
attr(results, "memory") <- do.call(rbind, mem_rows)
attr(results, "config") <- list(scatac_peaks = scatac_peaks, variants = variants,
                                n_annotation = n_annot, cores = cores)

saveRDS(results, "benchmark_results.rds")
cat("\nSaved benchmark_results.rds\n")
print(results, row.names = FALSE)
