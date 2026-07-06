## make_timings_table.R
##
## Helpers used by the "Benchmarking DuckDBGRanges" vignette to render PRECOMPUTED
## benchmark results, so the vignette does not run the (multi-minute, large-memory)
## full-scale benchmark at build time. Mirrors the HDF5Array performance vignette's
## precompute approach.
##
## The results are produced offline by run_vignette_benchmarks.R and saved as
## inst/scripts/benchmark_results.rds. To regenerate (see that script's header):
##
##   Rscript inst/scripts/run_vignette_benchmarks.R
##   # then copy benchmark_results.rds into inst/scripts/ and commit it
##
## benchmark_results.rds is a tidy data.frame with columns
##   Scenario (scATAC/variant), Operation, Backend (GRanges/DuckDBGRanges), Seconds
## and two attributes:
##   attr(., "memory"): data.frame(Scenario, Backend, MB) -- object footprint
##   attr(., "config"): list(scatac_peaks, variants, cores, ...) describing the run

## Load the shipped results, or NULL if absent (so the vignette degrades cleanly).
load_vignette_timings <- function() {
    path <- system.file("scripts", "benchmark_results.rds", package = "DuckDBGRanges")
    if (!nzchar(path) || !file.exists(path)) {
        return(NULL)
    }
    readRDS(path)
}

.BACKENDS <- c("GRanges", "DuckDBGRanges")

## Which backend wins each operation (lower seconds), as a readable label.
.winner <- function(gr, ddb) {
    if (is.na(gr) || is.na(ddb)) return("")
    if (ddb < gr) "DuckDBGRanges" else "GRanges"
}

## Pivot one scenario's timings to wide form (one row per operation) with a
## speedup column (GRanges / DuckDBGRanges; > 1 means DuckDB is faster).
timings_wide <- function(results, scenario) {
    df <- results[results$Scenario == scenario, , drop = FALSE]
    ops <- unique(df$Operation)
    out <- data.frame(Operation = ops, stringsAsFactors = FALSE)
    for (b in .BACKENDS) {
        out[[b]] <- vapply(ops, function(o) {
            v <- df$Seconds[df$Operation == o & df$Backend == b]
            if (length(v) == 0L) NA_real_ else v[[1L]]
        }, numeric(1))
    }
    out[["Speedup"]] <- round(out[["GRanges"]] / out[["DuckDBGRanges"]], 1)
    out[["Winner"]] <- mapply(.winner, out[["GRanges"]], out[["DuckDBGRanges"]])
    out
}

## Render one scenario's timings as a knitr::kable, or a short note if the
## precomputed results are absent from this build.
make_timings_table <- function(scenario, caption = NULL,
                               results = load_vignette_timings()) {
    if (is.null(results)) {
        cat("_Precomputed benchmark results are not bundled in this build. ",
            "Generate them with `inst/scripts/run_vignette_benchmarks.R` ",
            "(see that script's header)._\n", sep = "")
        return(invisible(NULL))
    }
    knitr::kable(timings_wide(results, scenario), digits = 3, caption = caption,
                 col.names = c("Operation", "GRanges (s)", "DuckDBGRanges (s)",
                               "Speedup (x)", "Winner"))
}

## Render the memory-footprint table (both scenarios), or NULL if absent.
make_memory_table <- function(caption = NULL, results = load_vignette_timings()) {
    mem <- if (is.null(results)) NULL else attr(results, "memory")
    if (is.null(mem)) {
        return(invisible(NULL))
    }
    scen <- unique(mem$Scenario)
    out <- data.frame(Scenario = scen, stringsAsFactors = FALSE)
    for (b in .BACKENDS) {
        out[[b]] <- vapply(scen, function(s) {
            v <- mem$MB[mem$Scenario == s & mem$Backend == b]
            if (length(v) == 0L) NA_real_ else v[[1L]]
        }, numeric(1))
    }
    out[["Ratio"]] <- round(out[["GRanges"]] / out[["DuckDBGRanges"]])
    knitr::kable(out, digits = 3, caption = caption,
                 col.names = c("Scenario", "GRanges (MB)", "DuckDBGRanges (MB)",
                               "Smaller by (x)"))
}

## Emit the run configuration as a Markdown paragraph (for reproducibility).
timings_config_note <- function(results = load_vignette_timings()) {
    cfg <- if (is.null(results)) NULL else attr(results, "config")
    if (is.null(cfg)) {
        return(invisible(NULL))
    }
    cat(sprintf(paste0(
        "Configuration: scATAC-seq scenario = %s peaks, variant scenario = %s ",
        "variants, %d-core budget. **GRanges:** in-memory, single-threaded. ",
        "**DuckDBGRanges:** on-disk Parquet; DuckDB autotunes threads up to the ",
        "core budget.\n"),
        format(cfg$scatac_peaks, big.mark = ","),
        format(cfg$variants, big.mark = ","), cfg$cores))
}
