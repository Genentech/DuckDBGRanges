# Tests the basic functions of a DuckDBGRanges.
# library(testthat); library(DuckDBGRanges); source("setup.R"); source("test-DuckDBGRanges.R")

library(GenomicRanges)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test data factories (adapted from GenomicRanges test_findOverlaps-methods.R)
###

make_test_subject <- function() {
    GRanges(Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
            IRanges(1:10, end = 10),
            Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
            seqinfo = Seqinfo(paste0("chr", 1:3)),
            score = 1:10, GC = seq(1, 0, length = 10))
}

make_test_query <- function() {
    # Adapted from GRangesList to GRanges for DuckDBGRanges compatibility
    GRanges(c("chr1:5-10:+",   # nomatch (different strand from chr1 ranges)
              "chr3:2-7:-",    # onematch
              "chr1:1-5:-"))   # twomatch
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor tests
###

test_that("DuckDBGRanges constructor works as expected", {
    seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")

    # start only
    expected <- GRanges(granges_df[["seqnames"]], ranges = IRanges(granges_df[["start"]]))
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", start = "start")
    checkDuckDBGRanges(object, expected)

    # start and end
    expected <- GRanges(granges_df[["seqnames"]], ranges = IRanges(start = granges_df[["start"]], end = granges_df[["end"]], names = granges_df[["id"]]))
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", start = "start", end = "end", keycol = "id")
    checkDuckDBGRanges(object, expected)

    # start and end with mcols
    expected <- GRanges(granges_df[["seqnames"]], ranges = IRanges(start = granges_df[["start"]], end = granges_df[["end"]], names = granges_df[["id"]]),
                        score = granges_df[["score"]], GC = granges_df[["GC"]])
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", start = "start", end = "end", mcols = c("score", "GC"), keycol = "id")
    checkDuckDBGRanges(object, expected)

    # start and width
    expected <- GRanges(granges_df[["seqnames"]], ranges = IRanges(start = granges_df[["start"]], width = granges_df[["width"]], names = granges_df[["id"]]))
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", start = "start", width = "width", keycol = "id")
    checkDuckDBGRanges(object, expected)

    # start and width with mcols
    expected <- GRanges(granges_df[["seqnames"]], ranges = IRanges(start = granges_df[["start"]], width = granges_df[["width"]], names = granges_df[["id"]]),
                        score = granges_df[["score"]], GC = granges_df[["GC"]])
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", start = "start", width = "width", mcols = c("score", "GC"), keycol = "id")
    checkDuckDBGRanges(object, expected)

    # end and width
    expected <- GRanges(granges_df[["seqnames"]], ranges = IRanges(end = granges_df[["end"]], width = granges_df[["width"]], names = granges_df[["id"]]))
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", end = "end", width = "width", keycol = "id")
    checkDuckDBGRanges(object, expected)

    # end and width with mcols
    expected <- GRanges(granges_df[["seqnames"]], ranges = IRanges(end = granges_df[["end"]], width = granges_df[["width"]], names = granges_df[["id"]]),
                        score = granges_df[["score"]], GC = granges_df[["GC"]])
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", end = "end", width = "width", mcols = c("score", "GC"), keycol = "id")
    checkDuckDBGRanges(object, expected)
})

test_that("coersion to a GRanges works for a DuckDBGRanges", {
    seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")

    # seqinfo
    expected <- GRanges(granges_df[["seqnames"]],
                        ranges = IRanges(start = granges_df[["start"]], width = granges_df[["width"]], names = granges_df[["id"]]),
                        seqinfo = seqinfo)
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", start = "start", width = "width",
                            keycol = list(id = granges_df[["id"]]), seqinfo = seqinfo)
    expect_identical(as(object, "GRanges"), expected)

    # strand and seqinfo
    expected <- GRanges(granges_df[["seqnames"]],
                        ranges = IRanges(start = granges_df[["start"]], width = granges_df[["width"]], names = granges_df[["id"]]),
                        strand = granges_df[["strand"]], seqinfo = seqinfo)
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", start = "start", width = "width", strand = "strand",
                            keycol = list(id = granges_df[["id"]]), seqinfo = seqinfo)
    expect_identical(as(object, "GRanges"), expected)

    # strand, mcols, and seqinfo
    expected <- GRanges(granges_df[["seqnames"]],
                        ranges = IRanges(start = granges_df[["start"]], width = granges_df[["width"]], names = granges_df[["id"]]),
                        strand = granges_df[["strand"]], score = granges_df[["score"]], GC = granges_df[["GC"]], seqinfo = seqinfo)
    object <- DuckDBGRanges(granges_tf, seqnames = "seqnames", start = "start", width = "width", strand = "strand",
                            mcols = c("score", "GC"), keycol = list(id = granges_df[["id"]]), seqinfo = seqinfo)
    expect_identical(as(object, "GRanges"), expected)
})

test_that("subsetting after logical filtering works correctly", {
    # This test verifies the fix for the subsetting bug where index-based
    # subsetting (e.g., ddb[1:3]) after logical filtering (e.g., ddb[seqnames(ddb) == "chr22"])
    # would return empty results because row_number wasn't recomputed after filtering.
    
    df <- data.frame(
        seqnames = c(rep("chr1", 5), rep("chr22", 5)),
        start = c(100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000),
        end = c(150, 250, 350, 450, 550, 1050, 2050, 3050, 4050, 5050),
        strand = "*"
    )
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    
    ddb <- DuckDBGRanges(tf, seqnames = "seqnames", start = "start", 
                         end = "end", strand = "strand")
    
    # Original should work
    expect_equal(length(ddb), 10L)
    expect_equal(length(ddb[1:3]), 3L)
    
    # Filter to chr22
    ddb_filtered <- ddb[seqnames(ddb) == "chr22"]
    expect_equal(length(ddb_filtered), 5L)
    
    # Subsetting after filtering should now work (was returning 0 before fix)
    ddb_subset <- ddb_filtered[1:3]
    expect_equal(length(ddb_subset), 3L)
    
    # Verify correct content
    gr_subset <- as(ddb_subset, "GRanges")
    expect_equal(start(gr_subset), c(1000, 2000, 3000))
    
    # Multiple chained filters should also work
    ddb_multi <- ddb[seqnames(ddb) == "chr22"]
    ddb_multi <- ddb_multi[start(ddb_multi) >= 2000]
    expect_equal(length(ddb_multi), 4L)
    gr_multi <- as(ddb_multi, "GRanges")
    expect_equal(start(gr_multi), c(2000, 3000, 4000, 5000))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Overlap tests (adapted from GenomicRanges test_findOverlaps-methods.R)
###

# Helper to create DuckDBGRanges from GRanges for testing
.gr_to_ddb <- function(gr, keycol = NULL) {
    df <- data.frame(
        seqnames = as.character(seqnames(gr)),
        start = start(gr),
        end = end(gr),
        strand = as.character(strand(gr)),
        stringsAsFactors = FALSE
    )
    if (!is.null(mcols(gr)) && ncol(mcols(gr)) > 0) {
        df <- cbind(df, as.data.frame(mcols(gr)))
    }
    if (!is.null(keycol)) {
        df$id <- keycol
    }
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    mcol_names <- if (!is.null(mcols(gr)) && ncol(mcols(gr)) > 0) names(mcols(gr)) else NULL
    if (!is.null(keycol)) {
        DuckDBGRanges(tf, seqnames = "seqnames", start = "start", end = "end",
                      strand = "strand", mcols = mcol_names, keycol = "id")
    } else {
        DuckDBGRanges(tf, seqnames = "seqnames", start = "start", end = "end",
                      strand = "strand", mcols = mcol_names)
    }
}

test_that("findOverlaps with no overlaps returns empty matches", {
    # Adapted from test_findOverlaps_no_overlaps_returns_empty_matches
    query <- make_test_query()
    subject <- make_test_subject()
    ranges(subject) <- shift(ranges(subject), 1000L)

    ddb_subject <- .gr_to_ddb(subject, keycol = seq_len(length(subject)))
    gr_subject <- as(ddb_subject, "GRanges")

    # findOverlaps select = "all"
    hits_ddb <- findOverlaps(query, ddb_subject, ignore.strand = TRUE)
    hits_gr <- findOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(length(hits_ddb), 0L)
    expect_equal(length(hits_ddb), length(hits_gr))

    # countOverlaps
    counts_ddb <- countOverlaps(query, ddb_subject, ignore.strand = TRUE)
    counts_gr <- countOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(counts_ddb, counts_gr)
    expect_equal(counts_ddb, rep(0L, length(query)))

    # subsetByOverlaps
    subset_ddb <- subsetByOverlaps(query, ddb_subject, ignore.strand = TRUE)
    subset_gr <- subsetByOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(subset_ddb, subset_gr)
    expect_equal(length(subset_ddb), 0L)

    # findOverlaps select = "first"
    hits_ddb_first <- findOverlaps(query, ddb_subject, ignore.strand = TRUE, select = "first")
    hits_gr_first <- findOverlaps(query, gr_subject, ignore.strand = TRUE, select = "first")
    expect_equal(hits_ddb_first, hits_gr_first)
    expect_equal(hits_ddb_first, rep(NA_integer_, length(query)))
})

test_that("findOverlaps with empty query", {
    # Adapted from test_findOverlaps_empty_query
    query <- GRanges()
    subject <- make_test_subject()

    ddb_subject <- .gr_to_ddb(subject, keycol = seq_len(length(subject)))
    gr_subject <- as(ddb_subject, "GRanges")

    # findOverlaps select = "all"
    hits_ddb <- findOverlaps(query, ddb_subject, ignore.strand = TRUE)
    hits_gr <- findOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(length(hits_ddb), length(hits_gr))
    expect_equal(queryLength(hits_ddb), queryLength(hits_gr))
    expect_equal(subjectLength(hits_ddb), subjectLength(hits_gr))

    # countOverlaps
    counts_ddb <- countOverlaps(query, ddb_subject, ignore.strand = TRUE)
    counts_gr <- countOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(counts_ddb, counts_gr)

    # subsetByOverlaps
    subset_ddb <- subsetByOverlaps(query, ddb_subject, ignore.strand = TRUE)
    subset_gr <- subsetByOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(subset_ddb, subset_gr)
})

test_that("findOverlaps with empty subject", {
    # Adapted from test_findOverlaps_empty_subject
    query <- make_test_query()
    subject <- GRanges()

    ddb_subject <- .gr_to_ddb(subject)
    gr_subject <- as(ddb_subject, "GRanges")

    # findOverlaps select = "all"
    hits_ddb <- findOverlaps(query, ddb_subject, ignore.strand = TRUE)
    hits_gr <- findOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(length(hits_ddb), length(hits_gr))
    expect_equal(queryLength(hits_ddb), queryLength(hits_gr))
    expect_equal(subjectLength(hits_ddb), subjectLength(hits_gr))

    # countOverlaps
    counts_ddb <- countOverlaps(query, ddb_subject, ignore.strand = TRUE)
    counts_gr <- countOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(counts_ddb, counts_gr)

    # subsetByOverlaps
    subset_ddb <- subsetByOverlaps(query, ddb_subject, ignore.strand = TRUE)
    subset_gr <- subsetByOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(subset_ddb, subset_gr)
})

test_that("findOverlaps with zero, one, and two matches", {
    # Adapted from test_findOverlaps_zero_one_two_matches
    query <- make_test_query()
    subject <- make_test_subject()

    ddb_subject <- .gr_to_ddb(subject, keycol = seq_len(length(subject)))
    gr_subject <- as(ddb_subject, "GRanges")

    # findOverlaps select = "all" with ignore.strand = TRUE
    hits_ddb <- findOverlaps(query, ddb_subject, ignore.strand = TRUE)
    hits_gr <- findOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(queryHits(hits_ddb), queryHits(hits_gr))
    expect_equal(subjectHits(hits_ddb), subjectHits(hits_gr))

    # countOverlaps
    counts_ddb <- countOverlaps(query, ddb_subject, ignore.strand = TRUE)
    counts_gr <- countOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(counts_ddb, counts_gr)

    # subsetByOverlaps
    subset_ddb <- subsetByOverlaps(query, ddb_subject, ignore.strand = TRUE)
    subset_gr <- subsetByOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(subset_ddb, subset_gr)

    # findOverlaps select = "first"
    hits_ddb_first <- findOverlaps(query, ddb_subject, ignore.strand = TRUE, select = "first")
    hits_gr_first <- findOverlaps(query, gr_subject, ignore.strand = TRUE, select = "first")
    expect_equal(hits_ddb_first, hits_gr_first)
})

test_that("findOverlaps respects strand matching", {
    # Adapted from test_findOverlaps_either_strand
    query <- make_test_query()
    subject <- make_test_subject()

    ddb_subject <- .gr_to_ddb(subject, keycol = seq_len(length(subject)))
    gr_subject <- as(ddb_subject, "GRanges")

    # With strand matching (default ignore.strand = FALSE)
    hits_ddb <- findOverlaps(query, ddb_subject)
    hits_gr <- findOverlaps(query, gr_subject)
    expect_equal(queryHits(hits_ddb), queryHits(hits_gr))
    expect_equal(subjectHits(hits_ddb), subjectHits(hits_gr))

    counts_ddb <- countOverlaps(query, ddb_subject)
    counts_gr <- countOverlaps(query, gr_subject)
    expect_equal(counts_ddb, counts_gr)

    subset_ddb <- subsetByOverlaps(query, ddb_subject)
    subset_gr <- subsetByOverlaps(query, gr_subject)
    expect_equal(subset_ddb, subset_gr)

    # With ignore.strand = TRUE
    hits_ddb_ignore <- findOverlaps(query, ddb_subject, ignore.strand = TRUE)
    hits_gr_ignore <- findOverlaps(query, gr_subject, ignore.strand = TRUE)
    expect_equal(queryHits(hits_ddb_ignore), queryHits(hits_gr_ignore))
    expect_equal(subjectHits(hits_ddb_ignore), subjectHits(hits_gr_ignore))

    # Verify that ignoring strand gives more (or equal) matches
    expect_true(length(hits_ddb_ignore) >= length(hits_ddb))
})

test_that("findOverlaps with minoverlap parameter", {
    # Adapted from test_findOverlaps_minoverlap_GRanges_GRangesList
    query <- make_test_query()
    subject <- make_test_subject()

    ddb_subject <- .gr_to_ddb(subject, keycol = seq_len(length(subject)))
    gr_subject <- as(ddb_subject, "GRanges")

    # Test various minoverlap values
    for (minoverlap in c(1L, 3L, 5L, 6L)) {
        hits_ddb <- findOverlaps(query, ddb_subject, minoverlap = minoverlap, ignore.strand = TRUE)
        hits_gr <- findOverlaps(query, gr_subject, minoverlap = minoverlap, ignore.strand = TRUE)
        expect_equal(queryHits(hits_ddb), queryHits(hits_gr),
                     info = paste("minoverlap =", minoverlap))
        expect_equal(subjectHits(hits_ddb), subjectHits(hits_gr),
                     info = paste("minoverlap =", minoverlap))
    }

    # countOverlaps with minoverlap
    counts_ddb <- countOverlaps(query, ddb_subject, minoverlap = 5L, ignore.strand = TRUE)
    counts_gr <- countOverlaps(query, gr_subject, minoverlap = 5L, ignore.strand = TRUE)
    expect_equal(counts_ddb, counts_gr)
})

test_that("findOverlaps with maxgap parameter", {
    subject <- GRanges("chr1", IRanges(100, 150), "+")
    query <- GRanges("chr1", IRanges(161, 200), "+")

    ddb_subject <- .gr_to_ddb(subject)
    gr_subject <- as(ddb_subject, "GRanges")

    # Default maxgap = -1 (must overlap)
    hits_ddb_default <- findOverlaps(query, ddb_subject)
    hits_gr_default <- findOverlaps(query, gr_subject)
    expect_equal(length(hits_ddb_default), length(hits_gr_default))
    expect_equal(length(hits_ddb_default), 0L)

    # maxgap = 10 (should find overlap since gap is 10)
    hits_ddb_10 <- findOverlaps(query, ddb_subject, maxgap = 10L)
    hits_gr_10 <- findOverlaps(query, gr_subject, maxgap = 10L)
    expect_equal(length(hits_ddb_10), length(hits_gr_10))
    expect_equal(length(hits_ddb_10), 1L)

    # maxgap = 5 (should not find overlap since gap is 10)
    hits_ddb_5 <- findOverlaps(query, ddb_subject, maxgap = 5L)
    hits_gr_5 <- findOverlaps(query, gr_subject, maxgap = 5L)
    expect_equal(length(hits_ddb_5), length(hits_gr_5))
    expect_equal(length(hits_ddb_5), 0L)
})

test_that("findOverlaps(GRanges, DuckDBGRanges) select options", {
    subject_df <- data.frame(
        seqnames = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100L, 140L, 150L, 300L),
        end = c(150L, 250L, 200L, 350L),
        strand = c("+", "-", "+", "-")
    )
    subject_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(subject_df, subject_tf)
    ddb_subject <- DuckDBGRanges(subject_tf, seqnames = "seqnames",
                                 start = "start", end = "end", strand = "strand")
    gr_subject <- as(ddb_subject, "GRanges")

    query <- GRanges(c("chr1:120-180:+", "chr2:160-220:+", "chr3:1-50:*"))

    # select = "all"
    hits_ddb_all <- findOverlaps(query, ddb_subject, ignore.strand = TRUE, select = "all")
    hits_gr_all <- findOverlaps(query, gr_subject, ignore.strand = TRUE, select = "all")
    expect_s4_class(hits_ddb_all, "Hits")
    expect_equal(queryHits(hits_ddb_all), queryHits(hits_gr_all))
    expect_equal(subjectHits(hits_ddb_all), subjectHits(hits_gr_all))

    # select = "first"
    hits_ddb_first <- findOverlaps(query, ddb_subject, ignore.strand = TRUE, select = "first")
    hits_gr_first <- findOverlaps(query, gr_subject, ignore.strand = TRUE, select = "first")
    expect_equal(hits_ddb_first, hits_gr_first)

    # select = "last"
    hits_ddb_last <- findOverlaps(query, ddb_subject, ignore.strand = TRUE, select = "last")
    hits_gr_last <- findOverlaps(query, gr_subject, ignore.strand = TRUE, select = "last")
    expect_equal(hits_ddb_last, hits_gr_last)

    # select = "arbitrary"
    hits_ddb_arb <- findOverlaps(query, ddb_subject, ignore.strand = TRUE, select = "arbitrary")
    hits_gr_arb <- findOverlaps(query, gr_subject, ignore.strand = TRUE, select = "arbitrary")
    # For arbitrary, just check that NA positions match
    expect_equal(is.na(hits_ddb_arb), is.na(hits_gr_arb))

    unlink(subject_tf)
})

test_that("findOverlaps(DuckDBGRanges, DuckDBGRanges) matches GRanges results", {
    query_df <- data.frame(
        seqnames = c("chr1", "chr2", "chr3"),
        start = c(120L, 160L, 500L),
        end = c(180L, 220L, 600L),
        strand = c("+", "+", "+")
    )
    query_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(query_df, query_tf)
    ddb_query <- DuckDBGRanges(query_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")
    gr_query <- as(ddb_query, "GRanges")

    subject_df <- data.frame(
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 140L, 150L),
        end = c(150L, 250L, 200L),
        strand = c("+", "-", "+")
    )
    subject_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(subject_df, subject_tf)
    ddb_subject <- DuckDBGRanges(subject_tf, seqnames = "seqnames",
                                 start = "start", end = "end", strand = "strand")
    gr_subject <- as(ddb_subject, "GRanges")

    # ignore.strand = TRUE
    hits_ddb <- findOverlaps(ddb_query, ddb_subject, ignore.strand = TRUE)
    hits_gr <- findOverlaps(gr_query, gr_subject, ignore.strand = TRUE)
    expect_s4_class(hits_ddb, "Hits")
    expect_equal(queryHits(hits_ddb), queryHits(hits_gr))
    expect_equal(subjectHits(hits_ddb), subjectHits(hits_gr))

    # ignore.strand = FALSE
    hits_ddb_strand <- findOverlaps(ddb_query, ddb_subject)
    hits_gr_strand <- findOverlaps(gr_query, gr_subject)
    expect_equal(queryHits(hits_ddb_strand), queryHits(hits_gr_strand))
    expect_equal(subjectHits(hits_ddb_strand), subjectHits(hits_gr_strand))

    unlink(c(query_tf, subject_tf))
})

test_that("findOverlaps(DuckDBGRanges, GRanges) matches GRanges results", {
    query_df <- data.frame(
        seqnames = c("chr1", "chr2"),
        start = c(120L, 160L),
        end = c(180L, 220L),
        strand = c("+", "+")
    )
    query_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(query_df, query_tf)
    ddb_query <- DuckDBGRanges(query_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")
    gr_query <- as(ddb_query, "GRanges")

    subject <- GRanges(c("chr1:100-150:+", "chr1:140-250:-", "chr2:150-200:+"))

    hits_ddb <- findOverlaps(ddb_query, subject, ignore.strand = TRUE)
    hits_gr <- findOverlaps(gr_query, subject, ignore.strand = TRUE)
    expect_s4_class(hits_ddb, "Hits")
    expect_equal(queryHits(hits_ddb), queryHits(hits_gr))
    expect_equal(subjectHits(hits_ddb), subjectHits(hits_gr))

    unlink(query_tf)
})

test_that("countOverlaps(GRanges, DuckDBGRanges) matches GRanges results", {
    subject_df <- data.frame(
        seqnames = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100L, 140L, 150L, 300L),
        end = c(150L, 250L, 200L, 350L),
        strand = c("+", "+", "+", "+")
    )
    subject_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(subject_df, subject_tf)
    ddb_subject <- DuckDBGRanges(subject_tf, seqnames = "seqnames",
                                 start = "start", end = "end", strand = "strand")
    gr_subject <- as(ddb_subject, "GRanges")

    query <- GRanges(c("chr1:120-180", "chr2:160-220", "chr3:1-100"))

    counts_ddb <- countOverlaps(query, ddb_subject)
    counts_gr <- countOverlaps(query, gr_subject)
    expect_equal(counts_ddb, counts_gr)

    unlink(subject_tf)
})

test_that("countOverlaps(DuckDBGRanges, DuckDBGRanges) matches GRanges results", {
    query_df <- data.frame(
        seqnames = c("chr1", "chr2", "chr3"),
        start = c(120L, 160L, 500L),
        end = c(180L, 220L, 600L),
        strand = c("+", "+", "+")
    )
    query_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(query_df, query_tf)
    ddb_query <- DuckDBGRanges(query_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")
    gr_query <- as(ddb_query, "GRanges")

    subject_df <- data.frame(
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 140L, 150L),
        end = c(150L, 250L, 200L),
        strand = c("+", "+", "+")
    )
    subject_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(subject_df, subject_tf)
    ddb_subject <- DuckDBGRanges(subject_tf, seqnames = "seqnames",
                                 start = "start", end = "end", strand = "strand")
    gr_subject <- as(ddb_subject, "GRanges")

    counts_ddb <- countOverlaps(ddb_query, ddb_subject)
    counts_gr <- countOverlaps(gr_query, gr_subject)
    expect_equal(counts_ddb, counts_gr)

    unlink(c(query_tf, subject_tf))
})

test_that("subsetByOverlaps(GRanges, DuckDBGRanges) matches GRanges results", {
    subject_df <- data.frame(
        seqnames = c("chr1", "chr2"),
        start = c(100L, 150L),
        end = c(150L, 200L),
        strand = c("+", "+")
    )
    subject_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(subject_df, subject_tf)
    ddb_subject <- DuckDBGRanges(subject_tf, seqnames = "seqnames",
                                 start = "start", end = "end", strand = "strand")
    gr_subject <- as(ddb_subject, "GRanges")

    query <- GRanges(c("chr1:120-180", "chr2:160-220", "chr3:1-100"))

    # invert = FALSE
    result_ddb <- subsetByOverlaps(query, ddb_subject)
    result_gr <- subsetByOverlaps(query, gr_subject)
    expect_s4_class(result_ddb, "GRanges")
    expect_equal(result_ddb, result_gr)

    # invert = TRUE
    result_ddb_inv <- subsetByOverlaps(query, ddb_subject, invert = TRUE)
    result_gr_inv <- subsetByOverlaps(query, gr_subject, invert = TRUE)
    expect_equal(result_ddb_inv, result_gr_inv)

    unlink(subject_tf)
})

test_that("subsetByOverlaps(DuckDBGRanges, DuckDBGRanges) returns DuckDBGRanges", {
    query_df <- data.frame(
        id = c("q1", "q2", "q3", "q4"),
        seqnames = c("chr1", "chr1", "chr2", "chr3"),
        start = c(100L, 500L, 150L, 1000L),
        end = c(150L, 550L, 200L, 1050L),
        strand = c("+", "+", "+", "+")
    )
    query_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(query_df, query_tf)
    ddb_query <- DuckDBGRanges(query_tf, seqnames = "seqnames",
                               start = "start", end = "end",
                               strand = "strand", keycol = "id")
    gr_query <- as(ddb_query, "GRanges")

    subject_df <- data.frame(
        seqnames = c("chr1", "chr2"),
        start = c(120L, 180L),
        end = c(180L, 250L),
        strand = c("+", "+")
    )
    subject_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(subject_df, subject_tf)
    ddb_subject <- DuckDBGRanges(subject_tf, seqnames = "seqnames",
                                 start = "start", end = "end", strand = "strand")
    gr_subject <- as(ddb_subject, "GRanges")

    # invert = FALSE
    object <- subsetByOverlaps(ddb_query, ddb_subject)
    expected <- subsetByOverlaps(gr_query, gr_subject)
    checkDuckDBGRanges(object, expected)

    # invert = TRUE
    object_inv <- subsetByOverlaps(ddb_query, ddb_subject, invert = TRUE)
    expected_inv <- subsetByOverlaps(gr_query, gr_subject, invert = TRUE)
    checkDuckDBGRanges(object_inv, expected_inv)

    unlink(c(query_tf, subject_tf))
})

test_that("subsetByOverlaps(DuckDBGRanges, GRanges) returns DuckDBGRanges", {
    query_df <- data.frame(
        id = c("q1", "q2", "q3"),
        seqnames = c("chr1", "chr2", "chr3"),
        start = c(100L, 150L, 500L),
        end = c(150L, 200L, 550L),
        strand = c("+", "+", "+")
    )
    query_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(query_df, query_tf)
    ddb_query <- DuckDBGRanges(query_tf, seqnames = "seqnames",
                               start = "start", end = "end",
                               strand = "strand", keycol = "id")
    gr_query <- as(ddb_query, "GRanges")

    subject <- GRanges(c("chr1:120-180:+", "chr2:180-250:+"))

    object <- subsetByOverlaps(ddb_query, subject)
    expected <- subsetByOverlaps(gr_query, subject)
    checkDuckDBGRanges(object, expected)

    # invert = TRUE
    object_inv <- subsetByOverlaps(ddb_query, subject, invert = TRUE)
    expected_inv <- subsetByOverlaps(gr_query, subject, invert = TRUE)
    checkDuckDBGRanges(object_inv, expected_inv)

    unlink(query_tf)
})

test_that("subsetByOverlaps works with row_number keycols", {
    query_df <- data.frame(
        seqnames = c("chr1", "chr1", "chr2", "chr3"),
        start = c(100L, 500L, 150L, 1000L),
        end = c(150L, 550L, 200L, 1050L),
        strand = c("+", "+", "+", "+")
    )
    query_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(query_df, query_tf)
    ddb_query <- DuckDBGRanges(query_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")
    gr_query <- as(ddb_query, "GRanges")

    subject_df <- data.frame(
        seqnames = c("chr1", "chr2"),
        start = c(120L, 180L),
        end = c(180L, 250L),
        strand = c("+", "+")
    )
    subject_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(subject_df, subject_tf)
    ddb_subject <- DuckDBGRanges(subject_tf, seqnames = "seqnames",
                                 start = "start", end = "end", strand = "strand")
    gr_subject <- as(ddb_subject, "GRanges")

    # invert = FALSE
    object <- subsetByOverlaps(ddb_query, ddb_subject)
    expected <- subsetByOverlaps(gr_query, gr_subject)
    checkDuckDBGRanges(object, expected)

    # invert = TRUE
    object_inv <- subsetByOverlaps(ddb_query, ddb_subject, invert = TRUE)
    expected_inv <- subsetByOverlaps(gr_query, gr_subject, invert = TRUE)
    checkDuckDBGRanges(object_inv, expected_inv)

    unlink(c(query_tf, subject_tf))
})

test_that("overlap methods work with granges_df test data from setup.R", {
    ddb_gr <- DuckDBGRanges(granges_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    query <- GRanges(c("chr1:1-5:+", "chr2:3-8:-", "chr3:7-12:*"))

    # findOverlaps
    hits_ddb <- findOverlaps(query, ddb_gr, ignore.strand = TRUE)
    hits_gr <- findOverlaps(query, gr, ignore.strand = TRUE)
    expect_s4_class(hits_ddb, "Hits")
    expect_equal(countQueryHits(hits_ddb), countQueryHits(hits_gr))
    expect_equal(length(hits_ddb), length(hits_gr))

    # countOverlaps
    counts_ddb <- countOverlaps(query, ddb_gr, ignore.strand = TRUE)
    counts_gr <- countOverlaps(query, gr, ignore.strand = TRUE)
    expect_equal(counts_ddb, counts_gr)

    # subsetByOverlaps
    subset_ddb <- subsetByOverlaps(query, ddb_gr, ignore.strand = TRUE)
    subset_gr <- subsetByOverlaps(query, gr, ignore.strand = TRUE)
    expect_equal(subset_ddb, subset_gr)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Intra-range method tests
###

test_that("shift works for DuckDBGRanges", {
    ddb_gr <- DuckDBGRanges(granges_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # shift by positive amount
    result_ddb <- shift(ddb_gr, shift = 10L)
    result_gr <- shift(gr, shift = 10L)
    checkDuckDBGRanges(result_ddb, result_gr)

    # shift by negative amount
    result_ddb_neg <- shift(ddb_gr, shift = -5L)
    result_gr_neg <- shift(gr, shift = -5L)
    checkDuckDBGRanges(result_ddb_neg, result_gr_neg)
})

test_that("narrow works for DuckDBGRanges", {
    ddb_gr <- DuckDBGRanges(granges_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # narrow by start
    result_ddb <- narrow(ddb_gr, start = 2)
    result_gr <- narrow(gr, start = 2)
    checkDuckDBGRanges(result_ddb, result_gr)

    # narrow by end
    result_ddb_end <- narrow(ddb_gr, end = -2)
    result_gr_end <- narrow(gr, end = -2)
    checkDuckDBGRanges(result_ddb_end, result_gr_end)

    # narrow by NEGATIVE start: counts back from the range end (base solveUserSEW),
    # not start + (start - 1). Oracle-checked against base IRanges::narrow.
    result_ddb_ns <- narrow(ddb_gr, start = -1)
    result_gr_ns <- narrow(gr, start = -1)
    checkDuckDBGRanges(result_ddb_ns, result_gr_ns)
})

test_that("resize works for DuckDBGRanges", {
    ddb_gr <- DuckDBGRanges(granges_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # resize anchored at start
    result_ddb <- resize(ddb_gr, width = 5L, fix = "start")
    result_gr <- resize(gr, width = 5L, fix = "start")
    checkDuckDBGRanges(result_ddb, result_gr)

    # resize anchored at end
    result_ddb_end <- resize(ddb_gr, width = 5L, fix = "end")
    result_gr_end <- resize(gr, width = 5L, fix = "end")
    checkDuckDBGRanges(result_ddb_end, result_gr_end)

    # resize with ignore.strand
    result_ddb_ign <- resize(ddb_gr, width = 5L, fix = "start", ignore.strand = TRUE)
    result_gr_ign <- resize(gr, width = 5L, fix = "start", ignore.strand = TRUE)
    checkDuckDBGRanges(result_ddb_ign, result_gr_ign)

    # resize anchored at center must match base IRanges for both parities: an odd
    result_ddb_c_odd <- resize(ddb_gr, width = 5L, fix = "center")
    result_gr_c_odd <- resize(gr, width = 5L, fix = "center")
    checkDuckDBGRanges(result_ddb_c_odd, result_gr_c_odd)

    result_ddb_c_even <- resize(ddb_gr, width = 4L, fix = "center")
    result_gr_c_even <- resize(gr, width = 4L, fix = "center")
    checkDuckDBGRanges(result_ddb_c_even, result_gr_c_even)
})

test_that("flank works for DuckDBGRanges", {
    ddb_gr <- DuckDBGRanges(granges_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # upstream flanking regions
    result_ddb <- flank(ddb_gr, width = 100L, start = TRUE)
    result_gr <- flank(gr, width = 100L, start = TRUE)
    checkDuckDBGRanges(result_ddb, result_gr)

    # downstream flanking regions
    result_ddb_down <- flank(ddb_gr, width = 100L, start = FALSE)
    result_gr_down <- flank(gr, width = 100L, start = FALSE)
    checkDuckDBGRanges(result_ddb_down, result_gr_down)

    # both sides
    result_ddb_both <- flank(ddb_gr, width = 50L, both = TRUE)
    result_gr_both <- flank(gr, width = 50L, both = TRUE)
    checkDuckDBGRanges(result_ddb_both, result_gr_both)
})

test_that("promoters works for DuckDBGRanges", {
    ddb_gr <- DuckDBGRanges(granges_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- promoters(ddb_gr, upstream = 500, downstream = 100)
    result_gr <- promoters(gr, upstream = 500, downstream = 100)
    checkDuckDBGRanges(result_ddb, result_gr)
})

test_that("terminators works for DuckDBGRanges", {
    ddb_gr <- DuckDBGRanges(granges_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- terminators(ddb_gr, upstream = 500, downstream = 100)
    result_gr <- terminators(gr, upstream = 500, downstream = 100)
    checkDuckDBGRanges(result_ddb, result_gr)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Inter-range method tests
###

test_that("range works for DuckDBGRanges", {
    ddb_gr <- DuckDBGRanges(granges_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # range with strand
    result_ddb <- range(ddb_gr)
    result_gr <- range(gr)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    expect_equal(length(result_ddb), length(result_gr))
    # Compare values after sorting since aggregated results may have different order
    expect_setequal(as.vector(start(result_ddb)), start(result_gr))
    expect_setequal(as.vector(end(result_ddb)), end(result_gr))

    # range ignoring strand
    result_ddb_ign <- range(ddb_gr, ignore.strand = TRUE)
    result_gr_ign <- range(gr, ignore.strand = TRUE)
    expect_equal(length(result_ddb_ign), length(result_gr_ign))
    expect_setequal(as.vector(start(result_ddb_ign)), start(result_gr_ign))
    expect_setequal(as.vector(end(result_ddb_ign)), end(result_gr_ign))
})

test_that("range works for empty DuckDBGRanges", {
    # Create empty DuckDBGRanges
    empty_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0)
    )
    empty_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(empty_df, empty_tf)
    ddb_empty <- DuckDBGRanges(empty_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")

    result <- range(ddb_empty)
    expect_s4_class(result, "DuckDBGRanges")
    expect_equal(length(result), 0L)

    unlink(empty_tf)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparison method tests
###

test_that("duplicated works for DuckDBGRanges", {
    # Create test data with duplicates
    dup_df <- data.frame(
        id = c("a", "b", "c", "d", "e"),
        seqnames = c("chr1", "chr1", "chr1", "chr2", "chr1"),
        start = c(100L, 100L, 200L, 100L, 100L),
        end = c(150L, 150L, 250L, 150L, 150L),
        strand = c("+", "+", "+", "+", "+")
    )
    dup_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(dup_df, dup_tf)
    ddb_gr <- DuckDBGRanges(dup_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # duplicated returns DuckDBColumn
    result_ddb <- duplicated(ddb_gr)
    result_gr <- duplicated(gr)
    expect_s4_class(result_ddb, "DuckDBColumn")
    expect_equal(unname(as.vector(result_ddb)), unname(result_gr))

    # fromLast = TRUE
    result_ddb_last <- duplicated(ddb_gr, fromLast = TRUE)
    result_gr_last <- duplicated(gr, fromLast = TRUE)
    expect_s4_class(result_ddb_last, "DuckDBColumn")
    expect_equal(unname(as.vector(result_ddb_last)), unname(result_gr_last))

    # Test ! operator on result
    not_dup <- !result_ddb
    expect_s4_class(not_dup, "DuckDBColumn")
    expect_equal(unname(as.vector(not_dup)), unname(!result_gr))

    unlink(dup_tf)
})

test_that("unique works for DuckDBGRanges", {
    # Create test data with duplicates
    dup_df <- data.frame(
        id = c("a", "b", "c", "d", "e"),
        seqnames = c("chr1", "chr1", "chr1", "chr2", "chr1"),
        start = c(100L, 100L, 200L, 100L, 100L),
        end = c(150L, 150L, 250L, 150L, 150L),
        strand = c("+", "+", "+", "+", "+")
    )
    dup_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(dup_df, dup_tf)
    ddb_gr <- DuckDBGRanges(dup_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- unique(ddb_gr)
    result_gr <- unique(gr)
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(dup_tf)
})

test_that("match works for DuckDBGRanges", {
    # Create two overlapping sets
    df1 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 200L, 100L),
        end = c(150L, 250L, 150L),
        strand = c("+", "+", "+")
    )
    df2 <- data.frame(
        id = c("x", "y", "z", "w"),
        seqnames = c("chr1", "chr1", "chr2", "chr3"),
        start = c(200L, 100L, 100L, 100L),
        end = c(250L, 150L, 150L, 150L),
        strand = c("+", "+", "+", "+")
    )
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    # match returns DuckDBColumn
    result_ddb <- match(ddb_gr1, ddb_gr2)
    result_gr <- match(gr1, gr2)
    expect_s4_class(result_ddb, "DuckDBColumn")
    expect_equal(unname(as.vector(result_ddb)), unname(result_gr))

    unlink(c(tf1, tf2))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Parallel set operation tests
###

test_that("punion works for DuckDBGRanges", {
    df1 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 200L, 100L),
        end = c(150L, 250L, 150L),
        strand = c("+", "+", "+")
    )
    df2 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(120L, 180L, 120L),
        end = c(180L, 300L, 200L),
        strand = c("+", "+", "+")
    )
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    # Convert to GRanges and ensure same order
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    # Ensure same order by matching on names
    gr2 <- gr2[names(gr1)]

    # DuckDBGRanges, DuckDBGRanges
    result_ddb <- punion(ddb_gr1, ddb_gr2)
    result_gr <- punion(gr1, gr2)
    checkDuckDBGRanges(result_ddb, result_gr)

    # Mixed-type tests: DuckDBGRanges with GRanges falls back to GRanges
    # These require matching order which DuckDB doesn't guarantee,
    # so we reorder GRanges to match DuckDBGRanges order before calling
    # Build GRanges from original data frames to get consistent ordering
    gr2_for_mix <- GRanges(df2$seqnames, IRanges(df2$start, df2$end, names = df2$id),
                           strand = df2$strand)
    gr1_for_mix <- GRanges(df1$seqnames, IRanges(df1$start, df1$end, names = df1$id),
                           strand = df1$strand)

    # DuckDBGRanges, GRanges - reorder gr2 to match ddb_gr1's order
    gr2_ordered <- gr2_for_mix[names(ddb_gr1)]
    result_mixed <- punion(ddb_gr1, gr2_ordered)
    result_gr_mixed <- punion(gr1, gr2_ordered[names(gr1)])
    checkDuckDBGRanges(result_mixed, result_gr_mixed)

    # GRanges, DuckDBGRanges - reorder gr1 to match ddb_gr2's order
    gr1_ordered <- gr1_for_mix[names(ddb_gr2)]
    result_mixed2 <- punion(gr1_ordered, ddb_gr2)
    result_gr_mixed2 <- punion(gr1_ordered[names(gr2)], gr2)
    checkDuckDBGRanges(result_mixed2, result_gr_mixed2)

    unlink(c(tf1, tf2))
})

test_that("pintersect works for DuckDBGRanges", {
    df1 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 200L, 100L),
        end = c(150L, 250L, 150L),
        strand = c("+", "+", "+")
    )
    df2 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(120L, 180L, 120L),
        end = c(180L, 300L, 200L),
        strand = c("+", "+", "+")
    )
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    # Convert to GRanges and ensure same order
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    # Ensure same order by matching on names
    gr2 <- gr2[names(gr1)]

    # DuckDBGRanges, DuckDBGRanges
    result_ddb <- pintersect(ddb_gr1, ddb_gr2)
    result_gr <- pintersect(gr1, gr2)
    checkDuckDBGRanges(result_ddb, result_gr)

    # Mixed-type: DuckDBGRanges with GRanges - now returns DuckDBGRanges
    # Build GRanges from original data to get consistent ordering
    gr2_for_mix <- GRanges(df2$seqnames, IRanges(df2$start, df2$end, names = df2$id),
                           strand = df2$strand)
    gr2_ordered <- gr2_for_mix[names(ddb_gr1)]
    result_mixed <- pintersect(ddb_gr1, gr2_ordered)
    result_gr_mixed <- pintersect(gr1, gr2_ordered[names(gr1)])
    checkDuckDBGRanges(result_mixed, result_gr_mixed)

    unlink(c(tf1, tf2))
})

test_that("psetdiff works for DuckDBGRanges", {
    # psetdiff requires y ranges to be at the start or end of x ranges
    # (not strictly inside), so we use valid test data
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(100L, 200L),
        end = c(200L, 300L),
        strand = c("+", "+")
    )
    df2 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(100L, 200L),
        end = c(150L, 250L),
        strand = c("+", "+")
    )
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    # Convert to GRanges and ensure same order
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    # Ensure same order by matching on names
    gr2 <- gr2[names(gr1)]

    # DuckDBGRanges, DuckDBGRanges
    result_ddb <- psetdiff(ddb_gr1, ddb_gr2)
    result_gr <- psetdiff(gr1, gr2)
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Distance method tests
###

test_that("distance works for DuckDBGRanges", {
    df1 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 200L, 100L),
        end = c(150L, 250L, 150L),
        strand = c("+", "+", "+")
    )
    df2 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(200L, 100L, 200L),
        end = c(250L, 150L, 250L),
        strand = c("+", "+", "+")
    )
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    # Convert to GRanges and ensure same order
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    # Ensure same order by matching on names
    gr2 <- gr2[names(gr1)]

    # DuckDBGRanges, DuckDBGRanges - returns integer vector
    result_ddb <- distance(ddb_gr1, ddb_gr2)
    result_gr <- distance(gr1, gr2)
    expect_type(result_ddb, "integer")
    expect_setequal(result_ddb, result_gr)

    # Mixed-type tests: requires matching order which DuckDB doesn't guarantee
    # Build GRanges from original data to get consistent ordering
    gr1_for_mix <- GRanges(df1$seqnames, IRanges(df1$start, df1$end, names = df1$id),
                           strand = df1$strand)
    gr2_for_mix <- GRanges(df2$seqnames, IRanges(df2$start, df2$end, names = df2$id),
                           strand = df2$strand)

    # DuckDBGRanges, GRanges - reorder gr2 to match ddb_gr1's order
    gr2_ordered <- gr2_for_mix[names(ddb_gr1)]
    result_mixed <- distance(ddb_gr1, gr2_ordered)
    expect_type(result_mixed, "integer")
    expect_equal(length(result_mixed), length(result_gr))

    # GRanges, DuckDBGRanges - reorder gr1 to match ddb_gr2's order
    gr1_ordered <- gr1_for_mix[names(ddb_gr2)]
    result_mixed2 <- distance(gr1_ordered, ddb_gr2)
    expect_type(result_mixed2, "integer")
    expect_equal(length(result_mixed2), length(result_gr))

    unlink(c(tf1, tf2))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ordering method tests
###

test_that("sort works for DuckDBGRanges", {
    # Create unsorted test data
    unsorted_df <- data.frame(
        id = c("a", "b", "c", "d", "e"),
        seqnames = c("chr2", "chr1", "chr1", "chr3", "chr1"),
        start = c(100L, 200L, 100L, 50L, 150L),
        end = c(150L, 250L, 150L, 100L, 200L),
        strand = c("+", "-", "+", "+", "+")
    )
    unsorted_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(unsorted_df, unsorted_tf)
    ddb_gr <- DuckDBGRanges(unsorted_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")

    # sort (default)
    result_ddb <- sort(ddb_gr)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    # Verify result is sorted by comparing starts in sequence
    sorted_starts <- as.vector(start(result_ddb))
    sorted_seqnames <- as.character(seqnames(result_ddb))
    # Check that chr1 ranges come before chr2, chr2 before chr3
    expect_true(all(sorted_seqnames[1:3] == "chr1"))
    expect_true(sorted_seqnames[4] == "chr2")
    expect_true(sorted_seqnames[5] == "chr3")

    # sort decreasing
    result_ddb_dec <- sort(ddb_gr, decreasing = TRUE)
    expect_s4_class(result_ddb_dec, "DuckDBGRanges")
    dec_seqnames <- as.character(seqnames(result_ddb_dec))
    # Decreasing order: chr3 first, then chr2, then chr1
    expect_true(dec_seqnames[1] == "chr3")
    expect_true(dec_seqnames[2] == "chr2")
    expect_true(all(dec_seqnames[3:5] == "chr1"))

    # sort with ignore.strand
    result_ddb_ign <- sort(ddb_gr, ignore.strand = TRUE)
    expect_s4_class(result_ddb_ign, "DuckDBGRanges")
    # Just verify it returns correct type
    expect_equal(length(result_ddb_ign), length(ddb_gr))

    unlink(unsorted_tf)
})

test_that("order works for DuckDBGRanges", {
    # Create test data
    test_df <- data.frame(
        id = c("a", "b", "c", "d"),
        seqnames = c("chr2", "chr1", "chr1", "chr1"),
        start = c(100L, 200L, 100L, 150L),
        end = c(150L, 250L, 150L, 200L),
        strand = c("+", "+", "+", "+")
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)
    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")

    # order - get ordering permutation
    result_ddb <- order(ddb_gr)
    expect_type(result_ddb, "integer")
    expect_equal(length(result_ddb), 4L)
    # Verify result is a valid permutation
    expect_setequal(result_ddb, 1:4)

    # order decreasing
    result_ddb_dec <- order(ddb_gr, decreasing = TRUE)
    expect_type(result_ddb_dec, "integer")
    expect_setequal(result_ddb_dec, 1:4)

    unlink(test_tf)
})

test_that("is.unsorted works for DuckDBGRanges", {
    # Create sorted data (chr1 before chr2, sorted by start)
    sorted_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 200L, 100L),
        end = c(150L, 250L, 150L),
        strand = c("+", "+", "+")
    )
    sorted_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(sorted_df, sorted_tf)
    ddb_sorted <- DuckDBGRanges(sorted_tf, seqnames = "seqnames",
                                start = "start", end = "end", strand = "strand",
                                keycol = "id")

    # Verify method returns correct type and value for sorted data
    result_sorted <- is.unsorted(ddb_sorted)
    expect_type(result_sorted, "logical")
    expect_false(result_sorted)

    # Create unsorted data (chr2 before chr1 - wrong order)
    unsorted_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr2", "chr1", "chr1"),
        start = c(100L, 200L, 100L),
        end = c(150L, 250L, 150L),
        strand = c("+", "+", "+")
    )
    unsorted_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(unsorted_df, unsorted_tf)
    ddb_unsorted <- DuckDBGRanges(unsorted_tf, seqnames = "seqnames",
                                  start = "start", end = "end", strand = "strand",
                                  keycol = "id")

    # Verify method returns correct type and value for unsorted data
    result_unsorted <- is.unsorted(ddb_unsorted)
    expect_type(result_unsorted, "logical")
    expect_true(result_unsorted)

    # Check strictly sorted (data with ties)
    tied_df <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(100L, 100L),
        end = c(150L, 150L),
        strand = c("+", "+")
    )
    tied_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(tied_df, tied_tf)
    ddb_tied <- DuckDBGRanges(tied_tf, seqnames = "seqnames",
                              start = "start", end = "end", strand = "strand",
                              keycol = "id")

    # Non-strict: ties are OK, so NOT unsorted
    result_tied <- is.unsorted(ddb_tied)
    expect_type(result_tied, "logical")
    expect_false(result_tied)

    # Strict: ties mean unsorted
    result_tied_strict <- is.unsorted(ddb_tied, strictly = TRUE)
    expect_type(result_tied_strict, "logical")
    expect_true(result_tied_strict)

    unlink(c(sorted_tf, unsorted_tf, tied_tf))
})

test_that("rank works for DuckDBGRanges", {
    test_df <- data.frame(
        id = c("a", "b", "c", "d"),
        seqnames = c("chr2", "chr1", "chr1", "chr1"),
        start = c(100L, 200L, 100L, 100L),
        end = c(150L, 250L, 150L, 200L),
        strand = c("+", "+", "+", "+")
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)
    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")

    # rank with ties.method = "first"
    result_ddb <- rank(ddb_gr, ties.method = "first")
    expect_type(result_ddb, "integer")
    expect_equal(length(result_ddb), 4L)
    # Verify ranks are valid (all 1 to n)
    expect_setequal(result_ddb, 1:4)

    # rank with ties.method = "min"
    result_ddb_min <- rank(ddb_gr, ties.method = "min")
    expect_type(result_ddb_min, "integer")
    expect_equal(length(result_ddb_min), 4L)

    unlink(test_tf)
})

test_that("rank works for empty DuckDBGRanges", {
    empty_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0)
    )
    empty_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(empty_df, empty_tf)
    ddb_empty <- DuckDBGRanges(empty_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")

    result <- rank(ddb_empty, ties.method = "first")
    expect_type(result, "integer")
    expect_equal(length(result), 0L)

    unlink(empty_tf)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Range restriction method tests
###

test_that("trim works for DuckDBGRanges", {
    # Create test data with seqinfo
    seqinfo <- Seqinfo(c("chr1", "chr2"), c(200L, 300L), NA, "test")

    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(150L, 1L, 250L),
        end = c(250L, 50L, 400L),  # Some ranges extend beyond seqlengths
        strand = c("+", "+", "+")
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)
    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id", seqinfo = seqinfo)
    gr <- as(ddb_gr, "GRanges")

    # trim
    result_ddb <- trim(ddb_gr)
    result_gr <- trim(gr)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("restrict works for DuckDBGRanges", {
    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(50L, 150L, 100L),
        end = c(150L, 250L, 200L),
        strand = c("+", "+", "+")
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)
    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # restrict with start bound
    result_ddb <- restrict(ddb_gr, start = 100L)
    result_gr <- restrict(gr, start = 100L)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    # restrict with end bound
    result_ddb_end <- restrict(ddb_gr, end = 200L)
    result_gr_end <- restrict(gr, end = 200L)
    checkDuckDBGRanges(result_ddb_end, result_gr_end)

    # restrict with both bounds
    result_ddb_both <- restrict(ddb_gr, start = 100L, end = 200L)
    result_gr_both <- restrict(gr, start = 100L, end = 200L)
    checkDuckDBGRanges(result_ddb_both, result_gr_both)

    # restrict with keep.all.ranges = TRUE
    result_ddb_keep <- restrict(ddb_gr, start = 100L, end = 120L, keep.all.ranges = TRUE)
    result_gr_keep <- restrict(gr, start = 100L, end = 120L, keep.all.ranges = TRUE)
    checkDuckDBGRanges(result_ddb_keep, result_gr_keep)

    unlink(test_tf)
})

test_that("sort works for empty DuckDBGRanges", {
    empty_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0)
    )
    empty_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(empty_df, empty_tf)
    ddb_empty <- DuckDBGRanges(empty_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")

    result <- sort(ddb_empty)
    expect_s4_class(result, "DuckDBGRanges")
    expect_equal(length(result), 0L)

    unlink(empty_tf)
})

test_that("order works for empty DuckDBGRanges", {
    empty_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0)
    )
    empty_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(empty_df, empty_tf)
    ddb_empty <- DuckDBGRanges(empty_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")

    result <- order(ddb_empty)
    expect_type(result, "integer")
    expect_equal(length(result), 0L)

    unlink(empty_tf)
})

test_that("restrict works for empty DuckDBGRanges", {
    empty_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0)
    )
    empty_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(empty_df, empty_tf)
    ddb_empty <- DuckDBGRanges(empty_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")

    result <- restrict(ddb_empty, start = 100L)
    expect_s4_class(result, "DuckDBGRanges")
    expect_equal(length(result), 0L)

    unlink(empty_tf)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reduce() tests
###

test_that("reduce works for DuckDBGRanges with overlapping ranges", {
    # Create overlapping ranges
    test_df <- data.frame(
        id = c("a", "b", "c", "d", "e"),
        seqnames = c("chr1", "chr1", "chr1", "chr2", "chr2"),
        start = c(100L, 120L, 200L, 50L, 100L),
        end = c(150L, 180L, 250L, 80L, 130L),
        strand = c("+", "+", "+", "+", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # Basic reduce
    result_ddb <- reduce(ddb_gr)
    result_gr <- reduce(gr)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("reduce works with ignore.strand=TRUE", {
    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(100L, 120L, 200L),
        end = c(150L, 180L, 250L),
        strand = c("+", "-", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- reduce(ddb_gr, ignore.strand = TRUE)
    result_gr <- reduce(gr, ignore.strand = TRUE)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("reduce works with min.gapwidth parameter", {
    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(100L, 160L, 200L),
        end = c(150L, 190L, 250L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # With min.gapwidth = 10, ranges [100,150] and [160,190] should merge
    # because gap is 9 (160 - 150 - 1 = 9) < 10
    result_ddb <- reduce(ddb_gr, min.gapwidth = 10L)
    result_gr <- reduce(gr, min.gapwidth = 10L)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("reduce works for empty DuckDBGRanges", {
    empty_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0)
    )
    empty_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(empty_df, empty_tf)
    ddb_empty <- DuckDBGRanges(empty_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")

    result <- reduce(ddb_empty)
    expect_s4_class(result, "DuckDBGRanges")
    expect_equal(length(result), 0L)

    unlink(empty_tf)
})

test_that("reduce handles non-overlapping ranges correctly", {
    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 300L, 100L),
        end = c(150L, 350L, 150L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- reduce(ddb_gr)
    result_gr <- reduce(gr)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### gaps() tests
###

test_that("gaps works for DuckDBGRanges", {
    # Note: DuckDBGRanges gaps() only returns gaps for seqnames/strands
    # that have ranges in the input (unlike GRanges which returns gaps
    # for all seqname/strand combinations). Use ignore.strand=TRUE for
    # exact matching behavior.
    seqinfo <- Seqinfo(seqnames = c("chr1", "chr2"),
                       seqlengths = c(1000L, 500L))

    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 300L, 100L),
        end = c(150L, 350L, 200L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id", seqinfo = seqinfo)
    gr <- as(ddb_gr, "GRanges")

    # Test with ignore.strand=TRUE for exact matching
    result_ddb <- gaps(ddb_gr, ignore.strand = TRUE)
    result_gr <- gaps(gr, ignore.strand = TRUE)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("gaps works with ignore.strand=TRUE", {
    seqinfo <- Seqinfo(seqnames = c("chr1"),
                       seqlengths = c(1000L))

    test_df <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(100L, 200L),
        end = c(150L, 250L),
        strand = c("+", "-"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id", seqinfo = seqinfo)
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- gaps(ddb_gr, ignore.strand = TRUE)
    result_gr <- gaps(gr, ignore.strand = TRUE)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("gaps works with custom start and end bounds", {
    # Note: Use ignore.strand=TRUE for exact match (see gaps() documentation)
    seqinfo <- Seqinfo(seqnames = c("chr1"),
                       seqlengths = c(1000L))

    test_df <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(100L, 300L),
        end = c(150L, 350L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id", seqinfo = seqinfo)
    gr <- as(ddb_gr, "GRanges")

    # Custom bounds with ignore.strand=TRUE for exact match
    result_ddb <- gaps(ddb_gr, start = 50L, end = 500L, ignore.strand = TRUE)
    result_gr <- gaps(gr, start = 50L, end = 500L, ignore.strand = TRUE)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("gaps works for empty DuckDBGRanges", {
    seqinfo <- Seqinfo(seqnames = c("chr1"),
                       seqlengths = c(1000L))

    empty_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0)
    )
    empty_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(empty_df, empty_tf)
    ddb_empty <- DuckDBGRanges(empty_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand",
                               seqinfo = seqinfo)

    # For empty ranges, gaps should return the full sequence
    result <- gaps(ddb_empty)
    expect_s4_class(result, "DuckDBGRanges")
    # Should have gaps for each seqname/strand combination with valid seqlengths
    expect_true(length(result) > 0L)

    unlink(empty_tf)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### disjoin() tests
###

test_that("disjoin works for DuckDBGRanges with overlapping ranges", {
    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(100L, 120L, 200L),
        end = c(150L, 180L, 250L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- disjoin(ddb_gr)
    result_gr <- disjoin(gr)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("disjoin works with ignore.strand=TRUE", {
    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(100L, 120L, 200L),
        end = c(150L, 180L, 250L),
        strand = c("+", "-", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- disjoin(ddb_gr, ignore.strand = TRUE)
    result_gr <- disjoin(gr, ignore.strand = TRUE)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("disjoin handles non-overlapping ranges correctly", {
    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr2"),
        start = c(100L, 300L, 100L),
        end = c(150L, 350L, 150L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- disjoin(ddb_gr)
    result_gr <- disjoin(gr)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

test_that("disjoin works for empty DuckDBGRanges", {
    empty_df <- data.frame(
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0)
    )
    empty_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(empty_df, empty_tf)
    ddb_empty <- DuckDBGRanges(empty_tf, seqnames = "seqnames",
                               start = "start", end = "end", strand = "strand")

    result <- disjoin(ddb_empty)
    expect_s4_class(result, "DuckDBGRanges")
    expect_equal(length(result), 0L)

    unlink(empty_tf)
})

test_that("disjoin produces correct breakpoints", {
    # Test case where we have overlapping ranges [1,10], [5,15], [12,20]
    # Expected disjoin result: [1,4], [5,10], [11,11], [12,15], [16,20]
    test_df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(1L, 5L, 12L),
        end = c(10L, 15L, 20L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    test_tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(test_df, test_tf)

    ddb_gr <- DuckDBGRanges(test_tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    result_ddb <- disjoin(ddb_gr)
    result_gr <- disjoin(gr)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(test_tf)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Set operations tests: union(), intersect(), setdiff()
###

test_that("union works for DuckDBGRanges", {
    # Create two sets of ranges with some overlap
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(100L, 300L),
        end = c(200L, 400L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y"),
        seqnames = c("chr1", "chr2"),
        start = c(150L, 100L),
        end = c(250L, 200L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    result_ddb <- union(ddb_gr1, ddb_gr2)
    result_gr <- union(gr1, gr2)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

test_that("union works with ignore.strand=TRUE", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(100L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr1"),
        start = c(150L),
        end = c(250L),
        strand = c("-"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    result_ddb <- union(ddb_gr1, ddb_gr2, ignore.strand = TRUE)
    result_gr <- union(gr1, gr2, ignore.strand = TRUE)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

test_that("union works with DuckDBGRanges and GRanges", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(100L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    tf1 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr2 <- GRanges("chr1:150-250:+")

    result_ddb <- union(ddb_gr1, gr2)
    expect_s4_class(result_ddb, "DuckDBGRanges")

    result_ddb2 <- union(gr2, ddb_gr1)
    expect_s4_class(result_ddb2, "DuckDBGRanges")

    unlink(tf1)
})

test_that("intersect works for DuckDBGRanges", {
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(100L, 300L),
        end = c(200L, 400L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y"),
        seqnames = c("chr1", "chr1"),
        start = c(150L, 350L),
        end = c(250L, 450L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    result_ddb <- intersect(ddb_gr1, ddb_gr2)
    result_gr <- intersect(gr1, gr2)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

test_that("intersect works with ignore.strand=TRUE", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(100L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr1"),
        start = c(150L),
        end = c(250L),
        strand = c("-"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    result_ddb <- intersect(ddb_gr1, ddb_gr2, ignore.strand = TRUE)
    result_gr <- intersect(gr1, gr2, ignore.strand = TRUE)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

test_that("intersect returns empty when no overlap", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(100L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr2"),
        start = c(100L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")

    result_ddb <- intersect(ddb_gr1, ddb_gr2)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    expect_equal(length(result_ddb), 0L)

    unlink(c(tf1, tf2))
})

test_that("setdiff works for DuckDBGRanges", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(100L),
        end = c(300L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr1"),
        start = c(150L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    result_ddb <- setdiff(ddb_gr1, ddb_gr2)
    result_gr <- setdiff(gr1, gr2)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

test_that("setdiff works with ignore.strand=TRUE", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(100L),
        end = c(300L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr1"),
        start = c(150L),
        end = c(200L),
        strand = c("-"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    result_ddb <- setdiff(ddb_gr1, ddb_gr2, ignore.strand = TRUE)
    result_gr <- setdiff(gr1, gr2, ignore.strand = TRUE)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

test_that("setdiff returns x when y is empty or no overlap", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(100L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr2"),
        start = c(100L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    result_ddb <- setdiff(ddb_gr1, ddb_gr2)
    result_gr <- setdiff(gr1, gr2)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    checkDuckDBGRanges(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

test_that("setdiff with y fully covering x returns empty", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(150L),
        end = c(200L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr1"),
        start = c(100L),
        end = c(300L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )

    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)

    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")

    result_ddb <- setdiff(ddb_gr1, ddb_gr2)
    result_gr <- setdiff(gr1, gr2)
    expect_s4_class(result_ddb, "DuckDBGRanges")
    expect_equal(length(result_ddb), length(result_gr))

    unlink(c(tf1, tf2))
})

### =========================================================================
### Phase 4: Nearest neighbor methods (precede, follow, nearest, distanceToNearest)
### =========================================================================

test_that("precede works correctly", {
    # Query ranges
    df1 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(10L, 50L, 100L),
        end = c(20L, 60L, 110L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    # Subject ranges
    df2 <- data.frame(
        id = c("x", "y", "z"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(30L, 70L, 200L),
        end = c(40L, 80L, 210L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- precede(ddb_gr1, ddb_gr2)
    result_gr <- precede(gr1, gr2)
    expect_equal(result_ddb, result_gr)
    
    unlink(c(tf1, tf2))
})

test_that("precede with ignore.strand=TRUE works", {
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(10L, 50L),
        end = c(20L, 60L),
        strand = c("+", "-"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y", "z"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(30L, 40L, 70L),
        end = c(35L, 45L, 80L),
        strand = c("+", "+", "-"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- precede(ddb_gr1, ddb_gr2, ignore.strand = TRUE)
    result_gr <- precede(gr1, gr2, ignore.strand = TRUE)
    expect_equal(result_ddb, result_gr)
    
    unlink(c(tf1, tf2))
})

test_that("precede with empty inputs works", {
    df1 <- data.frame(
        id = character(0),
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr1"),
        start = c(30L),
        end = c(40L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    
    result <- precede(ddb_gr1, ddb_gr2)
    expect_equal(result, integer(0))
    
    unlink(c(tf1, tf2))
})

test_that("follow works correctly", {
    df1 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(10L, 50L, 100L),
        end = c(20L, 60L, 110L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y", "z"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(30L, 70L, 200L),
        end = c(40L, 80L, 210L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- follow(ddb_gr1, ddb_gr2)
    result_gr <- follow(gr1, gr2)
    expect_equal(result_ddb, result_gr)
    
    unlink(c(tf1, tf2))
})

test_that("follow with ignore.strand=TRUE works", {
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(50L, 100L),
        end = c(60L, 110L),
        strand = c("+", "-"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y", "z"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(10L, 20L, 70L),
        end = c(15L, 30L, 80L),
        strand = c("+", "+", "-"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- follow(ddb_gr1, ddb_gr2, ignore.strand = TRUE)
    result_gr <- follow(gr1, gr2, ignore.strand = TRUE)
    expect_equal(result_ddb, result_gr)
    
    unlink(c(tf1, tf2))
})

test_that("nearest works correctly", {
    df1 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(10L, 50L, 100L),
        end = c(20L, 60L, 110L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y", "z"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(30L, 70L, 200L),
        end = c(40L, 80L, 210L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- nearest(ddb_gr1, ddb_gr2)
    result_gr <- nearest(gr1, gr2)
    expect_equal(result_ddb, result_gr)

    unlink(c(tf1, tf2))
})

test_that("nearest(x)/distanceToNearest(x) with no subject exclude self-hits", {
    # Base GenomicRanges nearest(x, missing) uses drop.self=TRUE: a range is never
    # its own nearest neighbour. The self-query methods must match the base oracle
    # (c(2, 1, 1) here, not the self-hits c(1, 2, 3)).
    df <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(10L, 100L, 500L),
        end = c(20L, 110L, 510L),
        strand = c("*", "*", "*"),
        stringsAsFactors = FALSE
    )
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    ddb_gr <- DuckDBGRanges(tf, seqnames = "seqnames", start = "start",
                            end = "end", strand = "strand", keycol = "id")
    gr <- as(ddb_gr, "GRanges")

    # No subject -> self excluded, matching base.
    expect_equal(nearest(ddb_gr), nearest(gr))

    # distanceToNearest: compare by components (base returns a SortedByQueryHits
    # subclass, the port a plain Hits — identical content, different label).
    d2n_ddb <- distanceToNearest(ddb_gr)
    d2n_gr <- distanceToNearest(gr)
    expect_equal(queryHits(d2n_ddb), queryHits(d2n_gr))
    expect_equal(subjectHits(d2n_ddb), subjectHits(d2n_gr))
    expect_equal(mcols(d2n_ddb)$distance, mcols(d2n_gr)$distance)

    # An explicit subject == x keeps self-hits (base does not drop self there).
    expect_equal(nearest(ddb_gr, ddb_gr), nearest(gr, gr))

    unlink(tf)
})

test_that("nearest with overlapping ranges returns correct index", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(50L),
        end = c(70L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y"),
        seqnames = c("chr1", "chr1"),
        start = c(60L, 100L),
        end = c(80L, 120L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- nearest(ddb_gr1, ddb_gr2)
    result_gr <- nearest(gr1, gr2)
    # Both should return 1 (overlapping has distance 0)
    expect_equal(result_ddb, result_gr)
    
    unlink(c(tf1, tf2))
})

test_that("nearest with ignore.strand=TRUE works", {
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(10L, 50L),
        end = c(20L, 60L),
        strand = c("+", "-"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y", "z"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(30L, 40L, 70L),
        end = c(35L, 45L, 80L),
        strand = c("+", "+", "-"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- nearest(ddb_gr1, ddb_gr2, ignore.strand = TRUE)
    result_gr <- nearest(gr1, gr2, ignore.strand = TRUE)
    expect_equal(result_ddb, result_gr)
    
    unlink(c(tf1, tf2))
})

test_that("distanceToNearest works correctly", {
    df1 <- data.frame(
        id = c("a", "b", "c"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(10L, 50L, 100L),
        end = c(20L, 60L, 110L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y", "z"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(30L, 70L, 200L),
        end = c(40L, 80L, 210L),
        strand = c("+", "+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- distanceToNearest(ddb_gr1, ddb_gr2)
    result_gr <- distanceToNearest(gr1, gr2)
    
    expect_s4_class(result_ddb, "Hits")
    expect_equal(queryHits(result_ddb), queryHits(result_gr))
    expect_equal(subjectHits(result_ddb), subjectHits(result_gr))
    expect_equal(mcols(result_ddb)$distance, mcols(result_gr)$distance)
    
    unlink(c(tf1, tf2))
})

test_that("distanceToNearest with overlapping ranges returns distance 0", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(50L),
        end = c(70L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y"),
        seqnames = c("chr1", "chr1"),
        start = c(60L, 100L),
        end = c(80L, 120L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- distanceToNearest(ddb_gr1, ddb_gr2)
    result_gr <- distanceToNearest(gr1, gr2)
    
    # Should return distance of 0 for overlapping range
    expect_equal(mcols(result_ddb)$distance, mcols(result_gr)$distance)
    expect_equal(subjectHits(result_ddb), subjectHits(result_gr))
    
    unlink(c(tf1, tf2))
})

test_that("distanceToNearest with empty input returns empty Hits", {
    df1 <- data.frame(
        id = character(0),
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr1"),
        start = c(30L),
        end = c(40L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    
    result <- distanceToNearest(ddb_gr1, ddb_gr2)
    expect_s4_class(result, "Hits")
    expect_equal(length(result), 0)
    
    unlink(c(tf1, tf2))
})

test_that("distanceToNearest with ignore.strand=TRUE works", {
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(10L, 50L),
        end = c(20L, 60L),
        strand = c("+", "-"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y", "z"),
        seqnames = c("chr1", "chr1", "chr1"),
        start = c(30L, 40L, 70L),
        end = c(35L, 45L, 80L),
        strand = c("+", "+", "-"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- distanceToNearest(ddb_gr1, ddb_gr2, ignore.strand = TRUE)
    result_gr <- distanceToNearest(gr1, gr2, ignore.strand = TRUE)
    
    expect_equal(queryHits(result_ddb), queryHits(result_gr))
    expect_equal(subjectHits(result_ddb), subjectHits(result_gr))
    expect_equal(mcols(result_ddb)$distance, mcols(result_gr)$distance)
    
    unlink(c(tf1, tf2))
})

test_that("precede/follow/nearest with mixed types work", {
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(10L, 50L),
        end = c(20L, 60L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y"),
        seqnames = c("chr1", "chr1"),
        start = c(30L, 70L),
        end = c(40L, 80L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    # Test DuckDBGRanges, GRanges
    expect_equal(precede(ddb_gr1, gr2), precede(gr1, gr2))
    expect_equal(follow(ddb_gr1, gr2), follow(gr1, gr2))
    expect_equal(nearest(ddb_gr1, gr2), nearest(gr1, gr2))
    
    # Test GRanges, DuckDBGRanges
    expect_equal(precede(gr1, ddb_gr2), precede(gr1, gr2))
    expect_equal(follow(gr1, ddb_gr2), follow(gr1, gr2))
    expect_equal(nearest(gr1, ddb_gr2), nearest(gr1, gr2))

    unlink(c(tf1, tf2))
})

test_that("precede/follow/nearest are strand-directional and '*'-compatible", {
    # A '-' strand query inverts precede/follow direction, and '*' is compatible
    # with any strand. Both used to diverge from base GenomicRanges (a strict
    # strand equi-join dropped every '*' pair, and precede/follow used a fixed
    # genomic direction regardless of strand). Oracle-check against base across
    # '+', '-' and '*' queries with mixed-strand subjects.
    df_q <- data.frame(
        id = c("qp", "qm", "qs"),
        seqnames = "chr1",
        start = c(100L, 100L, 100L),
        end = c(110L, 110L, 110L),
        strand = c("+", "-", "*"),
        stringsAsFactors = FALSE
    )
    df_s <- data.frame(
        id = c("lp", "rp", "lm", "rm"),
        seqnames = "chr1",
        start = c(10L, 200L, 20L, 300L),
        end = c(50L, 300L, 60L, 400L),
        strand = c("+", "+", "-", "-"),
        stringsAsFactors = FALSE
    )
    tfq <- tempfile(fileext = ".parquet")
    tfs <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df_q, tfq)
    arrow::write_parquet(df_s, tfs)
    q <- DuckDBGRanges(tfq, seqnames = "seqnames", start = "start", end = "end",
                       strand = "strand", keycol = "id")
    s <- DuckDBGRanges(tfs, seqnames = "seqnames", start = "start", end = "end",
                       strand = "strand", keycol = "id")
    gq <- as(q, "GRanges")
    gs <- as(s, "GRanges")

    expect_equal(precede(q, s), precede(gq, gs))
    expect_equal(follow(q, s), follow(gq, gs))
    expect_equal(nearest(q, s), nearest(gq, gs))

    unlink(c(tfq, tfs))
})

test_that("follow/precede tie-breaks match base select last/first", {
    # Identical-coordinate subjects: base precede uses select="first" (smallest
    # index), follow uses select="last" (largest index).
    df_q <- data.frame(id = "q", seqnames = "chr1", start = 100L, end = 110L,
                       strand = "+", stringsAsFactors = FALSE)
    df_left <- data.frame(id = c("a", "b"), seqnames = "chr1",
                          start = c(50L, 50L), end = c(90L, 90L),
                          strand = c("+", "+"), stringsAsFactors = FALSE)
    df_right <- data.frame(id = c("a", "b"), seqnames = "chr1",
                           start = c(200L, 200L), end = c(250L, 250L),
                           strand = c("+", "+"), stringsAsFactors = FALSE)
    tfq <- tempfile(fileext = ".parquet")
    tfl <- tempfile(fileext = ".parquet")
    tfr <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df_q, tfq)
    arrow::write_parquet(df_left, tfl)
    arrow::write_parquet(df_right, tfr)
    q <- DuckDBGRanges(tfq, seqnames = "seqnames", start = "start", end = "end",
                       strand = "strand", keycol = "id")
    sl <- DuckDBGRanges(tfl, seqnames = "seqnames", start = "start", end = "end",
                        strand = "strand", keycol = "id")
    sr <- DuckDBGRanges(tfr, seqnames = "seqnames", start = "start", end = "end",
                        strand = "strand", keycol = "id")

    # follow -> largest index (2); precede -> smallest index (1).
    expect_equal(follow(q, sl), follow(as(q, "GRanges"), as(sl, "GRanges")))
    expect_equal(precede(q, sr), precede(as(q, "GRanges"), as(sr, "GRanges")))

    unlink(c(tfq, tfl, tfr))
})

### =========================================================================
### Phase 5: Tiling methods (tile, slidingWindows, pgap)
### =========================================================================

test_that("tile with n works correctly", {
    df <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(1L, 101L),
        end = c(50L, 200L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    
    ddb_gr <- DuckDBGRanges(tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")
    
    result_ddb <- tile(ddb_gr, n = 5)
    result_gr <- tile(gr, n = 5)
    
    expect_s4_class(result_ddb, "GRangesList")
    expect_equal(length(result_ddb), length(result_gr))
    
    # Check each element has the same number of tiles
    for (i in seq_along(result_ddb)) {
        expect_equal(length(result_ddb[[i]]), length(result_gr[[i]]))
    }
    
    unlink(tf)
})

test_that("tile with width works correctly", {
    df <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(1L, 101L),
        end = c(50L, 200L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    
    ddb_gr <- DuckDBGRanges(tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")
    
    result_ddb <- tile(ddb_gr, width = 10)
    result_gr <- tile(gr, width = 10)
    
    expect_s4_class(result_ddb, "GRangesList")
    expect_equal(length(result_ddb), length(result_gr))
    
    # Check first element
    expect_equal(length(result_ddb[[1]]), length(result_gr[[1]]))
    
    unlink(tf)
})

test_that("tile with empty input returns empty GRangesList", {
    df <- data.frame(
        id = character(0),
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0),
        stringsAsFactors = FALSE
    )
    
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    
    ddb_gr <- DuckDBGRanges(tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    
    result <- tile(ddb_gr, n = 5)
    expect_s4_class(result, "GRangesList")
    expect_equal(length(result), 0)
    
    unlink(tf)
})

test_that("slidingWindows works correctly", {
    df <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(1L),
        end = c(100L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    
    ddb_gr <- DuckDBGRanges(tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")
    
    result_ddb <- slidingWindows(ddb_gr, width = 20, step = 10)
    result_gr <- slidingWindows(gr, width = 20, step = 10)
    
    expect_s4_class(result_ddb, "GRangesList")
    expect_equal(length(result_ddb), length(result_gr))
    expect_equal(length(result_ddb[[1]]), length(result_gr[[1]]))
    
    unlink(tf)
})

test_that("slidingWindows with step=1 works", {
    df <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(1L),
        end = c(20L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    
    ddb_gr <- DuckDBGRanges(tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    gr <- as(ddb_gr, "GRanges")
    
    result_ddb <- slidingWindows(ddb_gr, width = 5, step = 1)
    result_gr <- slidingWindows(gr, width = 5, step = 1)
    
    expect_s4_class(result_ddb, "GRangesList")
    expect_equal(length(result_ddb[[1]]), length(result_gr[[1]]))
    
    unlink(tf)
})

test_that("slidingWindows with empty input returns empty GRangesList", {
    df <- data.frame(
        id = character(0),
        seqnames = character(0),
        start = integer(0),
        end = integer(0),
        strand = character(0),
        stringsAsFactors = FALSE
    )
    
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)
    
    ddb_gr <- DuckDBGRanges(tf, seqnames = "seqnames",
                            start = "start", end = "end", strand = "strand",
                            keycol = "id")
    
    result <- slidingWindows(ddb_gr, width = 10, step = 5)
    expect_s4_class(result, "GRangesList")
    expect_equal(length(result), 0)
    
    unlink(tf)
})

test_that("pgap works correctly", {
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(10L, 100L),
        end = c(50L, 150L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y"),
        seqnames = c("chr1", "chr1"),
        start = c(60L, 160L),
        end = c(90L, 200L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- pgap(ddb_gr1, ddb_gr2)
    result_gr <- pgap(gr1, gr2)
    
    expect_s4_class(result_ddb, "GRanges")
    expect_equal(length(result_ddb), length(result_gr))
    
    unlink(c(tf1, tf2))
})

test_that("pgap with overlapping ranges works", {
    df1 <- data.frame(
        id = c("a"),
        seqnames = c("chr1"),
        start = c(10L),
        end = c(50L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x"),
        seqnames = c("chr1"),
        start = c(40L),
        end = c(80L),
        strand = c("+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    result_ddb <- pgap(ddb_gr1, ddb_gr2)
    result_gr <- pgap(gr1, gr2)
    
    expect_s4_class(result_ddb, "GRanges")
    # For overlapping ranges, the gap width should be 0
    expect_equal(width(result_ddb), width(result_gr))
    
    unlink(c(tf1, tf2))
})

test_that("pgap with mixed types works", {
    df1 <- data.frame(
        id = c("a", "b"),
        seqnames = c("chr1", "chr1"),
        start = c(10L, 100L),
        end = c(50L, 150L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    df2 <- data.frame(
        id = c("x", "y"),
        seqnames = c("chr1", "chr1"),
        start = c(60L, 160L),
        end = c(90L, 200L),
        strand = c("+", "+"),
        stringsAsFactors = FALSE
    )
    
    tf1 <- tempfile(fileext = ".parquet")
    tf2 <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df1, tf1)
    arrow::write_parquet(df2, tf2)
    
    ddb_gr1 <- DuckDBGRanges(tf1, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    ddb_gr2 <- DuckDBGRanges(tf2, seqnames = "seqnames",
                             start = "start", end = "end", strand = "strand",
                             keycol = "id")
    gr1 <- as(ddb_gr1, "GRanges")
    gr2 <- as(ddb_gr2, "GRanges")
    
    # Test DuckDBGRanges, GRanges
    result_mixed1 <- pgap(ddb_gr1, gr2)
    expect_s4_class(result_mixed1, "GRanges")
    
    # Test GRanges, DuckDBGRanges
    result_mixed2 <- pgap(gr1, ddb_gr2)
    expect_s4_class(result_mixed2, "GRanges")
    
    unlink(c(tf1, tf2))
})
