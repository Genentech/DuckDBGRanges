# Tests the basic functions of a DuckDBGRangesList.
# library(testthat); library(DuckDBGRanges); source("setup.R"); source("test-DuckDBGRangesList.R")

library(GenomicRanges)

test_that("basic methods work as expected for a DuckDBGRangesList", {
    object <- DuckDBGRangesList(grlist_tf, seqnames = "seqnames", start = "start", end = "end", strand = "strand",
                                mcols = c("score", "GC", "label", "description"),
                                keycol = list(id = c("gr1", "gr2", "gr3", "gr4")))

    expected <- grlist_grl

    checkDuckDBGRangesList(object, expected)
})

test_that("renaming list elements work for a DuckDBGRangesList", {
    object <- DuckDBGRangesList(grlist_tf, seqnames = "seqnames", start = "start", end = "end", strand = "strand",
                                mcols = c("score", "GC", "label", "description"),
                                keycol = list(id = c("gr1", "gr2", "gr3", "gr4")))
    names(object) <- head(letters, length(object))

    expected <- grlist_grl
    names(expected) <- head(letters, length(expected))

    checkDuckDBGRangesList(object, expected)
})

test_that("subscripting works for a DuckDBGRangesList", {
    object <- DuckDBGRangesList(grlist_tf, seqnames = "seqnames", start = "start", end = "end", strand = "strand",
                                mcols = c("score", "GC", "label", "description"),
                                keycol = list(id = c("gr1", "gr2", "gr3", "gr4")))

    expected <- grlist_grl

    checkDuckDBGRangesList(object[c(3, 1)], expected[c(3, 1)])

    checkDuckDBGRanges(object[[2]], expected[[2]])
})

test_that("head works for a DuckDBGRangesList", {
    object <- DuckDBGRangesList(grlist_tf, seqnames = "seqnames", start = "start", end = "end", strand = "strand",
                                mcols = c("score", "GC", "label", "description"),
                                keycol = list(id = c("gr1", "gr2", "gr3", "gr4")))

    expected <- grlist_grl

    checkDuckDBGRangesList(head(object, 0), head(expected, 0))
    checkDuckDBGRangesList(head(object, 2), head(expected, 2))
})

test_that("tail works for a DuckDBGRangesList", {
    object <- DuckDBGRangesList(grlist_tf, seqnames = "seqnames", start = "start", end = "end", strand = "strand",
                                mcols = c("score", "GC", "label", "description"),
                                keycol = list(id = c("gr1", "gr2", "gr3", "gr4")))

    expected <- grlist_grl

    checkDuckDBGRangesList(tail(object, 0), tail(expected, 0))
    checkDuckDBGRangesList(tail(object, 2), tail(expected, 2))
})

test_that("coersion to a GRangesList works for a DuckDBGRangesList", {
    seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")

    object <- DuckDBGRangesList(grlist_tf, seqnames = "seqnames", start = "start", end = "end", strand = "strand",
                                mcols = c("score", "GC", "label", "description"),
                                keycol = list(id = c("gr1", "gr2", "gr3", "gr4")),
                                seqinfo = seqinfo)
    object <- as(object, "CompressedGRangesList")

    expected <- grlist_grl
    expected@unlistData@seqinfo <- seqinfo

    for (i in names(object)) {
        object_i <- object[[i]]
        expected_i <- expected[[i]]
        expect_identical(object_i, expected_i)
    }
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Inter-range method tests
###

test_that("range works for DuckDBGRangesList", {
    object <- DuckDBGRangesList(grlist_tf, seqnames = "seqnames", start = "start", end = "end", strand = "strand",
                                mcols = c("score", "GC", "label", "description"),
                                keycol = list(id = c("gr1", "gr2", "gr3", "gr4")))

    expected <- grlist_grl

    # range with strand
    result_ddb <- range(object)
    result_gr <- range(expected)
    expect_s4_class(result_ddb, "DuckDBGRangesList")
    expect_equal(length(result_ddb), length(result_gr))
    expect_equal(names(result_ddb), names(result_gr))
    expect_identical(elementNROWS(result_ddb), elementNROWS(result_gr))

    # Check each list element has the same ranges (order may differ)
    for (i in names(result_ddb)) {
        expect_setequal(as.vector(start(result_ddb[[i]])), start(result_gr[[i]]))
        expect_setequal(as.vector(end(result_ddb[[i]])), end(result_gr[[i]]))
    }

    # range ignoring strand
    result_ddb_ign <- range(object, ignore.strand = TRUE)
    result_gr_ign <- range(expected, ignore.strand = TRUE)
    expect_equal(length(result_ddb_ign), length(result_gr_ign))
    expect_equal(names(result_ddb_ign), names(result_gr_ign))
    expect_identical(elementNROWS(result_ddb_ign), elementNROWS(result_gr_ign))
})

test_that("range works for DuckDBGRangesList with different seqnames per group", {
    # Create data with different chromosomes per group
    # Group A: chr1(2), chr2(3) | Group B: chr1(1), chr3(4)
    df <- data.frame(
        group = c("A", "B"),
        seqnames = I(list(
            c("chr1", "chr1", "chr2", "chr2", "chr2"),
            c("chr1", "chr3", "chr3", "chr3", "chr3")
        )),
        start = I(list(
            c(100L, 200L, 150L, 300L, 400L),
            c(50L, 100L, 200L, 250L, 300L)
        )),
        end = I(list(
            c(150L, 250L, 180L, 350L, 450L),
            c(80L, 150L, 250L, 300L, 350L)
        )),
        strand = I(list(
            c("+", "-", "+", "-", "+"),
            c("-", "+", "-", "+", "-")
        ))
    )
    tf <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, tf)

    ddb_grlist <- DuckDBGRangesList(tf, seqnames = "seqnames", start = "start", 
                                   end = "end", strand = "strand", 
                                   keycol = list(group = c("A", "B")))

    # Create expected result using standard GRangesList
    grlist <- GRangesList(
        A = GRanges(
            seqnames = c("chr1", "chr1", "chr2", "chr2", "chr2"),
            ranges = IRanges(start = c(100L, 200L, 150L, 300L, 400L),
                           end = c(150L, 250L, 180L, 350L, 450L)),
            strand = c("+", "-", "+", "-", "+")
        ),
        B = GRanges(
            seqnames = c("chr1", "chr3", "chr3", "chr3", "chr3"),
            ranges = IRanges(start = c(50L, 100L, 200L, 250L, 300L),
                           end = c(80L, 150L, 250L, 300L, 350L)),
            strand = c("-", "+", "-", "+", "-")
        )
    )

    result_ddb <- range(ddb_grlist)
    result_gr <- range(grlist)

    expect_s4_class(result_ddb, "DuckDBGRangesList")
    expect_equal(names(result_ddb), names(result_gr))
    expect_identical(elementNROWS(result_ddb), elementNROWS(result_gr))

    # For each group, check that ranges are correct (order may differ)
    for (i in names(result_ddb)) {
        expect_setequal(as.vector(start(result_ddb[[i]])), start(result_gr[[i]]))
        expect_setequal(as.vector(end(result_ddb[[i]])), end(result_gr[[i]]))
    }

    unlink(tf)
})

test_that("range works for empty DuckDBGRangesList", {
    object <- DuckDBGRangesList(grlist_tf, seqnames = "seqnames", start = "start", 
                                end = "end", strand = "strand",
                                mcols = c("score", "GC", "label", "description"),
                                keycol = list(id = c("gr1", "gr2", "gr3", "gr4")))

    # Get first 0 elements
    empty_list <- head(object, 0)

    result <- range(empty_list)
    expect_s4_class(result, "DuckDBGRangesList")
    expect_equal(length(result), 0L)
})
