#' Utilities for DuckDBGRangesList objects
#'
#' @description
#' Various utility methods for DuckDBGRangesList objects including inter-range
#' operations. These methods translate R operations into efficient SQL queries
#' executed directly in DuckDB.
#'
#' @details
#' The DuckDBGRangesList methods translate genomic range operations into SQL
#' queries, enabling efficient manipulation of large disk-backed genomic
#' annotations without loading them into memory.
#'
#' @section Inter-range Methods:
#' Inter-range methods operate across multiple ranges within each list element:
#' \describe{
#'   \item{\code{range(x, ..., with.revmap=FALSE, ignore.strand=FALSE, na.rm=FALSE)}:}{
#'     Returns the range (min start to max end) per list element and
#'     seqname/strand combination. Computed via SQL MIN/MAX aggregation
#'     with grouping by the partitioning column.
#'     Returns a DuckDBGRangesList object.
#'   }
#' }
#'
#' @param x A DuckDBGRangesList object.
#' @param with.revmap If \code{TRUE}, include reverse mapping in output
#'   (not supported, will warn).
#' @param ignore.strand If \code{TRUE}, strand information is ignored when
#'   computing operations. Default is \code{FALSE}.
#' @param na.rm Ignored.
#' @param ... Additional arguments passed to methods.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{DuckDBGRangesList-class}} for the main class
#'   \item \code{\link{DuckDBGRanges-utils}} for DuckDBGRanges utilities
#'   \item \code{\link[GenomicRanges]{GRangesList-class}} for the base class
#' }
#'
#' @examples
#' # Create example DuckDBGRangesList
#' df <- data.frame(
#'     group = c("A", "B"),
#'     seqnames = I(list(
#'         rep(c("chr1", "chr2"), c(2, 3)),
#'         rep(c("chr1", "chr2"), c(3, 2))
#'     )),
#'     start = I(list(
#'         c(100L, 200L, 150L, 300L, 400L),
#'         c(50L, 100L, 200L, 250L, 300L)
#'     )),
#'     end = I(list(
#'         c(150L, 250L, 180L, 350L, 450L),
#'         c(80L, 150L, 250L, 300L, 350L)
#'     )),
#'     strand = I(list(
#'         rep(c("+", "-"), c(3, 2)),
#'         rep(c("+", "-"), c(2, 3))
#'     ))
#' )
#' tf <- tempfile(fileext = ".parquet")
#' arrow::write_parquet(df, tf)
#' grlist <- DuckDBGRangesList(tf, seqnames = "seqnames", start = "start", 
#'                             end = "end", strand = "strand",
#'                             keycol = list(group = c("A", "B")))
#'
#' # Get range per list element
#' range(grlist)
#'
#' @return
#' \code{range()} returns a DuckDBGRangesList giving the range (minimum start to
#' maximum end) spanned within each list element, per seqname/strand.
#'
#' @aliases
#' range,DuckDBGRangesList-method
#'
#' @include DuckDBGRangesList-class.R
#'
#' @keywords utilities methods
#'
#' @name DuckDBGRangesList-utils
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Inter-range methods
###

#' @export
#' @importFrom DBI dbQuoteIdentifier dbQuoteLiteral
#' @importFrom dplyr collect count group_by mutate sql summarize
#' @importFrom DuckDBDataFrame DuckDBDataFrame tblconn
#' @importFrom Seqinfo seqinfo
#' @importFrom S4Vectors make_zero_col_DFrame new2
setMethod("range", "DuckDBGRangesList",
function(x, ..., with.revmap = FALSE, ignore.strand = FALSE, na.rm = FALSE)
{
    if (!identical(na.rm, FALSE))
        warning("'na.rm' argument is ignored")

    if (with.revmap)
        warning("'with.revmap' is not supported for DuckDBGRangesList, ignoring")

    args <- list(x, ...)
    if (length(args) > 1L) {
        stop("multiple arguments not yet supported for DuckDBGRangesList range()")
    }

    n <- length(x)
    if (n == 0L)
        return(x)

    # Extract metadata
    keycol <- names(x@frame@keycols)[1L]
    keyvals <- x@frame@keycols[[1L]]
    conn <- tblconn(x@frame)

    # Step 1: UNNEST LIST columns to flat rows
    unnest_expr <- list(seqnames = sql("UNNEST(seqnames)"),
                        start = sql("UNNEST(start)"),
                        width = sql("UNNEST(width)"),
                        strand = sql("UNNEST(strand)"))
    conn <- mutate(conn, !!!unnest_expr)

    # Step 2: Group and compute range per element/seqnames/strand
    # Compute end from start + width before grouping
    if (ignore.strand) {
        groups <- setNames(
            list(as.name(keycol), as.name("seqnames")),
            c(keycol, "seqnames")
        )
        aggr <- list(
            start = call("min", as.name("start"), na.rm = TRUE),
            end = call("max", call("+", call("-", call("+", as.name("start"), as.name("width")), 1L), 0L), na.rm = TRUE)
        )
        range_conn <- conn |>
            group_by(!!!groups) |>
            summarize(!!!aggr, .groups = "drop")
        width_mutate <- list(
            strand = "*",
            width = call("+", call("-", as.name("end"), as.name("start")), 1L)
        )
        range_conn <- mutate(range_conn, !!!width_mutate)
    } else {
        groups <- setNames(
            list(as.name(keycol), as.name("seqnames"), as.name("strand")),
            c(keycol, "seqnames", "strand")
        )
        aggr <- list(
            start = call("min", as.name("start"), na.rm = TRUE),
            end = call("max", call("+", call("-", call("+", as.name("start"), as.name("width")), 1L), 0L), na.rm = TRUE)
        )
        range_conn <- conn |>
            group_by(!!!groups) |>
            summarize(!!!aggr, .groups = "drop")
        width_mutate <- list(
            width = call("+", call("-", as.name("end"), as.name("start")), 1L)
        )
        range_conn <- mutate(range_conn, !!!width_mutate)
    }

    # Step 3: Re-aggregate into LIST[] columns grouped by keycol
    list_groups <- setNames(list(as.name(keycol)), keycol)
    list_aggr <- list(seqnames = sql("list(seqnames)"),
                      start = sql("list(start)"),
                      width = sql("list(width)"),
                      strand = sql("list(strand)"))
    list_conn <- range_conn |>
        group_by(!!!list_groups) |>
        summarize(!!!list_aggr, .groups = "drop")

    # Step 4: Add integer index column using CASE WHEN to preserve order
    # Build CASE WHEN expression that maps keyvals to their positions
    conn_db <- dbconn(x)
    quoted_keycol <- as.character(dbQuoteIdentifier(conn_db, keycol))
    case_expr <- sprintf("CASE %s %s END",
                         paste0(sprintf("WHEN %s = %s THEN %d",
                                        quoted_keycol,
                                        dbQuoteLiteral(conn_db, keyvals),
                                        seq_along(keyvals)),
                                collapse = " "),
                         sprintf("ELSE %d", length(keyvals) + 1L))

    idx_mutate <- list(.row_idx = sql(case_expr))
    list_conn <- mutate(list_conn, !!!idx_mutate)

    # Step 5: Build datacols with computed end column
    datacols <- as.expression(c(seqnames = as.name("seqnames"),
                                start = as.name("start"),
                                width = as.name("width"),
                                strand = as.name("strand")))
    datacols[["end"]] <- call("list_transform",
                              call("list_zip", as.name("start"), as.name("width")),
                              sql("x -> x[1] + x[2] - 1"))

    # Create new DuckDBDataFrame using .row_idx as keycol
    # The integer indices 1, 2, 3, ... will be naturally sorted
    new_frame <- DuckDBDataFrame(list_conn, datacols = datacols, keycol = ".row_idx")

    # Map row indices back to original names
    idx_vals <- new_frame@keycols[[".row_idx"]]
    new_frame@keycols[[".row_idx"]] <- setNames(idx_vals, keyvals[idx_vals])

    # Return new DuckDBGRangesList
    new2("DuckDBGRangesList",
         frame = new_frame,
         seqinfo = x@seqinfo,
         elementMetadata = make_zero_col_DFrame(length(x)),
         check = FALSE)
})
