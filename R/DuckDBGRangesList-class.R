#' DuckDBGRangesList objects
#'
#' @description
#' The DuckDBGRangesList class extends \linkS4class{GRangesList} to represent
#' a DuckDB table with LIST[] columns as a \linkS4class{GRangesList} object.
#'
#' @details
#' The DuckDBGRangesList class extends the \linkS4class{GRangesList} instead of
#' \linkS4class{GenomicRangesList} class because the \code{rowRanges} slot
#' accepts either a \linkS4class{GenomicRanges} object or a
#' \linkS4class{GRangesList} object.
#'
#' @section Constructor:
#' \describe{
#'   \item{\code{DuckDBGRangesList(conn, seqnames, start, width, strand = NULL,
#'     keycol = NULL, dimtbl = NULL, mcols = NULL, seqinfo = NULL,
#'     seqlengths = NULL)}:}{
#'     Creates a DuckDBGRangesList object from a data source with LIST[] columns.
#'     \describe{
#'       \item{\code{conn}}{
#'         Either a character vector containing the paths to parquet, csv, or
#'         gzipped csv data files; a string that defines a duckdb \code{read_*}
#'         data source; a DuckDBDataFrame object; or a tbl_duckdb_connection
#'         object.
#'       }
#'       \item{\code{seqnames}}{
#'         A string specifying the column from \code{conn} that contains the
#'         LIST of sequence names for each element.
#'       }
#'       \item{\code{start}}{
#'         A string specifying the column from \code{conn} that contains the
#'         LIST of start positions for each element.
#'       }
#'       \item{\code{width}}{
#'         A string specifying the column from \code{conn} that contains the
#'         LIST of widths for each element.
#'       }
#'       \item{\code{strand}}{
#'         Either \code{NULL} or a string specifying the column from \code{conn}
#'         that contains the LIST of strand values for each element.
#'       }
#'       \item{\code{keycol}}{
#'         An optional string specifying the column name from \code{conn} that
#'         will define the foreign key in the underlying table, or a named list
#'         containing a character vector where the name of the list element
#'         defines the foreign key and the character vector set the distinct
#'         values for that key. If missing, a \code{row_number} column is
#'         created as an identifier.
#'       }
#'       \item{\code{dimtbl}}{
#'         An optional named \code{DataFrameList} that specifies the dimension
#'         table associated with the \code{keycol}.
#'       }
#'       \item{\code{mcols}}{
#'         Optional character vector specifying the columns that define the
#'         element-level metadata columns.
#'       }
#'       \item{\code{seqinfo}}{
#'         Either \code{NULL}, or a \code{\link{Seqinfo}} object, or a character
#'         vector of unique sequence names (a.k.a. seqlevels), or a named
#'         numeric vector of sequence lengths.
#'       }
#'       \item{\code{seqlengths}}{
#'         Either \code{NULL}, or an integer vector named with
#'         \code{levels(seqnames)} and containing the \code{lengths}
#'         (or \code{NA}) for each level in \code{levels(seqnames)}.
#'       }
#'     }
#'   }
#' }
#'
#' @section Accessors:
#' In the code snippets below, \code{x} is a DuckDBGRangesList object:
#' \describe{
#'   \item{\code{length(x)}:}{
#'     Get the number of elements in \code{x}.
#'   }
#'   \item{\code{names(x)}, \code{names(x) <- value}:}{
#'     Get or set the names of the elements of \code{x}.
#'   }
#'   \item{\code{seqnames(x)}:}{
#'     Get the sequence names as a \linkS4class{DuckDBAtomicList}.
#'   }
#'   \item{\code{start(x)}:}{
#'     Get the start values as a \linkS4class{DuckDBAtomicList}.
#'   }
#'   \item{\code{end(x)}:}{
#'     Get the end values as a \linkS4class{DuckDBAtomicList}.
#'   }
#'   \item{\code{width(x)}:}{
#'     Get the width values as a \linkS4class{DuckDBAtomicList}.
#'   }
#'   \item{\code{strand(x)}:}{
#'     Get the strand values as a \linkS4class{DuckDBAtomicList}.
#'   }
#'   \item{\code{mcols(x)}, \code{mcols(x) <- value}:}{
#'      Get or set the element-level metadata columns.
#'   }
#'   \item{\code{seqinfo(x)}:}{
#'     Get the information about the underlying sequences.
#'   }
#'   \item{\code{elementNROWS(x)}:}{
#'     Get the number of ranges in each element.
#'   }
#' }
#'
#' @section Coercion:
#' In the code snippets below, \code{x} is a DuckDBGRangesList object:
#' \describe{
#'   \item{\code{as(from, "DuckDBDataFrame")}:}{
#'     Creates a \linkS4class{DuckDBDataFrame} object.
#'   }
#'   \item{\code{as(from, "CompressedGRangesList")}:}{
#'     Converts a DuckDBGRangesList object to a CompressedGRangesList object.
#'     This conversion uses SQL UNNEST to flatten the LIST[] columns into a
#'     flat table, creates a GRanges object, then splits it by element to
#'     reconstruct the CompressedGRangesList.
#'   }
#'   \item{\code{realize(x, BACKEND = getAutoRealizationBackend())}:}{
#'     Realize an object into memory or on disk using the equivalent of
#'     \code{realize(as(x, "CompressedGRangesList"), BACKEND)}.
#'   }
#' }
#'
#' @section Subsetting:
#' In the code snippets below, \code{x} is a DuckDBGRangesList object:
#' \describe{
#'   \item{\code{x[i]}:}{
#'     Returns a DuckDBGRangesList object containing the selected elements.
#'   }
#'   \item{\code{x[[i]]}:}{
#'     Return the selected DuckDBGRanges by \code{i}, where \code{i} is an
#'     numeric or character vector of length 1.
#'   }
#'   \item{\code{x$name}:}{
#'     Similar to \code{x[[name]]}, but \code{name} is taken literally as an
#'     element name.
#'   }
#'   \item{\code{head(x, n = 6L)}:}{
#'     If \code{n} is non-negative, returns the first n elements of \code{x}.
#'     If \code{n} is negative, returns all but the last \code{abs(n)} elements
#'     of \code{x}.
#'   }
#'   \item{\code{tail(x, n = 6L)}:}{
#'     If \code{n} is non-negative, returns the last n elements of \code{x}.
#'     If \code{n} is negative, returns all but the first \code{abs(n)} elements
#'     of \code{x}.
#'   }
#' }
#'
#' @author Patrick Aboyoun
#'
#' @examples
#' # Create an example data set with LIST[] columns:
#' df <- data.frame(group = c("gr1", "gr2", "gr3", "gr4"),
#'                  seqnames = I(list("chr2", c("chr2", "chr2", "chr2"),
#'                                    c("chr1", "chr1"),
#'                                    c("chr3", "chr3", "chr3", "chr3"))),
#'                  start = I(list(1L, 2:4, 5:6, 7:10)),
#'                  width = I(list(10L, 9:7, 6:5, 4:1)),
#'                  strand = I(list("-", c("+", "+", "*"), c("+", "+"),
#'                                  c("+", "-", "-", "-"))),
#'                  score = c(1.0, 2.5, 5.5, 8.5),
#'                  GC = c(1.0, 0.7, 0.4, 0.15))
#' tf <- tempfile(fileext = ".parquet")
#' on.exit(unlink(tf))
#' arrow::write_parquet(df, tf)
#'
#' # Create the DuckDBGRangesList object
#' seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")
#' grlist <- DuckDBGRangesList(tf, seqnames = "seqnames", start = "start",
#'                             width = "width", strand = "strand",
#'                             mcols = c("score", "GC"), seqinfo = seqinfo,
#'                             keycol = list(group = c("gr1", "gr2", "gr3", "gr4")))
#' grlist
#'
#' @return
#' The \code{DuckDBGRangesList()} constructor returns a DuckDBGRangesList object.
#' Accessors return the requested component (for example \code{seqnames()},
#' \code{start()}, \code{end()}, \code{width()}, \code{strand()},
#' \code{length()}, \code{names()}, and \code{seqinfo()}); \code{dbconn()} and
#' \code{tblconn()} return the backing DuckDB connection. Replacement methods
#' (such as \code{names<-}, \code{seqlengths<-}, \code{genome<-}, and
#' \code{dimtbls<-}) return the updated DuckDBGRangesList. Coercion methods
#' return the corresponding in-memory object (for example a
#' \link[GenomicRanges]{GRangesList}), and subsetting returns a
#' DuckDBGRangesList.
#'
#' @aliases DuckDBGRangesList-class
#'
#' @aliases updateObject,DuckDBGRangesList-method
#'
#' @aliases dbconn,DuckDBGRangesList-method
#' @aliases tblconn,DuckDBGRangesList-method
#' @aliases keycols,DuckDBGRangesList-method
#' @aliases has_row_number,DuckDBGRangesList-method
#' @aliases dimtbls,DuckDBGRangesList-method
#' @aliases dimtbls<-,DuckDBGRangesList-method
#' @aliases length,DuckDBGRangesList-method
#' @aliases names,DuckDBGRangesList-method
#' @aliases names<-,DuckDBGRangesList-method
#' @aliases seqinfo,DuckDBGRangesList-method
#' @aliases seqnames,DuckDBGRangesList-method
#' @aliases seqlengths<-,DuckDBGRangesList-method
#' @aliases genome<-,DuckDBGRangesList-method
#' @aliases start,DuckDBGRangesList-method
#' @aliases end,DuckDBGRangesList-method
#' @aliases width,DuckDBGRangesList-method
#' @aliases strand,DuckDBGRangesList-method
#' @aliases elementNROWS,DuckDBGRangesList-method
#' @aliases elementMetadata,DuckDBGRangesList-method
#' @aliases elementMetadata<-,DuckDBGRangesList-method
#'
#' @aliases DuckDBGRangesList
#' @aliases split,DuckDBGRanges,DuckDBColumn-method
#'
#' @aliases unlist,DuckDBGRangesList-method
#'
#' @aliases extractROWS,DuckDBGRangesList,ANY-method
#' @aliases getListElement,DuckDBGRangesList-method
#' @aliases head,DuckDBGRangesList-method
#' @aliases tail,DuckDBGRangesList-method
#'
#' @aliases coerce,DuckDBGRangesList,DuckDBDataFrame-method
#' @aliases coerce,DuckDBGRangesList,CompressedGRangesList-method
#' @aliases realize,DuckDBGRangesList-method
#'
#' @aliases show,DuckDBGRangesList-method
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{DuckDBGRanges-class}} for the single GRanges class
#'   \item \code{\link[GenomicRanges]{GRangesList}} for the base class
#' }
#'
#' @include DuckDBGRanges-class.R
#'
#' @keywords classes methods
#'
#' @name DuckDBGRangesList-class
NULL

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importClassesFrom GenomicRanges GRangesList
#' @importFrom S4Vectors new2
setClass("DuckDBGRangesList", contains = "GRangesList",
         slots = c(frame = "DuckDBDataFrame", seqinfo = "Seqinfo"),
         prototype = prototype(elementType = "DuckDBGRanges",
                               frame = new2("DuckDBDataFrame", datacols = .datacols_granges, check = FALSE)))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating
###

#' @export
setMethod("updateObject", "DuckDBGRangesList", function(object, ..., verbose = FALSE) {
    object
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

#' @export
setMethod("dbconn", "DuckDBGRangesList", function(x) callGeneric(x@frame))

#' @export
#' @importFrom DuckDBDataFrame tblconn
setMethod("tblconn", "DuckDBGRangesList", function(x, select = TRUE, filter = TRUE) {
    callGeneric(x@frame, select = select, filter = filter)
})

#' @export
#' @importFrom DuckDBDataFrame keycols
setMethod("keycols", "DuckDBGRangesList", function(x) callGeneric(x@frame))

#' @export
#' @importFrom DuckDBDataFrame has_row_number
setMethod("has_row_number", "DuckDBGRangesList", function(x) callGeneric(x@frame))

#' @export
setMethod("dimtbls", "DuckDBGRangesList", function(x, drop = TRUE) {
    callGeneric(x@frame, drop = drop)
})

#' @export
setReplaceMethod("dimtbls", "DuckDBGRangesList", function(x, value) callGeneric(x@frame, value))

#' @export
setMethod("length", "DuckDBGRangesList", function(x) nrow(x@frame))

#' @export
setMethod("names", "DuckDBGRangesList", function(x) rownames(x@frame))

#' @export
setReplaceMethod("names", "DuckDBGRangesList", function(x, value) {
    rownames(x@frame) <- value
    x
})

#' @export
#' @importFrom DuckDBDataFrame sql_call
#' @importFrom S4Vectors elementNROWS
setMethod("elementNROWS", "DuckDBGRangesList", function(x) {
    nrows <- as.integer(as.vector(sql_call(x@frame[["seqnames"]], "len")))
    setNames(nrows, names(x))
})

#' @export
#' @importFrom Seqinfo seqinfo
setMethod("seqinfo", "DuckDBGRangesList", function(x) x@seqinfo)

#' @export
#' @importFrom S4Vectors elementMetadata
setMethod("elementMetadata", "DuckDBGRangesList", function(x) {
    mcols <- x@elementMetadata
    if (!is.null(mcols) && nrow(mcols) > 0L) {
        rownames(mcols) <- names(x)
    }
    mcols
})

#' @export
#' @importFrom S4Vectors elementMetadata<-
setReplaceMethod("elementMetadata", "DuckDBGRangesList", function(x, value) {
    x@elementMetadata <- value
    x
})

#' @export
#' @importFrom Seqinfo seqnames
setMethod("seqnames", "DuckDBGRangesList", function(x) x@frame[["seqnames"]])

#' @export
#' @importFrom Seqinfo seqinfo seqlengths<-
setReplaceMethod("seqlengths", "DuckDBGRangesList", function (x, value) {
    info <- seqinfo(x)
    seqlengths(info) <- value
    replaceSlots(x, seqinfo = info, check = FALSE)
})

#' @export
#' @importFrom Seqinfo seqinfo genome<-
setReplaceMethod("genome", "DuckDBGRangesList", function (x, value) {
    info <- seqinfo(x)
    genome(info) <- value
    replaceSlots(x, seqinfo = info, check = FALSE)
})

#' @export
setMethod("start", "DuckDBGRangesList", function(x, ...) x@frame[["start"]])

#' @export
setMethod("end", "DuckDBGRangesList", function(x, ...) x@frame[["end"]])

#' @export
setMethod("width", "DuckDBGRangesList", function(x) x@frame[["width"]])

#' @export
setMethod("strand", "DuckDBGRangesList", function(x, ...) x@frame[["strand"]])

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
#' @importFrom DuckDBDataFrame DuckDBDataFrame tblconn
#' @importFrom DBI dbGetQuery
#' @importFrom dplyr sql
#' @importFrom dbplyr sql_render
#' @importFrom GenomicRanges seqinfo
#' @importFrom S4Vectors isSingleString new2
#' @importFrom stats setNames
DuckDBGRangesList <-
function(conn, seqnames, start = NULL, end = NULL, width = NULL, strand = NULL,
         keycol = NULL, dimtbl = NULL, mcols = NULL, seqinfo = NULL,
         seqlengths = NULL)
{
    datacols <- .datacols_granges

    stringAsName <- function(x) if (isSingleString(x)) as.name(x) else x

    nargs_range <- sum(c(!is.null(start), !is.null(end), !is.null(width)))
    if ((nargs_range == 1L) && !is.null(start)) {
        end <- start
        width <- 1L
    } else if (nargs_range != 2L) {
        stop("must provide exactly two of 'start', 'end', and 'width'")
    }

    datacols[["seqnames"]] <- stringAsName(seqnames)
    datacols[["start"]] <- stringAsName(start)
    datacols[["end"]] <- stringAsName(end)
    datacols[["width"]] <- stringAsName(width)
    if (is.null(strand)) {
        datacols[["strand"]] <-
            call("list_transform", datacols[["seqnames"]], sql("x -> '*'"))
    } else {
        datacols[["strand"]] <- stringAsName(strand)
    }

    if (is.null(datacols[["start"]])) {
        datacols[["start"]] <-
            call("list_transform",
                 call("list_zip", datacols[["end"]], datacols[["width"]]),
                 sql("x -> x[1] - x[2] + 1"))
    } else if (is.null(datacols[["end"]])) {
        datacols[["end"]] <-
            call("list_transform",
                 call("list_zip", datacols[["start"]], datacols[["width"]]),
                 sql("x -> x[1] + x[2] - 1"))
    } else if (is.null(datacols[["width"]])) {
        datacols[["width"]] <-
            call("list_transform",
                 call("list_zip", datacols[["end"]], datacols[["start"]]),
                 sql("x -> x[1] - x[2] + 1"))
    }

    datacols <- datacols[c("seqnames", "start", "end", "width", "strand")]
    ccols <- datacols
    if (length(mcols) == 0L) {
        mcols <- NULL
    } else {
        if (is.character(mcols)) {
            mcols <- sapply(mcols, as.name, simplify = FALSE)
        }
        mcols <- as.expression(mcols)
        ccols <- c(ccols, mcols)
    }

    comb <- DuckDBDataFrame(conn, datacols = ccols, keycol = keycol, dimtbl = dimtbl)
    frame <- comb[, names(datacols), drop = FALSE]
    if (is.null(mcols)) {
        mcols <- comb[, character(), drop = FALSE]
    } else {
        mcols <- comb[, names(mcols), drop = FALSE]
    }

    seqinfo <- GenomicRanges:::normarg_seqinfo2(seqinfo, seqlengths)
    if (is.null(seqinfo)) {
        query <- sprintf(
            "SELECT DISTINCT UNNEST(seqnames) AS seqname FROM (%s) ORDER BY seqname",
            sql_render(tblconn(frame, select = FALSE))
        )
        seqlevels <- dbGetQuery(dbconn(frame), query)[[1L]]
        seqinfo <- Seqinfo(seqlevels)
    } else {
        seqinfo <- GenomicRanges:::normarg_seqinfo1(seqinfo)
    }

    new2("DuckDBGRangesList", frame = frame, seqinfo = seqinfo,
         elementMetadata = mcols, check = FALSE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBColumn
#' @importFrom S4Vectors split
setMethod("split", c("DuckDBGRanges", "DuckDBColumn"), function(x, f, drop = FALSE, ...) {
    .NotYetImplemented()
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unlisting
###

#' @export
setMethod("unlist", "DuckDBGRangesList", function(x, recursive = TRUE, use.names = TRUE) {
    .NotYetImplemented()
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

#' @export
#' @importFrom S4Vectors extractROWS
setMethod("extractROWS", "DuckDBGRangesList", function(x, i) {
    if (missing(i)) {
        return(x)
    }

    frame <- callGeneric(x@frame, i = i)

    mcols <- x@elementMetadata
    if (NROW(mcols) > 0L) {
        mcols <- callGeneric(mcols, i = i)
    }

    replaceSlots(x, frame = frame, elementMetadata = mcols, check = FALSE)
})

#' @export
#' @importFrom DuckDBDataFrame DuckDBDataFrame tblconn
#' @importFrom S4Vectors endoapply getListElement new2
setMethod("getListElement", "DuckDBGRangesList", function(x, i, exact = TRUE) {
    frame <- extractROWS(x@frame, i)
    unnest <- names(.datacols_granges)
    frame@datacols <- endoapply(frame@datacols[unnest], function(x) call("UNNEST", x))
    frame <- DuckDBDataFrame(tblconn(frame), datacols = unnest)
    new2("DuckDBGRanges", frame = frame, seqinfo = x@seqinfo, check = FALSE)
})

#' @export
#' @importFrom S4Vectors head
setMethod("head", "DuckDBGRangesList", function(x, n = 6L) {
    if (n < 0L)
        n <- max(length(x) + n, 0L)
    else
        n <- min(n, length(x))
    x[seq_len(n)]
})

#' @export
#' @importFrom S4Vectors tail
setMethod("tail", "DuckDBGRangesList", function(x, n = 6L) {
    len <- length(x)
    if (n < 0L)
        n <- max(len + n, 0L)
    else
        n <- min(n, len)
    x[seq.int(to = len, length.out = n)]
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
setAs("DuckDBGRangesList", "DuckDBDataFrame", function(from) {
    df <- from@frame
    mcols <- from@elementMetadata
    if (!is.null(mcols) && is(mcols, "DuckDBDataFrame") && ncol(mcols) > 0L) {
        df <- cbind(df, mcols)
    }
    df
})

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle split mcols mcols<-
setAs("DuckDBGRangesList", "CompressedGRangesList", function(from) {
    df <- as.data.frame(as(from, "DuckDBDataFrame"))

    seqnames <- unlist(df[["seqnames"]], use.names = FALSE)
    start <- unlist(df[["start"]], use.names = FALSE)
    width <- unlist(df[["width"]], use.names = FALSE)
    strand <- unlist(df[["strand"]], use.names = FALSE)
    gr <- GRanges(seqnames = Rle(seqnames),
                  ranges = IRanges(start = start, width = width),
                  strand = Rle(strand),
                  seqinfo = from@seqinfo)

    group <- rep.int(rownames(df), elementNROWS(from))
    grl <- split(gr, factor(group, levels = names(from)))

    if (NCOL(from@elementMetadata) > 0L) {
        mcols(grl) <- df[, setdiff(colnames(df), names(.datacols_granges)), drop = FALSE]
    }

    grl
})

#' @export
#' @importClassesFrom GenomicRanges CompressedGRangesList
#' @importFrom DelayedArray getAutoRealizationBackend realize
setMethod("realize", "DuckDBGRangesList",
function(x, BACKEND = getAutoRealizationBackend()) {
    x <- as(x, "CompressedGRangesList")
    if (!is.null(BACKEND)) {
        x <- callGeneric(x, BACKEND = BACKEND)
    }
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

#' @export
#' @importFrom S4Vectors classNameForDisplay elementNROWS
setMethod("show", "DuckDBGRangesList", function(object) {
    x_len <- length(object)
    cat(classNameForDisplay(object), " object of length ", x_len, ":\n",
        sep = "")
    cumsumN <- cumsum(elementNROWS(object))
    N <- tail(cumsumN, 1L)
    if (x_len == 0L) {
        cat("<0 elements>\n")
    } else if (x_len <= 3L || (x_len <= 5L && N <= 20L)) {
        ## Display full object.
        show(as.list(object))
    } else {
        ## Display truncated object.
        if (cumsumN[[3L]] <= 20L) {
            showK <- 3L
        } else if (cumsumN[[2L]] <= 20L) {
            showK <- 2L
        } else {
            showK <- 1L
        }
        show(as.list(object[seq_len(showK)]))
        diffK <- x_len - showK
        cat("...\n",
            "<", diffK, " more ", ngettext(diffK, "element", "elements"),
            ">\n", sep = "")
    }
})
