#' scl_tbl
#'
#' Convert anything to a tibble
#'
#' @param xy A rectangular object containing the coordinates
#' @return A tibble-ified version of the input, with coordinate columns
#' identified and re-labelled "x" and "y"
#' @noRd
scl_tbl <- function (xy) {
    if (!inherits (xy, "data.frame")) {
        if (!is.numeric (xy)) {
            stop ("coordinates must be numeric")
        }
        if (is.vector (xy)) {
            stop ("coordinates require at least 2 columns")
        }
        xy <- data.frame (xy)
    }
    xi <- grep ("^x|^lon", names (xy), ignore.case = TRUE)
    yi <- grep ("^y|^lat", names (xy), ignore.case = TRUE)
    if (length (xi) == 1L && length (yi) == 1L) {
        names (xy) [xi] <- "x"
        names (xy) [yi] <- "y"
    } else if (ncol (xy) == 2) {
        colnames (xy) <- c ("x", "y")
    } else {
        stop ("Cannot determine unambiguous coordinate columns")
    }
    tibble::as_tibble (xy)
}

#' dmat_is_long
#'
#' Identify whether \code{dmat} is in long form, with three columns of
#' (from, to, value), rather than square with \code{n} rows and columns.
#' @param dmat Distance matrix, either square or long-form
#' @param n Number of points (rows of \code{xy})
#' @noRd
dmat_is_long <- function (dmat, n) {

    if (is.null (dim (dmat))) {
        return (FALSE)
    }
    # square matrices are always taken as such, even for n = 3
    if (nrow (dmat) == n && ncol (dmat) == n) {
        return (FALSE)
    }
    ncol (dmat) == 3L
}

#' dmat_long_to_square
#'
#' Convert a long-form distance matrix to square form. Long-form input has
#' three columns, (from, to, value), where \code{from} and \code{to} are either
#' integer indices into the rows of \code{xy}, or character values matching
#' \code{rownames (xy)}. Pairs may be given in one direction only, in which case
#' values are mirrored. Diagonal entries default to zero.
#' @param dmat Long-form structure (matrix, data.frame, tibble)
#' @param n Number of points (rows of \code{xy})
#' @param ids Optional character vector of point identifiers, used to match
#' non-numeric \code{from} and \code{to} columns.
#' @return Square numeric matrix with \code{n} rows and columns
#' @noRd
dmat_long_to_square <- function (dmat, n, ids = NULL) {

    dmat <- as.data.frame (dmat, stringsAsFactors = FALSE)
    if (ncol (dmat) != 3L) {
        stop ("long-form dmat must have exactly three columns (from, to, value)")
    }
    from <- dmat [[1]]
    to <- dmat [[2]]
    d <- dmat [[3]]
    if (!is.numeric (d)) {
        stop ("third column of long-form dmat must be numeric")
    }

    if (is.numeric (from) && is.numeric (to)) {
        from <- as.integer (from)
        to <- as.integer (to)
    } else {
        if (is.null (ids)) {
            stop (
                "non-numeric 'from' and 'to' values in long-form dmat ",
                "require named rows in 'xy'"
            )
        }
        from <- match (as.character (from), ids)
        to <- match (as.character (to), ids)
    }
    if (anyNA (from) || anyNA (to) ||
        any (c (from, to) < 1L) || any (c (from, to) > n)) {
        stop ("long-form dmat contains identifiers not matching rows of 'xy'")
    }

    res <- matrix (NA_real_, nrow = n, ncol = n)
    res [cbind (from, to)] <- d
    # mirror values given in one direction only
    res_t <- t (res)
    indx <- which (is.na (res) & !is.na (res_t))
    res [indx] <- res_t [indx]
    diag (res) [is.na (diag (res))] <- 0

    if (anyNA (res)) {
        stop ("long-form dmat does not specify values for all pairs of points")
    }

    return (res)
}

#' dmat_to_square
#'
#' Pre-process \code{dmat}, converting from long to square form if necessary.
#' @param dmat Distance matrix, either square or long-form
#' @param n Number of points (rows of \code{xy})
#' @param ids Optional point identifiers (original row names of \code{xy})
#' @noRd
dmat_to_square <- function (dmat, n, ids = NULL) {

    if (dmat_is_long (dmat, n)) {
        dmat <- dmat_long_to_square (dmat, n, ids)
    }
    return (dmat)
}

#' scl_linkage_type
#'
#' Convert \code{linkage} string arg to matching type
#' @param linkage Type of linkage
#' @return Strict match to one of three options
#' @noRd
scl_linkage_type <- function (linkage) {
    linkages <- c ("single", "average", "complete", "full")
    i <- grep (linkage, linkages, ignore.case = TRUE)
    if (length (i) == 0L) {
        stop ("linkage must be one of (single, average, complete, full)")
    }

    return (linkages [i])
}
