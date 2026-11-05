#' scl_edges_tri
#'
#' Generate triangulated nearest-neighbour edges between a set of input points.
#' Edges are always defined by, and sorted in ascending order of, **spatial**
#' distances; the user-specified \code{shortest} parameter applies only to the
#' subsequent clustering of \code{dmat} values.
#'
#' @inheritParams scl_redcap
#' @examples
#' set.seed (1)
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' edges <- scl_edges_tri (xy)
#' @noRd
scl_edges_tri <- function (xy) {

    if (!inherits (xy, "data.frame")) {
        xy <- as.data.frame (xy)
    }
    if (ncol (xy) > 2) {
        xy <- dplyr::select (xy, c (x, y))
    }
    names (xy) <- c ("x", "y")
    # Delaunay triangulation via 'deldir'; each segment is returned once, so
    # both directions are added to give a symmetric edge set.
    del <- deldir::deldir (xy$x, xy$y, suppressMsge = TRUE)$delsgs
    edges <- rbind (
        cbind (del$ind1, del$ind2),
        cbind (del$ind2, del$ind1)
    )

    edges <- tibble::tibble (
        from = edges [, 1],
        to = edges [, 2]
    )

    dxy <- as.matrix (stats::dist (xy))

    append_dist_to_edges (edges, dxy, shortest = TRUE)
}

#' scl_edges_cdt
#'
#' Generate triangulated nearest-neighbour edges between a set of input points,
#' constrained by a boundary. The points and the boundary rings are
#' triangulated together with a constrained Delaunay triangulation (via the
#' 'cdtr' package), triangles outside the boundary (and inside any holes) are
#' removed, and only edges joining two input points are retained. This avoids
#' the edges an unconstrained triangulation generates across concavities or
#' gaps in the boundary of the study region.
#'
#' Boundary vertices may take some edges which would otherwise join input
#' points. If the resulting graph is not connected, the edges of a minimal
#' spanning tree of spatial distances are added, as for \link{scl_edges_nn}.
#'
#' @param boundary Boundary rings, as accepted by \link{scl_redcap}.
#' @inheritParams scl_redcap
#' @return A `tibble` of `from` and `to` vertex indices in both directions,
#' along with corresponding spatial distances calculated from 'xy', sorted in
#' ascending order of distance.
#' @noRd
scl_edges_cdt <- function (xy, boundary) {

    if (!requireNamespace ("cdtr", quietly = TRUE)) {
        stop (
            "The 'cdtr' package is required for 'boundary'; ",
            "please install it with 'install.packages (\"cdtr\")'",
            call. = FALSE
        )
    }

    xy <- scl_tbl (xy)
    n <- nrow (xy)
    rings <- boundary_rings (boundary)

    bx <- unlist (lapply (rings, function (r) r [, 1]))
    by <- unlist (lapply (rings, function (r) r [, 2]))
    ring_len <- vapply (rings, nrow, integer (1L))
    ring_start <- n + cumsum (c (0L, ring_len [-length (ring_len)]))
    s0 <- unlist (lapply (seq_along (rings), function (i) {
        ring_start [i] + seq_len (ring_len [i])
    }))
    s1 <- unlist (lapply (seq_along (rings), function (i) {
        ring_start [i] + c (seq_len (ring_len [i]) [-1], 1L)
    }))

    # "holes" removes triangles outside the outer rings and inside any holes:
    tri <- cdtr::cdt_triangulate (
        c (xy$x, bx), c (xy$y, by),
        s0 = s0, s1 = s1,
        erase = "holes"
    )

    # Map triangulation vertices back to input points; boundary vertices (and
    # any Steiner vertices) map to NA.
    pt_index <- match (seq_len (nrow (tri$P)), tri$input_map [seq_len (n)])
    e <- cdtr::cdt_edges (tri)
    e <- cbind (pt_index [e [, 1]], pt_index [e [, 2]])
    e <- e [which (!is.na (e [, 1]) & !is.na (e [, 2])), , drop = FALSE]

    if (!all (seq_len (n) %in% e)) {
        stop (
            "Some points have no neighbours within 'boundary'; ",
            "all points must lie inside the boundary, and must not be ",
            "duplicated",
            call. = FALSE
        )
    }

    edges <- tibble::tibble (
        from = c (e [, 1], e [, 2]),
        to = c (e [, 2], e [, 1])
    )

    dxy <- as.matrix (stats::dist (xy [, c ("x", "y")]))

    if (n_components (edges, n) > 1L) {
        edges_all <- tibble::tibble (
            from = rep (seq_len (n), n),
            to = rep (seq_len (n), each = n),
            d = as.vector (dxy)
        ) |>
            dplyr::arrange (d, from, to) |>
            dplyr::filter (from != to)
        mst <- scl_spantree_ord1 (edges_all)
        edges <- rbind (
            edges,
            mst [, c ("from", "to")],
            tibble::tibble (from = mst$to, to = mst$from)
        )
        edges <- edges [which (!duplicated (edges)), ]
    }

    append_dist_to_edges (edges, dxy, shortest = TRUE)
}

#' boundary_rings
#'
#' Convert a boundary into a list of rings, each a two-column numeric matrix
#' with no repeated closing vertex. Accepts a two-column matrix or data frame
#' (one ring), a list of such rings, or (nested lists of) rings as in 'sf'
#' POLYGON or MULTIPOLYGON geometries, 'sfc' geometry columns, or 'sf' data
#' frames, without requiring the 'sf' package.
#' @noRd
boundary_rings <- function (boundary) {

    if (inherits (boundary, "sf")) {
        boundary <- boundary [[attr (boundary, "sf_column")]]
    }

    flatten <- function (b) {
        if (is.data.frame (b) || is.matrix (b)) {
            return (list (b))
        }
        if (is.list (b)) {
            return (do.call (c, lapply (unclass (b), flatten)))
        }
        stop ("'boundary' must be a matrix, data frame, or list of these",
            call. = FALSE
        )
    }
    rings <- flatten (boundary)

    rings <- lapply (rings, function (r) {
        r <- as.matrix (r)
        if (!is.numeric (r) || ncol (r) < 2L) {
            stop ("'boundary' rings must have at least 2 numeric columns",
                call. = FALSE
            )
        }
        r <- unname (r [, 1:2, drop = FALSE])
        r <- r [stats::complete.cases (r), , drop = FALSE]
        nr <- nrow (r)
        if (nr > 1L && all (r [1, ] == r [nr, ])) {
            r <- r [-nr, , drop = FALSE]
        }
        if (nrow (r) < 3L) {
            stop ("'boundary' rings must have at least 3 distinct vertices",
                call. = FALSE
            )
        }
        r
    })

    return (rings)
}

#' n_components
#'
#' Number of connected components of a graph defined by an edge table with
#' integer 'from' and 'to' columns, over vertices 1..n.
#' @noRd
n_components <- function (edges, n) {

    parent <- seq_len (n)
    find_root <- function (i) {
        while (parent [i] != i) {
            i <- parent [i]
        }
        i
    }
    for (k in seq_len (nrow (edges))) {
        a <- find_root (edges$from [k])
        b <- find_root (edges$to [k])
        if (a != b) {
            parent [a] <- b
        }
    }
    length (unique (vapply (seq_len (n), find_root, integer (1L))))
}

#' scl_edges_nn
#'
#' Generate distance-based nearest-neighbour edges between a set of input
#' points, ensuring that all edges connect to a single component. The minimal
#' spanning tree is constructed from **spatial** distances, not from distances
#' given in `dmat`. Neighbours are always the spatially nearest ones; the
#' \code{shortest} parameter of the clustering functions applies only to
#' \code{dmat} values, and has no influence here.
#' @param nnbs Number of nearest neighbours
#' @param xy Rectangular structure containing coordinates of points.
#'
#' @return A `tibble` of `from` and `to` vertex indices for the minimal spanning
#' tree edges, along with corresponding spatial distances calculated from 'xy'.
#'
#' @examples
#' set.seed (1)
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' scl_edges_nn (xy, nnbs = 3L)
#' @noRd
scl_edges_nn <- function (xy, nnbs) {

    # Initially construct with nnbs + 1, because the `d` matrix includes
    # self-distances of zero, which are subsequently removed.
    nnbs <- nnbs + 1
    d <- apply (
        as.matrix (stats::dist (xy)),
        2,
        function (i) order (i) [seq_len (nnbs)]
    )

    edges <- tibble::tibble (
        from = rep (as.integer (colnames (d)),
            each = nnbs
        ),
        to = as.vector (d)
    )
    # rm self-edges:
    edges <- edges [which (edges$from != edges$to), ]

    # then ensure that the minimal spanning tree is included, to ensure all
    # nearest neighbour edges are connected in a single component. The distances
    # used for this MST are spatial distances, not from `dmat`.
    dxy <- as.matrix (stats::dist (xy))

    n <- nrow (xy)
    edges_all <- tibble::tibble (
        from = rep (seq_len (n), n),
        to = rep (seq_len (n), each = n),
        d = as.vector (dxy)
    ) |>
        dplyr::arrange (d, from, to) |>
        dplyr::filter (from != to)

    mst <- scl_spantree_ord1 (edges_all)
    # duplicate all of those:
    mst <- rbind (
        mst [, c ("from", "to")],
        tibble::tibble (from = mst$to, to = mst$from)
    )

    edges <- rbind (edges, mst)
    edges <- edges [which (!duplicated (edges)), ]

    # Then append final distances from `dxy` to the return value:
    edges <- append_dist_to_edges (edges, dxy, shortest = TRUE)

    return (edges)
}

append_dist_to_edges <- function (edges, dmat, shortest) {
    index <- (edges$to - 1) * nrow (dmat) + edges$from
    edges$d <- dmat [index]

    if (shortest) {
        edges <- dplyr::arrange (edges, d)
    } else {
        edges <- dplyr::arrange (edges, dplyr::desc (d))
    }

    return (edges)
}

#' scl_edges_all
#'
#' Generate full set of edges between a set of input points
#'
#' @param dmat Either a spatial distance matrix generated from 'xy', or a
#' separate distance matrix passed as the 'dmat' parameter to the main
#' \link{scl_redcap} function.
#' @inheritParams scl_redcap
#' @noRd
scl_edges_all <- function (xy, dmat, shortest = TRUE) {

    n <- nrow (dmat)
    edges <- tibble::tibble (
        from = rep (seq (n), times = n),
        to = rep (seq (n), each = n),
        d = as.vector (dmat)
    )
    edges <- na.omit (edges)

    if (shortest) {
        edges <- dplyr::arrange (edges, d)
    } else {
        edges <- dplyr::arrange (edges, dplyr::desc (d))
    }

    return (edges)
}
