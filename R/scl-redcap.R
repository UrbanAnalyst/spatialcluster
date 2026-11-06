#' scl_redcap
#'
#' Cluster spatial data with REDCAP (REgionalization with Dynamically
#' Constrained Agglomerative clustering and Partitioning) routines.
#'
#' Clustering is in two stages. First, a spanning tree is constructed from the
#' **spatial** relationships between points in \code{xy} alone, using the
#' nearest spatial neighbours (see \code{nnbs}), and independent of both
#' \code{dmat} and \code{shortest}. Second, that tree is cut into \code{ncl}
#' clusters by evaluating the values of \code{dmat} along its edges.
#'
#' @param xy Rectangular structure (matrix, data.frame, tibble), containing
#' coordinates of points to be clustered.
#' @param dmat Structure (matrix, data.frame, tibble) containing distances or
#' equivalent metrics between all points in \code{xy}. If \code{xy} has \code{n}
#' rows, then \code{dmat} may either be square, with \code{n} rows and \code{n}
#' columns, or in long form, with three columns of (\code{from}, \code{to},
#' \code{value}). In long form, \code{from} and \code{to} are either integer row
#' indices of \code{xy} or character values matching \code{rownames (xy)}.
#' Values may be given for one direction of each pair only, and are then
#' mirrored; values for all pairs must be present.
#' @param ncl Desired number of clusters. See description of `ncl_iterate`
#' parameter for conditions under which actual number may be less than this
#' value.
#' @param full_order If \code{FALSE}, build spanning trees from first-order
#' relationships only, otherwise build from full-order relationships (see Note).
#' @param linkage One of \code{"single"}, \code{"average"}, or
#' \code{"complete"}; see Note.
#' @param shortest If \code{TRUE}, the \code{dmat} is interpreted as distances
#' such that lower values indicate stronger relationships, and the tree is cut
#' to minimise the within-cluster sum of squared deviations of \code{dmat}
#' values. If \code{FALSE}, then higher values of \code{dmat} are interpreted
#' to indicate stronger relationships, as is the case for example with
#' covariances, and the tree is cut to maximise the mean \code{dmat} value
#' within clusters. This parameter only affects how the tree is cut; it has no
#' effect on construction of the spanning tree, which is always based on
#' spatial distances.
#' @param nnbs Number of spatially nearest neighbours to be used in calculating
#' clustering trees, irrespective of \code{shortest}. Triangulation will be
#' used if \code{nnbs <= 0}. Ignored if \code{boundary} is given.
#' @param boundary Optional boundary of the study region. If given, spatial
#' neighbours are taken from a constrained Delaunay triangulation of the points
#' within this boundary (requires the 'cdtr' package), so that no neighbour
#' edges cross outside it, for example across concavities in the region. May be
#' a two-column matrix or data frame of the vertices of a single ring, a list of
#' such rings (holes are recognised by nesting), or an 'sf' polygon or
#' multipolygon object. All points in \code{xy} must lie inside the boundary.
#' @param iterate_ncl Actual numbers of clusters found may be less than the
#' specified value of `ncl`, because clusters formed from < 3 edges are removed.
#' If `iterate_ncl = FALSE` (the default), the value is returned with whatever
#' number of actual clusters is found. Setting this parameter to `TRUE` forces
#' the algorithm to iterate until the exact number of clusters has been found.
#' For large data sets, this may result in considerable longer calculation
#' times.
#' @param quiet If `FALSE` (default), display progress information on screen.
#'
#' @return A object of class \code{scl} with \code{tree} containing the
#' clustering scheme, and \code{xy} the original coordinate data of the
#' clustered points. An additional component, \code{tree_rest}, enables the tree
#' to be re-cut to a different number of clusters via \link{scl_recluster},
#' rather than calculating clusters anew.
#'
#' @note Please refer to the original REDCAP paper ('Regionalization with
#' dynamically constrained agglomerative clustering and partitioning (REDCAP)',
#' by D. Guo (2008), Int.J.Geo.Inf.Sci 22:801-823) for details of the
#' \code{full_order} and \code{linkage} parameters. Note that, unlike the
#' original, the spanning trees constructed here (for all values of
#' \code{linkage}) are based on spatial distances between points in \code{xy},
#' with \code{dmat} used only when cutting the trees into clusters. This paper
#' clearly
#' demonstrates the general inferiority of spanning trees constructed from
#' first-order relationships. It is therefore strongly recommended that the
#' default \code{full_order = TRUE} be used at all times.
#'
#' @family clustering_fns
#'
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n^2), ncol = n)
#' scl <- scl_redcap (xy, dmat, ncl = 4)
#' # By default, the spatial tree is cut to minimise within-cluster variation
#' # of (\code{shortest}) values of \code{dmat}. Setting \code{shortest = FALSE}
#' # instead cuts the same tree to maximise mean within-cluster values of
#' # \code{dmat}, as for covariances, and will generally give different
#' # clusters:
#' scl <- scl_redcap (xy, dmat, ncl = 4, shortest = FALSE)
#' # using 'full_order = FALSE' constructs clusters from first-order
#' # relationships only; not recommended, but possible nevertheless:
#' scl <- scl_redcap (xy, dmat, ncl = 4, full_order = FALSE)
#'
#' @export
scl_redcap <- function (xy,
                        dmat,
                        ncl,
                        full_order = TRUE,
                        linkage = "single",
                        shortest = TRUE,
                        nnbs = 6L,
                        boundary = NULL,
                        iterate_ncl = FALSE,
                        quiet = FALSE) {

    linkage <- scl_linkage_type (linkage)

    if (methods::is (xy, "scl")) {

        if (!identical (xy$pars$method, "redcap")) {
            stop (
                "scl_redcap can pass to scl_recluster only for scl objects",
                " previously generated with scl_redcap"
            )
        }

        message (
            "scl_redcap is for initial cluster construction; ",
            "passing to scl_recluster"
        )

        scl_recluster_redcap (xy, ncl = ncl, shortest = shortest)

    } else {

        scl_redcap_new (
            xy, dmat, ncl, full_order, linkage, shortest, nnbs,
            boundary, iterate_ncl, quiet
        )
    }
}

scl_redcap_new <- function (xy, dmat, ncl, full_order, linkage, shortest, nnbs,
                            boundary, iterate_ncl, quiet) {

    ids <- rownames (xy)
    xy <- scl_tbl (xy)
    dmat <- dmat_to_square (dmat, nrow (xy), ids)

    # The spatial structure (neighbours and spanning trees) is always built
    # from the nearest spatial neighbours. `shortest` only determines
    # whether lower or higher values of `dmat` are preferentially merged
    # when cutting the tree.
    if (!is.null (boundary)) {
        edges_nn <- scl_edges_cdt (xy, boundary)
    } else if (nnbs <= 0) {
        edges_nn <- scl_edges_tri (xy)
    } else {
        edges_nn <- scl_edges_nn (xy, nnbs = nnbs)
    }

    tree_full <- redcap_full_tree (xy, edges_nn, full_order, linkage, quiet)

    # Then the critical stage of changing the distance metric on 'edges_nn'
    # from the spatial distances of 'd_xy' to the data-based distances in
    # 'dmat':
    edges_nn <- append_dist_to_edges (edges_nn, dmat, shortest = shortest)

    tree <- scl_cuttree (
        tree_full,
        edges_nn,
        ncl,
        shortest = shortest,
        iterate_ncl = iterate_ncl,
        quiet = quiet
    )

    # meta-data:
    clo <- c ("single", "full") [match (full_order, c (FALSE, TRUE))]
    pars <- list (
        method = "redcap",
        ncl = ncl,
        cl_order = clo,
        linkage = linkage
    )

    res <- structure (
        list (
            tree = tree,
            nodes = dplyr::bind_cols (tree_nodes (tree), xy),
            pars = pars
        ),
        class = "scl"
    )

    return (scl_statistics (res))
}

# Match cluster numbers in edge tree to actual nodes
tree_nodes <- function (tree) {

    node <- NULL # suppress no visible binding note

    res <- tibble::tibble (
        node = c (tree$from, tree$to),
        cluster = rep (tree$cluster, 2)
    ) |>
        dplyr::distinct () |>
        dplyr::arrange (node) |>
        dplyr::filter (!is.na (cluster))

    # remove clusters with < 3 members:
    res$cluster [res$cluster %in% which (table (res$cluster) < 3)] <- NA

    return (res)
}

#' scl_recluster
#'
#' Re-cut a spatial cluster tree (\code{scl}) at a different number of clusters.
#'
#' @param scl An \code{scl} object returned from \link{scl_redcap}.
#' @inheritParams scl_redcap
#'
#' @return Modified \code{scl} object in which \code{tree} is re-cut into
#' \code{ncl} clusters.
#' @family clustering_fns
#'
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n^2), ncol = n)
#' scl <- scl_redcap (xy, dmat, ncl = 4)
#' plot (scl)
#' scl <- scl_recluster (scl, ncl = 5)
#' plot (scl)
#'
#' @export
scl_recluster <- function (scl, ncl, shortest = TRUE, quiet = FALSE) {

    if (!methods::is (scl, "scl")) {
        stop (
            "scl_recluster can only be applied to 'scl' objects ",
            "returned from scl_redcap"
        )
    } else if (identical (scl$pars$method, "redcap")) {
        scl_recluster_redcap (scl = scl, ncl = ncl, shortest = shortest)
    } else if (identical (scl$pars$method, "full")) {
        scl_recluster_full (scl = scl, ncl = ncl)
    }
}

scl_recluster_redcap <- function (scl, ncl, shortest = TRUE, quiet = FALSE) {

    from <- to <- d <- NULL # no visible binding messages

    tree_full <- scl$tree |> dplyr::select (from, to, d)

    if (shortest) {
        tree_full <- dplyr::arrange (tree_full, d)
    } else {
        tree_full <- dplyr::arrange (tree_full, dplyr::desc (d))
    }

    tree_full$cluster <- rcpp_cut_tree (tree_full, ncl,
        shortest = shortest,
        quiet = quiet
    ) + 1

    pars <- scl$pars
    pars$ncl <- ncl

    structure (
        list (
            tree = tree_full,
            nodes = dplyr::bind_cols (
                tree_nodes (tree_full),
                scl$nodes [, c ("x", "y")]
            ),
            pars = pars
        ),
        class = "scl"
    )
}

redcap_full_tree <- function (xy, edges_nn, full_order, linkage, quiet) {

    if (full_order) {

        if (linkage == "average") {

            tree_full <- scl_spantree_alk (edges_nn, shortest = TRUE)

        } else {

            d_xy <- as.matrix (stats::dist (xy))
            edges_all <- scl_edges_all (xy, d_xy, shortest = TRUE)

            if (linkage == "single") {

                tree_full <- scl_spantree_slk (
                    edges_all,
                    edges_nn,
                    shortest = TRUE,
                    quiet = quiet
                )

            } else if (linkage == "complete") {

                tree_full <- scl_spantree_clk (
                    edges_all,
                    edges_nn,
                    shortest = TRUE,
                    quiet = quiet
                )

            } else {

                stop (
                    "linkage must be one of ",
                    "(single, average, complete)"
                )
            }
        }

    } else {

        tree_full <- scl_spantree_ord1 (edges_nn) [, c ("from", "to")]

    }

    return (tree_full)
}
