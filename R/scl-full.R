#' scl_full
#'
#' Full spatially-constrained clustering.
#'
#' @param linkage One of \code{"single"}, \code{"average"}, or
#' \code{"complete"}. For covariance clustering, use \code{shortest = FALSE};
#' \code{"complete"} linkage then uses the smallest connecting value between
#' clusters. As for \link{scl_redcap}, initial clusters are always formed from
#' spatially nearest neighbours; \code{shortest} and \code{linkage} determine
#' only how these are then merged according to \code{dmat}.
#' @inheritParams scl_redcap
#'
#' @return A object of class \code{scl} with \code{tree} containing the
#' clustering scheme, and \code{xy} the original coordinate data of the
#' clustered points. An additional component, \code{tree_rest}, enables the tree
#' to be re-cut to a different number of clusters via \link{scl_recluster}.
#'
#' @family clustering_fns
#' @export
#' @examples
#' n <- 100
#' xy <- matrix (runif (2 * n), ncol = 2)
#' dmat <- matrix (runif (n^2), ncol = n)
#' scl <- scl_full (xy, dmat, ncl = 4)
scl_full <- function (xy,
                      dmat,
                      ncl,
                      linkage = "single",
                      shortest = TRUE,
                      nnbs = 6L) {

    linkage <- match.arg (
        tolower (linkage),
        c ("single", "average", "complete")
    )

    if (methods::is (xy, "scl")) {
        message (
            "scl_full is for initial cluster construction; ",
            "passing to scl_recluster"
        )
        scl_recluster_full (xy, ncl = ncl)
    } else {
        scl_full_new (xy, dmat, ncl, linkage, shortest, nnbs)
    }
}

scl_full_new <- function (xy, dmat, ncl, linkage, shortest, nnbs) {

    ids <- rownames (xy)
    xy <- scl_tbl (xy)
    dmat <- dmat_to_square (dmat, nrow (xy), ids)

    # The spatial neighbour graph and initial clusters are always built from
    # the nearest spatial neighbours. `shortest` only determines whether
    # the distances in `dmat` are minimised or maximised when merging.
    if (nnbs <= 0) {
        edges <- scl_edges_tri (xy)
    } else {
        edges <- scl_edges_nn (xy, nnbs = nnbs)
    }

    edges <- full_edges_add_clusters (edges)

    # Then replace the spatial distance in the edges table with the distance
    # from the data to use that as the basis for merging:
    edges <- append_dist_to_edges (edges, dmat, shortest = shortest)

    merges <- full_merges (edges, linkage, shortest)

    nodes <- full_nodes_min_size (edges, merges, ncl)
    tree <- full_realign_tree (edges, nodes)

    pars <- list (
        method = "full",
        ncl = ncl,
        linkage = linkage
    )

    res <- structure (
        list (
            tree = tree,
            merges = merges,
            ord = order_merges (merges),
            nodes = dplyr::bind_cols (nodes, xy),
            pars = pars
        ),
        class = "scl"
    )

    return (scl_statistics (res))
}

full_merges <- function (edges, linkage, shortest) {

    merges <- rcpp_full_merge (
        edges,
        linkage = linkage,
        shortest = shortest
    ) |> data.frame ()

    return (tibble::tibble (
        from = as.integer (merges$from),
        to = as.integer (merges$to),
        dist = merges$dist
    ))
}

# tree at that point has initial cluster numbers which must be re-aligned with
# clusters from the nodal merges:
full_realign_tree <- function (edges, nodes) {

    tree <- edges |> dplyr::select (from, to, d, cluster)
    tree$cluster <- tree$cl_fr <-
        nodes$cluster [match (tree$from, nodes$node)]
    tree$cl_to <- nodes$cluster [match (tree$to, nodes$node)]
    tree$cluster [tree$cl_fr != tree$cl_to] <- NA

    return (dplyr::select (tree, c (from, to, d, cluster)))
}

#' order_merges
#'
#' Order merges so they can be plotted as dendrogram
#' @param merges output from rcpp_full_merge
#' @noRd
order_merges <- function (merges) {

    merges <- as.matrix (merges)
    nodes <- merges [nrow (merges), c ("from", "to")]
    for (i in rev (seq_len (nrow (merges))) [-1]) {
        ii <- which (nodes == merges [i, 2])
        n1 <- n2 <- NULL
        if (ii > 1) {
            n1 <- nodes [1:(ii - 1)]
        }
        if (ii <= length (nodes)) {
            n2 <- nodes [ii:length (nodes)]
        }
        nodes <- c (n1, merges [i, 1], n2)
    }
    return (as.numeric (nodes))
}

#' full_cluster_nodes
#'
#' Transform edge and merge data into rectangle of nodes and cluster IDs
#' @noRd
full_cluster_nodes <- function (edges, merges, ncl) {

    edges$cluster [edges$cluster < 0] <- NA
    ncl_full <- length (unique (edges$cluster))
    merge_tree <- merges [1:(ncl_full - ncl - 1), ]
    for (i in seq_len (nrow (merge_tree))) {
        edges$cluster [edges$cluster == merge_tree$from [i]] <-
            merge_tree$to [i]
    }

    node <- cluster <- NULL # rm undefined variable note
    all_nodes <- unique (c (edges$from, edges$to))
    nodes <- tibble::tibble (
        node = c (edges$from, edges$to),
        cluster = rep (edges$cluster, 2)
    ) |>
        dplyr::distinct () |>
        dplyr::arrange (node) |>
        dplyr::filter (!is.na (cluster))

    # nodes can still be in multiple clusters, so these are set to NA
    dup_nodes <- unique (nodes$node [which (duplicated (nodes$node))])
    nodes$cluster [nodes$node %in% dup_nodes] <- NA_integer_
    nodes <- nodes [which (!duplicated (nodes)), ]

    # Plus nodes entirely in NA clusters can then be removed, and need to be
    # re-inserted:
    na_nodes <- all_nodes [which (!all_nodes %in% nodes$node)]
    if (length (na_nodes) > 0L) {
        na_nodes <- tibble::tibble (
            node = na_nodes,
            cluster = rep (NA_integer_, length (na_nodes))
        )
        nodes <- rbind (nodes, na_nodes)
    }

    # re-order cluster numbers by frequencies
    nt <- sort (table (nodes$cluster), decreasing = TRUE)
    nodes$cluster <- as.integer (names (nt) [match (nodes$cluster, names (nt))])

    return (nodes)
}

#' scl_recluster_full
#'
#' @noRd
scl_recluster_full <- function (scl, ncl = ncl) {

    xy <- scl$nodes |> dplyr::select (x, y)
    num_clusters <- 0
    ncl_trial <- ncl

    while (num_clusters < ncl) {

        scl$nodes <- full_cluster_nodes (scl$tree, scl$merges, ncl_trial)
        num_clusters <- length (which (table (scl$nodes$cluster) > 2))
        ncl_trial <- ncl_trial + 1
        if (ncl_trial >= nrow (scl$nodes)) {
            break
        }
    }
    nt <- sort (table (scl$nodes$cluster), decreasing = TRUE)
    n <- as.integer (names (nt) [which (nt <= 2)])
    scl$nodes$cluster [scl$nodes$cluster %in% n] <- NA

    scl$nodes <- dplyr::bind_cols (scl$nodes, xy)

    return (scl)
}

full_edges_add_clusters <- function (edges) {

    # cluster numbers can be joined with edges through either from or to:
    cl <- as.integer (rcpp_full_initial (edges, TRUE) + 1)

    # make 3 vectors of cluster numbers:
    #   1. cl = cluster number for intra-cluster edges only;
    #   2. cl_from = Num of origin cluster for inter-cluster edges only; and
    #   3. cl_to = Num of destination cluster for inter-cluster edges only.
    from_cl <- cl [edges$from]
    to_cl <- cl [edges$to]
    indx <- which (from_cl == to_cl)
    cl_in <- cl_join_from <- cl_join_to <- rep (NA, nrow (edges))
    cl_in [indx] <- from_cl [indx]
    indx <- which (from_cl != to_cl)
    cl_join_from [indx] <- from_cl [indx]
    cl_join_to [indx] <- to_cl [indx]

    edges$cluster <- cl_in - 1L # convert back to C++ 0-indexed values
    edges$cl_from <- cl_join_from - 1L
    edges$cl_to <- cl_join_to - 1L
    edges$cluster [is.na (edges$cluster)] <- -1L
    edges$cl_from [is.na (edges$cl_from)] <- -1L
    edges$cl_to [is.na (edges$cl_to)] <- -1L

    return (edges)
}

full_nodes_min_size <- function (edges, merges, ncl) {

    # full_cluster_nodes just auto-merges the tree to the specified number,
    # but some of these may be clusters with only 2 members. These are
    # excluded here by iterating until the desired number is achieved in
    # which each cluster has >= 3 members:
    num_clusters <- 0
    ncl_trial <- ncl
    while (num_clusters < ncl) {

        nodes <- full_cluster_nodes (edges, merges, ncl_trial)
        num_clusters <- length (which (table (nodes$cluster) > 2))
        ncl_trial <- ncl_trial + 1
        if (ncl_trial >= nrow (nodes)) {
            break
        }
    }
    nt <- sort (table (nodes$cluster), decreasing = TRUE)
    n <- as.integer (names (nt) [which (nt <= 2)])
    nodes$cluster [nodes$cluster %in% n] <- NA

    return (nodes)
}
