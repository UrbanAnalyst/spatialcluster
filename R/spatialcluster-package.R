#' Spatially Contiguous Clustering with the 'Redcap' Algorithm
#'
#' Description: Groups spatial data into clusters of similar observations,
#' while requiring that each cluster forms a single spatially contiguous
#' region. Contiguity is enforced by first building a spanning tree from the
#' spatial relationships between points, and then cutting that tree into
#' clusters according to a matrix of distances or other measures of similarity.
#' Trees may be built from first-order relationships only (direct neighbours),
#' yielding a minimum spanning tree, or from full-order relationships between
#' entire clusters, using single linkage (nearest members), average linkage
#' (mean of all pairs), or complete linkage (furthest members). The single and
#' average linkage methods implement the 'Redcap' algorithm ('Regionalization'
#' with dynamically constrained agglomerative clustering and partitioning).
#' Results can be visualised as maps of clusters, cluster hulls, or trees.
#'
#' @name spatialcluster
#' @family package
#' @importFrom Rcpp evalCpp
#' @useDynLib spatialcluster, .registration = TRUE
"_PACKAGE"
