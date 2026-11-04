context ("full")

test_that ("structure", {
    set.seed (1)
    n <- 100
    xy <- matrix (runif (2 * n), ncol = 2)
    dmat <- matrix (runif (n^2), ncol = n)
    ncl <- 4
    scl <- scl_full (xy, dmat, ncl = ncl)
    expect_is (scl, "scl")
    expect_true (scl$pars$ncl == ncl)
    expect_true (all (names (scl) %in%
        c (
            "tree", "merges", "ord", "nodes",
            "pars", "statistics"
        )))
    cl <- scl$nodes$cluster [!is.na (scl$nodes$cluster)]
    expect_true (length (unique (cl)) == ncl)
})

test_that ("methods", {
    set.seed (1)
    n <- 100
    xy <- matrix (runif (2 * n), ncol = 2)
    dmat <- matrix (runif (n^2), ncol = n)
    ncl <- 8
    scl1 <- scl_full (xy, dmat, ncl = ncl, linkage = "single")
    scl2 <- scl_full (xy, dmat, ncl = ncl, linkage = "average")
    expect_true (!identical (scl1, scl2))
    cl1 <- scl1$nodes$cluster [!is.na (scl1$nodes$cluster)]
    expect_equal (length (unique (cl1)), ncl)
    cl2 <- scl2$nodes$cluster [!is.na (scl2$nodes$cluster)]
    expect_equal (length (unique (cl2)), ncl)

    scl3 <- scl_full (xy, dmat, ncl = ncl, linkage = "complete")
    expect_true (!identical (scl1$merges, scl3$merges))
    expect_true (!identical (scl2$merges, scl3$merges))
    cl3 <- scl3$nodes$cluster [!is.na (scl3$nodes$cluster)]
    expect_equal (length (unique (cl3)), ncl)
})

test_that ("complete linkage merges", {
    set.seed (1)
    n <- 100
    xy <- matrix (runif (2 * n), ncol = 2)
    dmat <- matrix (runif (n^2), ncol = n)
    ncl <- 4

    # shortest: merge distances are non-decreasing under complete linkage
    scl <- scl_full (xy, dmat, ncl = ncl, linkage = "complete")
    expect_equal (scl$pars$linkage, "complete")
    m <- scl$merges
    expect_s3_class (m, "data.frame")
    expect_equal (names (m), c ("from", "to", "dist"))
    expect_true (nrow (m) > 0L)
    expect_true (all (is.finite (m$dist)))
    expect_false (is.unsorted (m$dist))
    # each cluster is merged away (as "from") at most once
    expect_false (any (duplicated (m$from)))
})

# Hand-built case for rcpp_full_merge: 4 initial clusters (0-3), each with one
# intra-cluster edge, joined by inter-cluster edges with a duplicated pair (0-1).
full_merge_test_edges <- function () {
    gr <- data.frame (
        from = 0:8L,
        to = 10:18,
        d = c (0.1, 0.1, 0.1, 0.1, 0.9, 0.3, 0.5, 0.8, 0.2),
        cluster = c (0:3, rep (-1L, 5)),
        cl_from = c (rep (-1L, 4), 0L, 0L, 1L, 2L, 0L),
        cl_to = c (rep (-1L, 4), 1L, 1L, 2L, 3L, 3L)
    )
    gr [order (gr$d), ]
}

test_that ("rcpp_full_merge linkages", {
    gr <- full_merge_test_edges ()

    ms <- rcpp_full_merge (gr, "single", TRUE)
    mc <- rcpp_full_merge (gr, "complete", TRUE)
    mm <- rcpp_full_merge (gr, "max", TRUE)
    expect_equal (colnames (mc), c ("from", "to", "dist"))
    expect_equal (nrow (mc), 3L)
    expect_equal (mc, mm)
    expect_equal (ms [, "dist"], c (0.2, 0.3, 0.5))
    # duplicated 0-1 edge: complete keeps the worst (0.9), single the best (0.3)
    expect_equal (mc [, "dist"], c (0.2, 0.5, 0.9))
    expect_true (all (mc [, "dist"] >= ms [, "dist"]))
    expect_equal (anyDuplicated (mc [, "from"]), 0L)

    # covariance: larger values are better, so complete uses the smallest
    mc <- rcpp_full_merge (gr, "complete", FALSE)
    expect_equal (nrow (mc), 3L)
    expect_equal (mc [, "dist"], c (0.8, 0.5, 0.2))
    expect_equal (anyDuplicated (mc [, "from"]), 0L)

    expect_error (rcpp_full_merge (gr, "nonsense", TRUE), "linkage not found")
})

test_that ("shortest = FALSE (covariance) clustering", {
    set.seed (1)
    n <- 100
    xy <- matrix (runif (2 * n), ncol = 2)
    dmat <- matrix (runif (n^2), ncol = n)
    ncl <- 4
    for (linkage in c ("single", "average", "complete")) {
        scl <- scl_full (xy, dmat, ncl = ncl, linkage = linkage, shortest = FALSE)
        expect_is (scl, "scl")
        cl <- scl$nodes$cluster [!is.na (scl$nodes$cluster)]
        expect_equal (length (unique (cl)), ncl, info = linkage)
    }
    scl <- scl_full (xy, dmat, ncl = ncl, linkage = "complete", shortest = FALSE)
    expect_true (all (is.finite (scl$merges$dist)))
    expect_false (is.unsorted (rev (scl$merges$dist)))
})

test_that ("complete linkage bounded by single linkage", {
    set.seed (2)
    n <- 100
    xy <- matrix (runif (2 * n), ncol = 2)
    dmat <- matrix (runif (n^2), ncol = n)
    ncl <- 4
    ms <- scl_full (xy, dmat, ncl = ncl, linkage = "single")$merges
    mc <- scl_full (xy, dmat, ncl = ncl, linkage = "complete")$merges
    # final complete-linkage height can never be below the single-linkage one
    expect_true (max (mc$dist) >= max (ms$dist))
})

test_that ("invalid linkage", {
    set.seed (1)
    n <- 20
    xy <- matrix (runif (2 * n), ncol = 2)
    dmat <- matrix (runif (n^2), ncol = n)
    expect_error (scl_full (xy, dmat, ncl = 2, linkage = "nonsense"))
})

test_that ("recluster", {
    set.seed (1)
    n <- 100
    xy <- matrix (runif (2 * n), ncol = 2)
    dmat <- matrix (runif (n^2), ncol = n)
    scl <- scl_full (xy, dmat, ncl = 4)
    scl1 <- scl_full (scl, ncl = 3)
    scl2 <- scl_recluster (scl, ncl = 3)
    expect_identical (scl1, scl2)
    expect_error (
        scl3 <- scl_redcap (scl, ncl = 3),
        "scl_redcap can pass to scl_recluster only"
    )
})
