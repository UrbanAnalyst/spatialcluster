# C-shaped region; an unconstrained triangulation joins points across the
# gap, a boundary-constrained one does not.
c_region <- function () {
    cbind (
        c (0, 10, 10, 3, 3, 10, 10, 0),
        c (0, 0, 3, 3, 7, 7, 10, 10)
    )
}

c_points <- function (n = 200) {
    b <- c_region ()
    set.seed (5)
    x <- runif (n * 3, 0, 10)
    y <- runif (n * 3, 0, 10)
    in_gap <- x > 3 & y > 3 & y < 7
    cbind (x = x [!in_gap] [seq_len (n)], y = y [!in_gap] [seq_len (n)])
}

# Number of edges with midpoints inside the gap of the C:
n_across_gap <- function (edges, xy) {
    mx <- (xy [edges$from, 1] + xy [edges$to, 1]) / 2
    my <- (xy [edges$from, 2] + xy [edges$to, 2]) / 2
    sum (mx > 3 & my > 3 & my < 7)
}

test_that ("boundary edges do not cross outside the region", {
    skip_if_not_installed ("cdtr")
    xy <- c_points ()
    e_tri <- scl_edges_tri (xy)
    e_cdt <- scl_edges_cdt (xy, c_region ())

    expect_gt (n_across_gap (e_tri, xy), 0L)
    expect_identical (n_across_gap (e_cdt, xy), 0L)

    # symmetric, sorted, and connected:
    expect_identical (nrow (e_cdt) %% 2L, 0L)
    expect_false (is.unsorted (e_cdt$d))
    expect_identical (n_components (e_cdt, nrow (xy)), 1L)
    expect_true (all (seq_len (nrow (xy)) %in% e_cdt$from))
})

test_that ("boundary forms", {
    skip_if_not_installed ("cdtr")
    xy <- c_points ()
    b <- c_region ()
    e0 <- scl_edges_cdt (xy, b)
    # closed ring, data frame, list of rings, and sf-style nested list:
    expect_identical (scl_edges_cdt (xy, rbind (b, b [1, ])), e0)
    expect_identical (scl_edges_cdt (xy, data.frame (b)), e0)
    expect_identical (scl_edges_cdt (xy, list (b)), e0)
    expect_identical (scl_edges_cdt (xy, list (list (rbind (b, b [1, ])))), e0)

    expect_error (scl_edges_cdt (xy, b [1:2, ]), "at least 3 distinct")
    expect_error (scl_edges_cdt (xy, "a"), "must be a matrix")
    expect_error (
        scl_edges_cdt (rbind (xy, c (5, 5)), b),
        "all points must lie inside"
    )
})

test_that ("clustering with boundary", {
    skip_if_not_installed ("cdtr")
    xy <- c_points (100)
    n <- nrow (xy)
    set.seed (1)
    dmat <- matrix (runif (n^2), ncol = n)
    b <- c_region ()

    scl <- scl_redcap (xy, dmat, ncl = 4, boundary = b, quiet = TRUE)
    expect_is (scl, "scl")
    expect_identical (n_across_gap (scl$tree, xy), 0L)

    # 'scl_full' trees retain neighbour edges, so an unconstrained
    # triangulation leaves edges across the gap:
    scl0 <- scl_full (xy, dmat, ncl = 4, nnbs = 0L)
    expect_gt (n_across_gap (scl0$tree, xy), 0L)
    scl <- scl_full (xy, dmat, ncl = 4, boundary = b)
    expect_is (scl, "scl")
    expect_identical (n_across_gap (scl$tree, xy), 0L)
})
