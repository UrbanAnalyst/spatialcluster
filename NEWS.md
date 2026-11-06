# spatialcluster 0.2.0.00X (current dev version)

- New `boundary` parameter for `scl_redcap()` and `scl_full()`: spatial
  neighbours are taken from a constrained Delaunay triangulation within the
  boundary of the study region (via the 'cdtr' package), so that no neighbour
  edges cross outside it (#5).
- Allow 'dmat' to be long-form as well as square (#8)
- Finish main vignette
- Fix `shortest` parameter behaviour for `scl_full()` and `scl_redcap()`:
  initial clusters are always formed from spatially nearest neighbours, and
  `shortest` determines only how clusters are then merged (#28).
- `scl_full()` no longer collapses to single shortest edges when
  `shortest = FALSE` (#28).
- Document return value of `scl_full()` and add example for `plot_merges()`.
- Refactor long functions (`plot.scl()`, `scl_full()`, `scl_redcap()`) into
  smaller internal helpers; no change in behaviour.
