# spatialcluster 0.2.0.00X (current dev version)

- Fix `shortest` parameter behaviour for `scl_full()` and `scl_redcap()`:
  initial clusters are always formed from spatially nearest neighbours, and
  `shortest` determines only how clusters are then merged (#28).
- `scl_full()` no longer collapses to single shortest edges when
  `shortest = FALSE` (#28).
- Document return value of `scl_full()` and add example for `plot_merges()`.
- Refactor long functions (`plot.scl()`, `scl_full()`, `scl_redcap()`) into
  smaller internal helpers; no change in behaviour.
