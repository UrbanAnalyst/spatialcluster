# CRAN notes for spatialcluster_0.3.0 submission

This is an initial submission that generates only one NOTE on some CRAN win-builder machines regarding possibly invalid URL. That NOTE nevertheless clearly states it returns a 202, rather than 200, with "Message: Accepted". The URL resolves as expected, but returns 202 to any automated queries.

Other than these, this submission generates no additional notes, and no warnings on:

* Ubuntu: R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
