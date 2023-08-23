partitionComparison
===================

An R package which implements many measures for (graph) partition comparison.
For more details, see the package details on CRAN (https://cran.r-project.org/package=partitionComparison)

Changes
-------

0.2.6
  - Fixed a test that failed since a change in R-devel (2023-08-16)

0.2.5
  - Changed maintainer information and updated URL/BugReports links

0.2.4
  - Fixed some encoding issues

0.2.3
  - Replaced ``==`` comparisons by tolerance permissive ``isTRUE(all.equal(x, y))``

0.2.2
  - First version published on CRAN
  - Finalized description, documentation, references
  - Added a ``compareAll`` method to run comparison with all measures

0.2.1
  - Dropped igraph dependency in favor of lpSolve
  - Added some tests for method signature registration
  - normalizedMutualInformation now with default "min"
  - Rewritten registration of measure methods for plain partition vectors

0.2.0
  - First finished version
