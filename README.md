partitionComparison
===================

An R package which implements many measures for (graph) partition comparison.
For more details, see the package details on CRAN (https://cran.r-project.org/package=partitionComparison)

Changes
-------

0.2.2
  - First version published on CRAN
  - Finalized description, documentation, references
  - Added a 'compareAll' method to run comparison with all measures

0.2.1
  - Dropped igraph dependency in favor of lpSolve
  - Added some tests for method signature registration
  - normalizedMutualInformation now with default "min"
  - Rewritten registration of measure methods for plain partition vectors

0.2.0
  - First finished version
