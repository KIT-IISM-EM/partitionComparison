
#' Compute the four coefficients \eqn{N_{11}}, \eqn{N_{10}},
#'  \eqn{N_{01}}, \eqn{N_{00}}
#' 
#' Given two object partitions P and Q, of same length n,
#' each of them described as a vector of cluster ids,
#' compute the four coefficients (\eqn{N_{11}}, \eqn{N_{10}},
#' \eqn{N_{01}}, \eqn{N_{00}})
#' all of the pair comparison measures are based on.
#' 
#' @param p The partition \eqn{P}
#' @param q The partition \eqn{Q}
#'
#' @template author
#' 
#' @examples 
#' pc <- computePairCoefficients(new("Partition", c(0, 0, 0, 1, 1)), 
#'                               new("Partition", c(0, 0, 1, 1, 1)))
#' N11(pc) == 2
#' N10(pc) == 2
#' N01(pc) == 2
#' N00(pc) == 4
#' 
#' @export
computePairCoefficients <- function(p, q) {
  if(length(p) != length(q))
    stop("Both partitions must be of the same set")
  
  N11 <- 0
  N10 <- 0
  N01 <- 0
  N00 <- 0
  
  for (i in 1:(length(p)-1)) {
    for (j in (i+1):length(p)) {
      if (p[i] == p[j] & q[i] == q[j]) {
        N11 <- N11 + 1
      } else if (p[i] == p[j] & q[i] != q[j]) {
        N10 <- N10 + 1
      } else if (p[i] != p[j] & q[i] == q[j]) {
        N01 <- N01 + 1
      } else {
        N00 <- N00 + 1
      }
    }
  }
  
  # Assert: the sum of the coefficients must be n choose 2
  stopifnot((N00 + N10 + N01 + N11 == choose(length(p), 2)))
  
  new("PairCoefficients", N00=N00, N10=N10, N01=N01, N11=N11)
}


#' Rand Index
#' 
#' Compute the Rand index
#' \deqn{\frac{N_{11} + N_{00}}{N}}
#' 
#' @references
#' \insertRef{Rand1971}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name randIndex
#' @examples 
#' randIndex(new("Partition", c(0, 0, 0, 1, 1)), 
#'           new("Partition", c(0, 0, 1, 1, 1))) == 0.6
#' 
#' @export
setGeneric("randIndex", function(p, q) standardGeneric("randIndex"))

#' @describeIn randIndex Compute given two partitions
setMethod("randIndex", signature(p="Partition", q="Partition"),
          function(p, q) randIndex(computePairCoefficients(p, q)))

#' @describeIn randIndex Compute given the pair coefficients
setMethod("randIndex", signature(p="PairCoefficients", q="missing"), 
          function(p, q=NULL) (N11(p) + N00(p))/N(p))


#' Adjusted Rand Index
#' 
#' Compute the Adjusted Rand Index (ARI)
#' \deqn{\frac{2(N_{00}N_{11} - N_{10}N_{01})}{N'_{01}N_{12} + N'_{10}N_{21}}}
#' 
#' @references
#' \insertRef{Hubert1985}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name adjustedRandIndex
#' @examples 
#' adjustedRandIndex(new("Partition", c(0, 0, 0, 1, 1)), 
#'                   new("Partition", c(0, 0, 1, 1, 1))) == 1/6
#' 
#' @export
setGeneric("adjustedRandIndex", function(p, q) standardGeneric("adjustedRandIndex"))

#' @describeIn adjustedRandIndex Compute given two partitions
setMethod("adjustedRandIndex", signature(p="Partition", q="Partition"),
          function(p, q) adjustedRandIndex(computePairCoefficients(p, q)))

#' @describeIn adjustedRandIndex Compute given the pair coefficients
setMethod("adjustedRandIndex", signature(p="PairCoefficients", q="missing"), 
          function(p, q=NULL) {
            2*(N11(p)*N00(p) - N10(p)*N01(p)) / (N01p(p)*N12(p) + N10p(p)*N21(p))
          })


#' Jaccard Coefficient
#' 
#' Compute the Jaccard coefficient
#' \deqn{\frac{N_{11}}{N_{11} + N_{10} + N_{01}}}
#' 
#' @references
#' \insertRef{Jaccard1908}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name jaccardCoefficient
#' @examples 
#' jaccardCoefficient(new("Partition", c(0, 0, 0, 1, 1)), 
#'                    new("Partition", c(0, 0, 1, 1, 1))) == 1/3
#' 
#' @export
setGeneric("jaccardCoefficient", function(p, q) standardGeneric("jaccardCoefficient"))

#' @describeIn jaccardCoefficient Compute given two partitions
setMethod("jaccardCoefficient", signature(p="Partition", q="Partition"),
          function(p, q) jaccardCoefficient(computePairCoefficients(p, q)))

#' @describeIn jaccardCoefficient Compute given the pair coefficients
setMethod("jaccardCoefficient", signature(p="PairCoefficients", q="missing"), 
          function(p, q=NULL) N11(p)/(N11(p) + N10(p) + N01(p)))


#' Wallace I
#' 
#' Compute Wallace' index I
#' \deqn{\frac{N_{11}}{N_{21}}}
#' 
#' @references
#' \insertRef{Wallace1983}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name wallaceI
#' @examples 
#' wallaceI(new("Partition", c(0, 0, 0, 1, 1)), new("Partition", c(0, 0, 1, 1, 1))) == 0.5
#' 
#' @seealso \code{\link{folwkesMallowsIndex}}
#' @export
setGeneric("wallaceI", function(p, q) standardGeneric("wallaceI"))

#' @describeIn wallaceI Compute given two partitions
setMethod("wallaceI", signature(p="Partition", q="Partition"),
          function(p, q) wallaceI(computePairCoefficients(p, q)))

#' @describeIn wallaceI Compute given the pair coefficients
setMethod("wallaceI", signature(p="PairCoefficients", q="missing"), 
          function(p, q=NULL) N11(p)/(N11(p) + N10(p)))


#' Wallace II
#' 
#' Compute Wallace' index II
#' \deqn{\frac{N_{11}}{N_{12}}}
#' 
#' @references
#' \insertRef{Wallace1983}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name wallaceII
#' @examples 
#' wallaceII(new("Partition", c(0, 0, 0, 1, 1)), 
#'           new("Partition", c(0, 0, 1, 1, 1))) == 0.5
#' 
#' @seealso \code{\link{folwkesMallowsIndex}}
#' @export
setGeneric("wallaceII", function(p, q) standardGeneric("wallaceII"))

#' @describeIn wallaceII Compute given two partitions
setMethod("wallaceII", signature(p="Partition", q="Partition"),
          function(p, q) wallaceII(computePairCoefficients(p, q)))

#' @describeIn wallaceII Compute given the pair coefficients
setMethod("wallaceII", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) N11(p)/(N11(p) + N01(p)))


#' Folwkes & Mallows Index
#' 
#' Compute the index of Folwkes and Mallows
#' \deqn{\sqrt{\frac{N_{11}}{N_{21}} \frac{N_{11}}{N_{12}}}}
#' which is a combination of the two Wallace indices.
#' 
#' @references
#' \insertRef{Fowlkes1983}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name folwkesMallowsIndex
#' @examples 
#' folwkesMallowsIndex(new("Partition", c(0, 0, 0, 1, 1)), 
#'                     new("Partition", c(0, 0, 1, 1, 1))) == 0.5
#' 
#' @seealso \code{\link{wallaceI}} \code{\link{wallaceII}}
#' @export
setGeneric("folwkesMallowsIndex", function(p, q) standardGeneric("folwkesMallowsIndex"))

#' @describeIn folwkesMallowsIndex Compute given two partitions
setMethod("folwkesMallowsIndex", signature(p="Partition", q="Partition"),
          function(p, q) folwkesMallowsIndex(computePairCoefficients(p, q)))

#' @describeIn folwkesMallowsIndex Compute given the pair coefficients
setMethod("folwkesMallowsIndex", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) sqrt(wallaceI(p) * wallaceII(p)))


#' RV Coefficient
#' 
#' Compute the RV coefficient
#' \deqn{\frac{n + 2N_{11}(p)}{\sqrt{(2N_{21}(p) + n) (2N_{12}(p) + n)}}}
#' 
#' @references
#' \insertRef{Robert1976}{partitionComparison}
#' 
#' \insertRef{Youness2004}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name rvCoefficient
#' @examples 
#' rvCoefficient(new("Partition", c(0, 0, 0, 1, 1)), 
#'               new("Partition", c(0, 0, 1, 1, 1))) == 9 / 13
#' 
#' @export
setGeneric("rvCoefficient", function(p, q) standardGeneric("rvCoefficient"))

#' @describeIn rvCoefficient Compute the RV coefficient given two partitions
setMethod("rvCoefficient", signature(p="Partition", q="Partition"),
          function(p, q) rvCoefficient(p=computePairCoefficients(p, q)))

#' @describeIn rvCoefficient Compute the RV coefficient given the pair coefficients
setMethod("rvCoefficient", 
          signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) {
            n <- (1 + sqrt(1 + 8 * N(p))) / 2  # Inverse of n*(n-1)/2 = N11+...+N00
            
            (n + 2 * N11(p)) / sqrt((2 * N21(p)+ n) * (2 * N12(p)+ n))
          })


#' Mirkin Metric
#' 
#' Compute the Mirkin metric
#' \deqn{2(N_{10} + N_{01})}
#' 
#' @references
#' \insertRef{Mirkin1970}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name mirkinMetric
#' @examples 
#' mirkinMetric(new("Partition", c(0, 0, 0, 1, 1)), 
#'              new("Partition", c(0, 0, 1, 1, 1))) == 8
#' 
#' @export
setGeneric("mirkinMetric", function(p, q) standardGeneric("mirkinMetric"))

#' @describeIn mirkinMetric Compute given two partitions
setMethod("mirkinMetric", signature(p="Partition", q="Partition"),
          function(p, q) mirkinMetric(computePairCoefficients(p, q)))

#' @describeIn mirkinMetric Compute given the pair coefficients
setMethod("mirkinMetric", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) 2 * (N01(p) + N10(p)))


#' Minkowski Measure
#' 
#' Compute the Minkowski measure
#' \deqn{\sqrt{ \frac{N_{10} + N_{01}}{N_{11} + N_{10}} }}
#' 
#' @references
#' \insertRef{Minkowski1911}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name minkowskiMeasure
#' @examples 
#' minkowskiMeasure(new("Partition", c(0, 0, 0, 1, 1)), 
#'                  new("Partition", c(0, 0, 1, 1, 1))) == 1
#' 
#' @export
setGeneric("minkowskiMeasure", function(p, q) standardGeneric("minkowskiMeasure"))

#' @describeIn minkowskiMeasure Compute given two partitions
setMethod("minkowskiMeasure", signature(p="Partition", q="Partition"),
          function(p, q) minkowskiMeasure(computePairCoefficients(p, q)))

#' @describeIn minkowskiMeasure Compute given the pair coefficients
setMethod("minkowskiMeasure", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) sqrt((N10(p) + N01(p)) / (N10(p) + N11(p))))


#' Gamma Statistics
#' 
#' Compute the Gamma statistics
#' \deqn{\frac{N_{11}N_{00} - N_{10}N_{01}}{\sqrt{ N_{21}N_{12}N'_{10}N'_{01} }}}
#' 
#' @references
#' \insertRef{Yule1900}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name gammaStatistics
#' @examples 
#' gammaStatistics(new("Partition", c(0, 0, 0, 1, 1)), 
#'                 new("Partition", c(0, 0, 1, 1, 1))) == 1/6
#' 
#' @export
setGeneric("gammaStatistics", function(p, q) standardGeneric("gammaStatistics"))

#' @describeIn gammaStatistics Compute given two partitions
setMethod("gammaStatistics", signature(p="Partition", q="Partition"),
          function(p, q) gammaStatistics(computePairCoefficients(p, q)))

#' @describeIn gammaStatistics Compute given the pair coefficients
setMethod("gammaStatistics", signature(p="PairCoefficients", q="missing"),
          # function(p, q=NULL) {
          #   (N(p) * N11(p) - N21(p) * N12(p)) / 
          #     sqrt(N21(p) * N12(p) * (N(p) - N12(p))  * (N(p) - N21(p)))
          #   })
          function(p, q=NULL) {
            (N11(p) * N00(p) - N10(p) * N01(p)) / 
              (sqrt(N21(p)) * sqrt(N12(p)) * sqrt(N10p(p)) * sqrt(N01p(p)))
          })


#' Hamann Coefficient
#' 
#' Compute the Hamann coefficient
#' \deqn{\frac{(N_{11} + N_{00}) - (N_{10} + N_{01})}{N}}
#' 
#' @references
#' \insertRef{Hamann1961}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name hamann
#' @examples 
#' hamann(new("Partition", c(0, 0, 0, 1, 1)), new("Partition", c(0, 0, 1, 1, 1))) == 0.2
#' 
#' @export
setGeneric("hamann", function(p, q) standardGeneric("hamann"))

#' @describeIn hamann Compute given two partitions
setMethod("hamann", signature(p="Partition", q="Partition"),
          function(p, q) hamann(computePairCoefficients(p, q)))

#' @describeIn hamann Compute given the pair coefficients
setMethod("hamann", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) (N11(p) + N00(p) - N01(p) - N10(p)) / N(p))


#' Czekanowski Index
#' 
#' Compute the Czekanowski index
#' \deqn{\frac{2N_{11}}{2N_{11} + N_{10} + N_{01}}}
#' 
#' @references
#' \insertRef{Czekanowski1932}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name czekanowski
#' @examples 
#' czekanowski(new("Partition", c(0, 0, 0, 1, 1)), 
#'             new("Partition", c(0, 0, 1, 1, 1))) == 0.5
#' 
#' @export
setGeneric("czekanowski", function(p, q) standardGeneric("czekanowski"))

#' @describeIn czekanowski Compute given two partitions
setMethod("czekanowski", signature(p="Partition", q="Partition"),
          function(p, q) czekanowski(computePairCoefficients(p, q)))

#' @describeIn czekanowski Compute given the pair coefficients
setMethod("czekanowski", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) 2 * N11(p) / (2 * N11(p) + N01(p) + N10(p)))


#' Kulczynski Index
#' 
#' Compute the Kulczynski index
#' \deqn{\frac{1}{2} \left(\frac{N_{11}}{N_{21}} + \frac{N_{11}}{N_{12}} \right)}
#' 
#' @references
#' \insertRef{Kulczynski1927}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name kulczynski
#' @examples 
#' kulczynski(new("Partition", c(0, 0, 0, 1, 1)), 
#'            new("Partition", c(0, 0, 1, 1, 1))) == 0.5
#' 
#' @export
setGeneric("kulczynski", function(p, q) standardGeneric("kulczynski"))

#' @describeIn kulczynski Compute given two partitions
setMethod("kulczynski", signature(p="Partition", q="Partition"),
          function(p, q) kulczynski(computePairCoefficients(p, q)))

#' @describeIn kulczynski Compute given the pair coefficients
setMethod("kulczynski", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) (N11(p) / N21(p) + N11(p) / N12(p)) / 2)


#' McConnaughey Index
#' 
#' Compute the McConnaughey index
#' \deqn{\frac{N_{11}^2 - N_{10}N_{01}}{N_{21}N_{12}}}
#' 
#' @references
#' \insertRef{McConnaughey1964}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name mcconnaughey
#' @examples 
#' mcconnaughey(new("Partition", c(0, 0, 0, 1, 1)), 
#'              new("Partition", c(0, 0, 1, 1, 1))) == 0
#' 
#' @export
setGeneric("mcconnaughey", function(p, q) standardGeneric("mcconnaughey"))

#' @describeIn mcconnaughey Compute given two partitions
setMethod("mcconnaughey", signature(p="Partition", q="Partition"),
          function(p, q) mcconnaughey(computePairCoefficients(p, q)))

#' @describeIn mcconnaughey Compute given the pair coefficients
setMethod("mcconnaughey", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) (N11(p)**2 - N10(p) * N01(p)) / 
            (N21(p) * N12(p)))


#' Peirce Index
#' 
#' Compute the Peirce index
#' \deqn{\frac{N_{11}N_{00} - N_{10}N_{01}}{N_{21}N'_{01}}}
#' 
#' @references
#' \insertRef{Peirce1884}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name peirce
#' @examples 
#' peirce(new("Partition", c(0, 0, 0, 1, 1)), new("Partition", c(0, 0, 1, 1, 1))) == 1/6
#' 
#' @export
setGeneric("peirce", function(p, q) standardGeneric("peirce"))

#' @describeIn peirce Compute given two partitions
setMethod("peirce", signature(p="Partition", q="Partition"),
          function(p, q) peirce(computePairCoefficients(p, q)))

#' @describeIn peirce Compute given the pair coefficients
setMethod("peirce", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) (N11(p) * N00(p) - N10(p) * N01(p)) / (N21(p) * N01p(p)))


#' Sokal & Sneath Index 1
#' 
#' Compute the index 1 of Sokal and Sneath
#' \deqn{
#' \frac{1}{4} \left( \frac{N_{11}}{N_{21}} + \frac{N_{11}}{N_{12}} + 
#' \frac{N_{00}}{N'_{10}} + \frac{N_{00}}{N'_{01}} \right)
#' }
#' 
#' @references
#' \insertRef{Sokal1963}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name sokalSneath1
#' @examples 
#' sokalSneath1(new("Partition", c(0, 0, 0, 1, 1)), 
#'              new("Partition", c(0, 0, 1, 1, 1))) == 7/12
#' 
#' @export
setGeneric("sokalSneath1", function(p, q) standardGeneric("sokalSneath1"))

#' @describeIn sokalSneath1 Compute given two partitions
setMethod("sokalSneath1", signature(p="Partition", q="Partition"),
          function(p, q) sokalSneath1(computePairCoefficients(p, q)))

#' @describeIn sokalSneath1 Compute given the pair coefficients
setMethod("sokalSneath1", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) (N11(p) / N21(p) + N11(p) / N12(p) + 
                                 N00(p) / N10p(p) + N00(p) / N01p(p)) / 4)


#' Sokal & Sneath Index 2
#' 
#' Compute the index 2 of Sokal and Sneath
#' \deqn{\frac{N_{11}}{N_{11} + 2(N_{10} + N_{01})}}
#' 
#' @references
#' \insertRef{Sokal1963}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name sokalSneath2
#' @examples 
#' sokalSneath2(new("Partition", c(0, 0, 0, 1, 1)), 
#'              new("Partition", c(0, 0, 1, 1, 1))) == 0.2
#' 
#' @export
setGeneric("sokalSneath2", function(p, q) standardGeneric("sokalSneath2"))

#' @describeIn sokalSneath2 Compute given two partitions
setMethod("sokalSneath2", signature(p="Partition", q="Partition"),
          function(p, q) sokalSneath2(computePairCoefficients(p, q)))

#' @describeIn sokalSneath2 Compute given the pair coefficients
setMethod("sokalSneath2", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) N11(p) / (N11(p) + 2 * N01(p) + 2 * N10(p)))


#' Sokal & Sneath Index 3
#' 
#' Compute the index 3 of Sokal and Sneath
#' \deqn{\frac{N_{11}N_{00}}{\sqrt{N_{21}N_{12}N'_{01}N'_{10}}}}
#' 
#' @references
#' \insertRef{Sokal1963}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name sokalSneath3
#' @examples 
#' sokalSneath3(new("Partition", c(0, 0, 0, 1, 1)), 
#'              new("Partition", c(0, 0, 1, 1, 1))) == 1/3
#' 
#' @export
setGeneric("sokalSneath3", function(p, q) standardGeneric("sokalSneath3"))

#' @describeIn sokalSneath3 Compute given two partitions
setMethod("sokalSneath3", signature(p="Partition", q="Partition"),
          function(p, q) sokalSneath3(computePairCoefficients(p, q)))

#' @describeIn sokalSneath3 Compute given the pair coefficients
setMethod("sokalSneath3", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) (N11(p) * N00(p)) / 
            (sqrt(N21(p)) * sqrt(N12(p)) * sqrt(N10p(p)) * sqrt(N01p(p))))


#' Baulieu Index 1
#' 
#' Compute the index 1 of Baulieu
#' \deqn{\frac{ N^2 - N(N_{10} + N_{01}) + (N_{10} - N_{01})^2 }{ N^2 }}
#' 
#' @references
#' \insertRef{Baulieu1989}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name baulieu1
#' @examples 
#' baulieu1(new("Partition", c(0, 0, 0, 1, 1)), 
#'          new("Partition", c(0, 0, 1, 1, 1))) == 0.76
#' 
#' @export
setGeneric("baulieu1", function(p, q) standardGeneric("baulieu1"))

#' @describeIn baulieu1 Compute given two partitions
setMethod("baulieu1", signature(p="Partition", q="Partition"),
          function(p, q) baulieu1(computePairCoefficients(p, q)))

#' @describeIn baulieu1 Compute given the pair coefficients
setMethod("baulieu1", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) (N(p)**2 - N(p) * (N10(p) + N01(p)) + 
                                 (N10(p) + N01(p))**2) / N(p)**2)


#' Baulieu Index 2
#' 
#' Compute the index 2 of Baulieu
#' \deqn{\frac{ N_{11}N_{00} - N_{10}N_{01} }{ N^2 }}
#' 
#' @references
#' \insertRef{Baulieu1989}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name baulieu2
#' @examples 
#' baulieu2(new("Partition", c(0, 0, 0, 1, 1)), 
#'          new("Partition", c(0, 0, 1, 1, 1))) == 0.04
#' 
#' @export
setGeneric("baulieu2", function(p, q) standardGeneric("baulieu2"))

#' @describeIn baulieu2 Compute given two partitions
setMethod("baulieu2", signature(p="Partition", q="Partition"),
          function(p, q) baulieu2(computePairCoefficients(p, q)))

#' @describeIn baulieu2 Compute given the pair coefficients
setMethod("baulieu2", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) (N11(p) * N00(p) - N10(p) * N01(p)) / N(p)**2)


#' Russel & Rao Index
#' 
#' Compute the index of Russel and Rao
#' \deqn{\frac{N_{11}}{N}}
#' 
#' @references
#' \insertRef{Russel1940}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name russelRao
#' @examples 
#' russelRao(new("Partition", c(0, 0, 0, 1, 1)), 
#'           new("Partition", c(0, 0, 1, 1, 1))) == 0.2
#' 
#' @export
setGeneric("russelRao", function(p, q) standardGeneric("russelRao"))

#' @describeIn russelRao Compute given two partitions
setMethod("russelRao", signature(p="Partition", q="Partition"),
          function(p, q) russelRao(computePairCoefficients(p, q)))

#' @describeIn russelRao Compute given the pair coefficients
setMethod("russelRao", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) N11(p) / N(p))


#' Fager & McGowan Index
#' 
#' Compute the index of Fager and McGowan
#' \deqn{\frac{N_{11}}{\sqrt{N_{21}N_{12}}} - \frac{1}{2\sqrt{N_{21}}}}
#' 
#' @references
#' \insertRef{Fager1963}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name fagerMcGowan
#' @examples 
#' fagerMcGowan(new("Partition", c(0, 0, 0, 1, 1)), 
#'              new("Partition", c(0, 0, 1, 1, 1))) == 0.25
#' 
#' @export
setGeneric("fagerMcGowan", function(p, q) standardGeneric("fagerMcGowan"))

#' @describeIn fagerMcGowan Compute given two partitions
setMethod("fagerMcGowan", signature(p="Partition", q="Partition"),
          function(p, q) fagerMcGowan(computePairCoefficients(p, q)))

#' @describeIn fagerMcGowan Compute given the pair coefficients
setMethod("fagerMcGowan", signature(p="PairCoefficients", q="missing"),
          function(p, q=NULL) N11(p) / sqrt(N21(p) * N12(p)) - 1 / (2 * sqrt(N21(p))))


#' Pearson Index
#' 
#' Compute the Pearson index
#' \deqn{\frac{N_{11}N_{00} - N_{10}N_{01}}{N_{21}N_{12}N'_{01}N'_{10}}}
#' 
#' @references
#' \insertRef{Pearson1926}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name pearson
#' @examples 
#' pearson(new("Partition", c(0, 0, 0, 1, 1)), 
#'         new("Partition", c(0, 0, 1, 1, 1))) == 1/144
#' 
#' @export
setGeneric("pearson", function(p, q) standardGeneric("pearson"))

#' @describeIn pearson Compute given two partitions
setMethod("pearson", signature(p="Partition", q="Partition"),
          function(p, q) pearson(computePairCoefficients(p, q)))

#' @describeIn pearson Compute given the pair coefficients
setMethod("pearson", signature(p="PairCoefficients", q="missing"),
          function(p, q) (N11(p) * N00(p) - N10(p) * N01(p)) / 
            (N21(p) * N12(p) * N10p(p) * N01p(p)))


#' Gower & Legendre Index
#' 
#' Compute the index of Gower and Legendre
#' \deqn{
#' \frac{N_{11} + N_{00}}{N_{11} + \frac{1}{2}\left(N_{10} + N_{01}\right) + N_{00}}
#' }
#' 
#' @references
#' \insertRef{Gower1986}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name gowerLegendre
#' @examples 
#' gowerLegendre(new("Partition", c(0, 0, 0, 1, 1)), 
#'               new("Partition", c(0, 0, 1, 1, 1))) == 0.75
#' 
#' @export
setGeneric("gowerLegendre", function(p, q) standardGeneric("gowerLegendre"))

#' @describeIn gowerLegendre Compute given two partitions
setMethod("gowerLegendre", signature(p="Partition", q="Partition"),
          function(p, q) gowerLegendre(computePairCoefficients(p, q)))

#' @describeIn gowerLegendre Compute given the pair coefficients
setMethod("gowerLegendre", signature(p="PairCoefficients", q="missing"),
          function(p, q) (N11(p) + N00(p)) / 
            (N11(p) + N00(p) + N10(p) / 2 + N01(p) / 2))


#' Rogers & Tanimoto Index
#' 
#' Compute the index of Rogers and Tanimoto
#' \deqn{\frac{N_{11} + N_{00}}{N_{11} + 2(N_{10} + N_{01}) + N_{00}}}
#' 
#' @references
#' \insertRef{Rogers1960}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name rogersTanimoto
#' @examples 
#' rogersTanimoto(new("Partition", c(0, 0, 0, 1, 1)), 
#'                new("Partition", c(0, 0, 1, 1, 1))) == 3/7
#' 
#' @export
setGeneric("rogersTanimoto", function(p, q) standardGeneric("rogersTanimoto"))

#' @describeIn rogersTanimoto Compute given two partitions
setMethod("rogersTanimoto", signature(p="Partition", q="Partition"),
          function(p, q) rogersTanimoto(computePairCoefficients(p, q)))

#' @describeIn rogersTanimoto Compute given the pair coefficients
setMethod("rogersTanimoto", signature(p="PairCoefficients", q="missing"),
          function(p, q) (N11(p) + N00(p)) / 
            (N11(p) + N00(p) + N10(p) * 2 + N01(p) * 2))


#' Goodman & Kruskal Index
#' 
#' Compute the index of Goodman and Kruskal
#' \deqn{\frac{N_{11}N_{00} - N_{10}N_{01}}{N_{11}N_{00} + N_{10}N_{01}}}
#' 
#' @references
#' \insertRef{Goodman1954}{partitionComparison}
#'
#' @template pair_comp_params
#' @template author
#' @name goodmanKruskal
#' @examples 
#' goodmanKruskal(new("Partition", c(0, 0, 0, 1, 1)), 
#'                new("Partition", c(0, 0, 1, 1, 1))) == 1/3
#' 
#' @export
setGeneric("goodmanKruskal", function(p, q) standardGeneric("goodmanKruskal"))

#' @describeIn goodmanKruskal Compute given two partitions
setMethod("goodmanKruskal", signature(p="Partition", q="Partition"),
          function(p, q) goodmanKruskal(computePairCoefficients(p, q)))

#' @describeIn goodmanKruskal Compute given the pair coefficients
setMethod("goodmanKruskal", signature(p="PairCoefficients", q="missing"),
          function(p, q) (N11(p) * N00(p) - N10(p) * N01(p)) / 
            (N11(p) * N00(p) + N10(p) * N01(p)))

#' Lerman Index
#' 
#' Compute the Lerman index
#' \deqn{\frac{N_{11} - E(N_{11})}{\sqrt{\sigma^2(N_{11})}}}
#' 
#' @references
#' \insertRef{Lerman1988}{partitionComparison}
#' 
#' \insertRef{Hubert1985}{partitionComparison}
#' 
#' \insertRef{Deneud2006}{partitionComparison}
#' 
#' @param p The partition \eqn{P}
#' @param q The partition \eqn{Q}
#' @param c \linkS4class{PairCoefficients} or NULL
#' @template author
#' @name lermanIndex
#' @examples 
#' lermanIndex(new("Partition", c(0, 0, 0, 1, 1)), 
#'             new("Partition", c(0, 0, 1, 1, 1))) == 2/sqrt(21)
#' 
#' @seealso \code{\link{normalizedLermanIndex}}
#' @export
setGeneric("lermanIndex", function(p, q, c=NULL) standardGeneric("lermanIndex"))

#' @describeIn lermanIndex Compute given two partitions
setMethod("lermanIndex", signature(p="Partition", q="Partition", c="missing"),
          function(p, q, c=NULL) lermanIndex(p, q, computePairCoefficients(p, q)))

#' @describeIn lermanIndex Compute given the partitions and pair coefficients
setMethod("lermanIndex", signature(p="Partition", q="Partition", c="PairCoefficients"),
          function(p, q, c) {
            cluster_sizes_p <- table(p)
            cluster_sizes_q <- table(q)
            n <- length(p)
            v1 <- function(x) sum(sapply(x, function(i) i * (i - 1)))
            v2 <- function(x) sum(sapply(x, function(i) i * (i - 1) * (i - 2)))
            v3 <- function(x) v1(x)**2 - 2 * sum(sapply(x, function(i) i * (i - 1) * (2 * i - 3)))
            
            expectation_N11 <- sum(sapply(cluster_sizes_p, function(i) choose(i, 2))) * 
              sum(sapply(cluster_sizes_q, function(j) choose(j, 2))) / N(c)
            variance_N11 <- v1(cluster_sizes_p) * v1(cluster_sizes_q) / (2 * n * (n - 1)) +
              v2(cluster_sizes_p) * v2(cluster_sizes_q) / (n * (n - 1) * (n - 2)) +
              v3(cluster_sizes_p) * v3(cluster_sizes_q) / (4 * n * (n - 1) * (n - 2) * (n - 3)) -
              (v1(cluster_sizes_p) * v1(cluster_sizes_q) / (2 * n * (n - 1)))**2
            
            (N11(c) - expectation_N11) / sqrt(variance_N11)
          })


#' Normalized Lerman Index
#' 
#' Compute the normalized Lerman index
#' \deqn{L(P, Q) / \sqrt{L(P, P)L(Q, Q)}}
#' where \eqn{L} is the Lerman index.
#' 
#' @references
#' \insertRef{Lerman1988}{partitionComparison}
#' 
#' \insertRef{Hubert1985}{partitionComparison}
#' 
#' @param p The partition \eqn{P}
#' @param q The partition \eqn{Q}
#' @param c \linkS4class{PairCoefficients} or NULL
#' @template author
#' @name normalizedLermanIndex
#' @examples 
#' normalizedLermanIndex(new("Partition", c(0, 0, 0, 1, 1)), 
#'                       new("Partition", c(0, 0, 1, 1, 1))) == 1/6
#' 
#' @seealso \code{\link{lermanIndex}}
#' @export
setGeneric("normalizedLermanIndex", 
           function(p, q, c=NULL) standardGeneric("normalizedLermanIndex"))

#' @describeIn normalizedLermanIndex Compute given two partitions
setMethod("normalizedLermanIndex", signature(p="Partition", q="Partition", c="missing"),
          function(p, q, c=NULL) normalizedLermanIndex(p, q, computePairCoefficients(p, q)))

#' @describeIn normalizedLermanIndex Compute given the partitions and pair coefficients
setMethod("normalizedLermanIndex", 
          signature(p="Partition", q="Partition", c="PairCoefficients"),
          function(p, q, c) {
            lp <- lermanIndex(p, p, computePairCoefficients(p, p))
            lq <- lermanIndex(q, q, computePairCoefficients(q, q))
                             
            lermanIndex(p, q, c) / sqrt(lp * lq)
          })
