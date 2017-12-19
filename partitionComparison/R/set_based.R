
#' Dongen's Metric
#' 
#' Compute Dongen's metric
#' \deqn{
#' 2n - \sum_{C \in P} \max_{D \in Q} |C \cap D| - \sum_{D \in Q} \max_{C \in P} |C \cap D|
#' }
#' 
#' @references
#' \insertRef{vanDongen2000}{partitionComparison}
#'
#' @template params
#' @template author
#' @name dongensMetric
#' @examples 
#' isTRUE(all.equal(dongensMetric(new("Partition", c(0, 0, 0, 1, 1)), 
#'                                new("Partition", c(0, 0, 1, 1, 1))), 2))
#'               
#' @seealso \code{\link{projectionNumber}}
#' 
#' @export
setGeneric("dongensMetric", function(p, q) standardGeneric("dongensMetric"))

#' @describeIn dongensMetric Compute given two partitions
setMethod("dongensMetric", signature(p="Partition", q="Partition"),
          function(p, q) 2 * length(p) - projectionNumber(p, q) - projectionNumber(q, p))


#' Larsen & Aone Measure
#' 
#' Compute the measure of Larsen and Aone
#' \deqn{
#' \frac{1}{|\cal{P}|} 
#' \sum_{C \in \cal{P}}{\max_{D \in \cal{Q}}{\frac{2|C \cap D|}{|C| + |D|}}}
#' }
#' 
#' @references
#' \insertRef{Larsen1999}{partitionComparison}
#'
#' @template params
#' @template author
#' @name larsenAone
#' @examples 
#' isTRUE(all.equal(larsenAone(new("Partition", c(0, 0, 0, 1, 1)), 
#'                             new("Partition", c(0, 0, 1, 1, 1))), 0.8))
#' 
#' @export
setGeneric("larsenAone", function(p, q) standardGeneric("larsenAone"))

#' @describeIn larsenAone Compute given two partitions
setMethod("larsenAone", signature(p="Partition", q="Partition"),
          function(p, q) {
            clusters_p <- unique(p)
            clusters_q <- unique(q)
            # Contingency table => cluster sizes
            c_sizes_p <- table(p)
            c_sizes_q <- table(q)
            
            sum(sapply(clusters_p, 
                       function(i) max(sapply(clusters_q, 
                                              function(j) 2 * setOverlap(p, q, i, j) / 
                                                (c_sizes_p[[as.character(i)]] + 
                                                   c_sizes_q[[as.character(j)]])
                                              )
                                       )
                       )
                ) / length(clusters_p)
          })


#' Classification Error Distance
#' 
#' Compute the classification error distance
#' \deqn{1 - \frac{1}{n} \max_{\sigma}{\sum_{C \in \cal{P}}{|C \cap \sigma(C)|}}}
#' with \eqn{\sigma} a weighted matching between the clusters of both partitions.
#' The nodes are the classes of each partition, the weights are the overlap of objects.
#' 
#' @section Hint:
#' This measure is implemented using \code{\link[lpSolve]{lp.assign}} from
#' the \code{lpSolve} package to compute the maxmimal matching of a 
#' weighted bipartite graph.
#' 
#' @references
#' \insertRef{Meila2001}{partitionComparison}
#' 
#' \insertRef{Meila2005}{partitionComparison}
#'
#' @template params
#' @template author
#' @name classificationErrorDistance
#' @examples 
#' isTRUE(all.equal(classificationErrorDistance(new("Partition", c(0, 0, 0, 1, 1)), 
#'                                              new("Partition", c(0, 0, 1, 1, 1))), 0.2))
#' 
#' @import lpSolve
#' @export
setGeneric("classificationErrorDistance", 
           function(p, q) standardGeneric("classificationErrorDistance"))

#' @describeIn classificationErrorDistance Compute given two partitions
setMethod("classificationErrorDistance", signature(p="Partition", q="Partition"),
          function(p, q) {
            clusters_p <- unique(p)
            clusters_q <- unique(q)
            # This creates a matrix with the classes of p on the rows,
            # the classes of q on the columns.
            # Entries are the overlap.
            w <- sapply(clusters_q, 
                        function(j) sapply(clusters_p, 
                                           function(i) setOverlap(p, q, i, j )))
            
            add_rows <- length(clusters_q) - length(clusters_p)
            # Fill matrix with rows of zeros to get an quadratic matrix 
            # (needed for lp.assign)
            if (add_rows > 0)
              w <- rbind(w, matrix(rep(0, add_rows * length(clusters_q)), 
                                   nrow = add_rows, ncol = length(clusters_q)))
            
            stopifnot(nrow(w) == ncol(w))
            
            matching <- lpSolve::lp.assign(w, "max")
            
            1 - matching$objval / length(p)
          })
