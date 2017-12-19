
#' Mutual Information
#' 
#' Compute the mutual information
#' \deqn{
#' \sum_{C \in P} \sum_{D \in Q} {\frac{|C \cap D|}{n} \log n\frac{|C \cap D|}{|C| |D|}}
#' }
#' 
#' @references
#' \insertRef{Vinh2010}{partitionComparison}
#'
#' @template params
#' @template author
#' @name mutualInformation
#' @examples 
#' isTRUE(all.equal(mutualInformation(new("Partition", c(0, 0, 0, 1, 1)), 
#'                  new("Partition", c(0, 0, 1, 1, 1))), 4/5*log(5/3) + 1/5*log(5/9)))
#' 
#' @seealso \code{\link{normalizedMutualInformation}}
#' @export
setGeneric("mutualInformation", function(p, q) standardGeneric("mutualInformation"))

#' @describeIn mutualInformation Compute given two partitions
setMethod("mutualInformation", signature(p="Partition", q="Partition"),
          function(p, q) {
            clusters_p = unique(p)
            clusters_q = unique(q)
            # Contingency table => cluster sizes
            c_sizes_p <- table(p)
            c_sizes_q <- table(q)
            n <- length(p)
            
            sum(sapply(clusters_p, 
                       function(i) sum(sapply(clusters_q, 
                                              function(j) {
                                                if (setOverlap(p, q, i, j) > 0) {
                                                  setOverlap(p, q, i, j) / n * 
                                                    log(n * setOverlap(p, q, i, j) / 
                                                          (c_sizes_p[[as.character(i)]] * 
                                                             c_sizes_q[[as.character(j)]]))
                                                } else {
                                                  0
                                                }
                                              })
                                       )
                       )
                )
          })


#' Normalized Mutual Information
#' 
#' Compute the mutual information (\eqn{MI}) which is normalized either by the
#' minimum/maximum partition entropy (\eqn{H})
#' \deqn{\frac{MI(P, Q)}{\varphi(H(P), H(Q))},\ \varphi \in \{\min, \max\}}
#' or the sum
#' \deqn{\frac{2 \cdot MI(P, Q)}{H(P) + H(Q)}}
#' 
#' @references
#' \insertRef{Kvalseth1987}{partitionComparison}
#'
#' @template params
#' @param type One of "min" (default), "max" or "sum"
#' @template author
#' @name normalizedMutualInformation
#' @examples
#' isTRUE(all.equal(normalizedMutualInformation(
#'                    new("Partition", c(0, 0, 0, 1, 1)),
#'                    new("Partition", c(0, 0, 1, 1, 1)), "min"),
#'                  normalizedMutualInformation(
#'                    new("Partition", c(0, 0, 0, 1, 1)), 
#'                    new("Partition", c(0, 0, 1, 1, 1)), "max")
#'                  ))
#' 
#' @seealso \code{\link{mutualInformation}}, \code{\link{entropy}}
#' @export
setGeneric("normalizedMutualInformation", function(p, q, type=c("min", "max", "sum")) 
  standardGeneric("normalizedMutualInformation"))

#' @describeIn normalizedMutualInformation Compute given two partitions
setMethod("normalizedMutualInformation", 
          signature(p="Partition", q="Partition", type="character"),
          function(p, q, type) {
            res<-switch(type,
                        min = mutualInformation(p, q) / min(entropy(p), entropy(q)),
                        max = mutualInformation(p, q) / max(entropy(p), entropy(q)),
                        sum = 2 * mutualInformation(p, q) / sum(entropy(p), entropy(q))
                        )
            
            if (is.null(res))
              stop(paste("Not a valid type: ", type))
            
            res
          })

#' @describeIn normalizedMutualInformation Compute given two partitions with \code{type="min"}
setMethod("normalizedMutualInformation", 
          signature(p="Partition", q="Partition", type="missing"),
          function(p, q, type=NULL) normalizedMutualInformation(p, q, "min"))


#' Variation of Information
#' 
#' Compute the variation of information 
#' \deqn{H(P) + H(Q) - 2MI(P, Q)}
#' where \eqn{MI} is the mutual information, \eqn{H} the partition entropy
#' 
#' @references
#' \insertRef{Meila2003}{partitionComparison}
#' 
#' \insertRef{Meila2007}{partitionComparison}
#'
#' @template params
#' @template author
#' @name variationOfInformation
#' @examples 
#' isTRUE(all.equal(variationOfInformation(new("Partition", c(0, 0, 0, 1, 1)),
#'                                         new("Partition", c(0, 0, 1, 1, 1))),
#'                                         0.763817))
#' 
#' @seealso \code{\link{mutualInformation}}, \code{\link{entropy}}
#' @export
setGeneric("variationOfInformation", 
           function(p, q) standardGeneric("variationOfInformation"))

#' @describeIn variationOfInformation Compute given two partitions
setMethod("variationOfInformation", signature(p="Partition", q="Partition"),
          function(p, q) {
            entropy(p) + entropy(q) - 2 * mutualInformation(p, q)
          })
