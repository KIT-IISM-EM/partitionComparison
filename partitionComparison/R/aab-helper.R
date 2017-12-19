
#' Compute the length of the overlap of two sets
#' 
#' Given two partitions (p, q) represented as vectors of cluster ids,
#' compute the overlap between clusters i and j.
#' 
#' @param p Partition \eqn{P}
#' @param q Partition \eqn{Q}
#' @param i Cluster id of \eqn{P}
#' @param j Cluster id of \eqn{Q}
#' 
#' @template author
#' 
#' @noRd
setOverlap <- function(p, q, i, j) sum(p == i & q == j)


#' Compute the projection number of two partitions
#' 
#' Given two partitions (p, q) represented as vectors of cluster ids,
#' compute the projection number which is the sum of maximum
#' cluster overlaps for all clusters of \eqn{P} to any cluster of \eqn{Q}.
#' 
#' @param p Partition \eqn{P}
#' @param q Partition \eqn{Q}
#' 
#' @template author
#' @examples 
#' isTRUE(all.equal(projectionNumber(c(0, 0, 0, 1, 1), c(0, 0, 1, 1, 1)), 4))
#' 
#' @seealso \code{\link{dongensMetric}}
#' 
#' @export
projectionNumber <- function(p, q) {
  # Gets the (arbitrary) cluster ids
  clusters_p <- unique(p)
  clusters_q <- unique(q)
  
  # Computes the sum over the maximum overlaps
  sum(sapply(clusters_p, 
             function(i) max(sapply(clusters_q, 
                                    function(j) setOverlap(p, q, i, j)))))
}


#' Entropy
#' 
#' Compute the Shannon entropy
#' \deqn{-\sum_{i} p_i \log_b p_i}
#' 
#' @section Hint:
#' This method is used internally for measures based on information theory
#' 
#' @param x A probability distribution
#' @param log_base Optional base of the logarithm (default: \eqn{e})
#' @template author
#' @name entropy
#' @examples 
#' isTRUE(all.equal(entropy(c(.5, .5)), log(2)))
#' isTRUE(all.equal(entropy(c(.5, .5), 2), 1))
#' isTRUE(all.equal(entropy(c(.5, .5), 4), .5))
#' 
#' @export
setGeneric("entropy", function(x, log_base) standardGeneric("entropy"))

#' @rdname entropy
setMethod("entropy", signature(x="numeric", log_base="numeric"), 
          function(x, log_base) {
            if (!isTRUE(all.equal(sum(x), 1)))
              stop("x is not a distribution that sums to 1")
            
            -sum(sapply(x[x > 0], function(p) p * log(p, base = log_base)))
          })

#' @examples 
#' # Entropy of a partition
#' isTRUE(all.equal(entropy(new("Partition", c(0, 0, 1, 1, 1))), entropy(c(2/5, 3/5))))
#' 
#' @describeIn entropy Entropy of a partition represented by \code{x}
setMethod("entropy", signature(x="Partition", log_base="numeric"), 
          function(x, log_base) entropy(as.vector(table(x) / length(x)), log_base))

#' @rdname entropy
setMethod("entropy", signature(x="ANY", log_base="missing"), 
          function(x, log_base=exp(1)) entropy(x, log_base=exp(1)))


# Internal function to create a 'proxy' function for the given
# function. The result is a new function with exactly the same
# signature but the first two arguments are replaced by instances
# of Partition.
.createProxyFunction <- function(f) {
  f_formals <- formals(f)  # Get the formals for the method,...
  formal_names <- names(f_formals)  # ...and its names.
  # Create a list of symbols from the formals (to be used in 'do.call')
  formal_symbol_list <- lapply(formal_names, as.symbol)
  
  # Create the new proxy function
  if (length(formal_symbol_list) > 2) {
    func <- function() {
      # Cast to partition
      x <- new("Partition", get(formal_names[1]))
      y <- new("Partition", get(formal_names[2]))
      do.call(f, c(list(x, y), formal_symbol_list[3:length(formal_symbol_list)]))
    }
  } else {
    func <- function() {
      # Cast to partition
      x <- new("Partition", get(formal_names[1]))
      y <- new("Partition", get(formal_names[2]))
      f(x, y)
    }
  }
  
  # Finally, set the 'original' formals to the proxy function
  formals(func) <- f_formals
  
  func
}


#' Make comparison measures usable with any vectors
#' 
#' The comparison measures are defined to use the class \linkS4class{Partition} 
#' as parameters. If you do not want to explicitly convert an arbitrary vector 
#' of class labels (probably as a result from another package's algorithm) into
#' a \linkS4class{Partition} instance, calling this function will create methods
#' for all measures that allow "ANY" input which is implicitly converted to
#' \linkS4class{Partition}.
#' 
#' @param e The environment to register the methods in 
#' (mostly \code{environment()} is fine)
#' 
#' @template author
#' @examples 
#' library(partitionComparison)
#' randIndex(new("Partition", c(0, 0, 0, 1, 1)), new("Partition", c(0, 0, 1, 1, 1)))
#' # [1] 0.6
#' \dontrun{randIndex(c(0, 0, 0, 1, 1), c(0, 0, 1, 1, 1))}
#' # Error in (function (classes, fdef, mtable) :
#' # unable to find an inherited method for function 'randIndex' for signature '"numeric", "numeric"'
#' registerPartitionVectorSignatures(environment())
#' randIndex(c(0, 0, 0, 1, 1), c(0, 0, 1, 1, 1))
#' # [1] 0.6
#' 
#' @export
registerPartitionVectorSignatures <- function(e) {
  # Get all S4 method names defined for "Partition"
  method_info <- attr(.S4methods(class="Partition"), "info")
  partition_methods <- unique(method_info[method_info$from == "partitionComparison",]$generic)

  for (method_name in partition_methods) {
    # Get all methods that have a Partition class in their signature and the given name
    methods_to_overload <- findMethods(method_name, classes = "Partition")

    for (method_definition in methods_to_overload) {
      method_signature <- method_definition@target
      
      # If the first two parameters are of type "Partition", we assume it is
      # a measure and register a method with an altered signature 
      # ("Partition" is replaced by "vector") and wrap the call to the 'original'
      # method.
      if (length(method_signature) >= 2 & 
          method_signature[1] == "Partition" & 
          method_signature[2] == "Partition") {
        method_signature[1] <- method_signature[2] <- "vector"

        setMethod(method_name, method_signature, 
                  .createProxyFunction(method_definition), where = e)
      }
    }
  }
}


#' Compare two partitions with all measures
#' 
#' Compute the comparison between two partitions for all available measures.
#' 
#' 
#' @param p The partition \eqn{P}
#' @param q The partition \eqn{Q}
#' 
#' @return 
#' Instance of \code{\link{data.frame}} with columns \code{measure} and \code{value}
#' 
#' @section Warning:
#' This method will identify every generic S4 method that has a signature 
#' \code{"Partition", "Partition"} (including signatures with following \code{"missing"} 
#' parameters, e.g. \code{"Partition", "Partition", "missing"}) as a partition 
#' comparison measure, \strong{except} this method itself (otherwise: infinite 
#' recursion). This means one has to take care when defining other methods with the same
#' signature in order not to produce unwanted side-effects!
#' 
#' @template author
#' @name compareAll
#' 
#' @examples 
#' compareAll(new("Partition", c(0, 0, 0, 1, 1)), new("Partition", c(0, 0, 1, 1, 1)))
#' \dontrun{
#'                         measure       value
#'  1            adjustedRandIndex 0.166666667
#'  2                     baulieu1 0.760000000
#'  3                     baulieu2 0.040000000
#'  4  classificationErrorDistance 0.200000000
#'  5                  czekanowski 0.500000000
#'  6                dongensMetric 2.000000000
#'  7                 fagerMcGowan 0.250000000
#'  8          folwkesMallowsIndex 0.500000000
#'  9              gammaStatistics 0.166666667
#'  10              goodmanKruskal 0.333333333
#'  11               gowerLegendre 0.750000000
#'  12                      hamann 0.200000000
#'  13          jaccardCoefficient 0.333333333
#'  14                  kulczynski 0.500000000
#'  15                  larsenAone 0.800000000
#'  16                 lermanIndex 0.436435780
#'  17                mcconnaughey 0.000000000
#'  18            minkowskiMeasure 1.000000000
#'  19                mirkinMetric 8.000000000
#'  20           mutualInformation 0.291103166
#'  21       normalizedLermanIndex 0.166666667
#'  22 normalizedMutualInformation 0.432538068
#'  23                     pearson 0.006944444
#'  24                      peirce 0.166666667
#'  25                   randIndex 0.600000000
#'  26              rogersTanimoto 0.428571429
#'  27                   russelRao 0.200000000
#'  28               rvCoefficient 0.692307692
#'  29                sokalSneath1 0.583333333
#'  30                sokalSneath2 0.200000000
#'  31                sokalSneath3 0.333333333
#'  32      variationOfInformation 0.763817002
#'  33                    wallaceI 0.500000000
#'  34                   wallaceII 0.500000000
#' }
#' 
#' @export
setGeneric("compareAll", function(p, q) standardGeneric("compareAll"))

#' @describeIn compareAll Compare given two \code{\linkS4class{Partition}} instances
setMethod("compareAll", signature(p="Partition", q="Partition"), 
          function(p, q) {
            result.df <- data.frame(measure=character(), value=numeric(), stringsAsFactors = FALSE)
  
            method_info <- attr(.S4methods(class="Partition"), "info")
            partition_methods <- unique(method_info[method_info$from == "partitionComparison",]$generic)
            
            for (method_name in partition_methods) {
              if (method_name == "compareAll") next
              # Get all methods that have a Partition class in their signature and 
              # the given name
              methods_to_overload <- findMethods(method_name, classes = "Partition")
              
              for (method_definition in methods_to_overload) {
                method_signature <- method_definition@target
                
                # If the first two parameters are of type "Partition", we assume it is
                # a measure and register a method with an altered signature 
                # ("Partition" is replaced by "vector") and wrap the call to the 'original'
                # method.
                if ((length(method_signature) == 2) & 
                    method_signature[1] == "Partition" & 
                    method_signature[2] == "Partition") {
                  result.df[nrow(result.df)+1,] <- list(method_name, 
                                                        do.call(method_name, list(p, q)))
                } else if ((length(method_signature) > 2) & 
                           method_signature[1] == "Partition" & 
                           method_signature[2] == "Partition" &
                           all(method_signature[3:length(method_signature)] == "missing")) {
                  result.df[nrow(result.df)+1,] <- list(method_name, 
                                                        do.call(method_name, list(p, q)))
                }
              }
            }
            
            result.df
          })
