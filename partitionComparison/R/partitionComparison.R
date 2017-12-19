#' @details 
#' This package provides a large collection of measures to compare two 
#' partitions. Some survey articles for these measures are cited below, 
#' the seminal papers for each individual measure is provided with the 
#' function definition.
#' 
#' Most functionality is implemented as S4 classes and methods so that an 
#' adoption is easily possible for special needs and specifications. 
#' The main class is \code{\linkS4class{Partition}} which merely wraps an atomic 
#' vector of length \eqn{n} for storing the class label of each object.
#' The computation of all measures is designed to work on vectors
#' of class labels. 
#' 
#' All partition comparison methods can be called in the 
#' same way: \code{<measure method>(p, q)} with \code{p, q} being the two 
#' partitions (as \code{\linkS4class{Partition}} instances).
#' One often does not explicitly want to transform the vector of class labels
#' (as output of another package's function/algorithm) into 
#' \code{\linkS4class{Partition}} instances before using measures from this 
#' package. For convenience, the function 
#' \code{\link{registerPartitionVectorSignatures}} exists which dynamically creates
#' versions of all measures that will directly work with plain R vectors.
#' 
#' @examples 
#' # Generate some data
#' set.seed(42)
#' data <- cbind(x=c(rnorm(50), rnorm(30, mean=5)), y=c(rnorm(50), rnorm(30, mean=5)))
#' # Run k-means with two/three centers
#' data.km2 <- kmeans(data, 2)
#' data.km3 <- kmeans(data, 3)
#' 
#' # Load this library
#' library(partitionComparison)
#' # Register the measures to take ANY input
#' registerPartitionVectorSignatures(environment())
#' # Compare the clusters
#' randIndex(data.km2$cluster, data.km3$cluster)
#' # [1] 0.8101266
#' 
#' @references 
#' \insertRef{Albatineh2006}{partitionComparison}
#' 
#' \insertRef{Meila2007}{partitionComparison}
#'
#' @importFrom Rdpack reprompt
"_PACKAGE"
