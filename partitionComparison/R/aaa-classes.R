
#' Simple S4 class to represent a partition of objects as vector of class labels.
#' 
#' This class is a wrapper around a vector but allows only the atomic vectors
#' logical, numeric, integer, complex, character, raw.
#' The reason for this is that only those types seem to make sense as class labels.
#' Furthermore, class labels are immutable.
#' 
#' @examples 
#' p <- new("Partition", c(0, 0, 1, 1, 1))
#' q <- new("Partition", c("a", "a", "b", "b", "b"))
#' 
#' \dontrun{
#' # This won't work:
#' new("Partition", c(list("a"), "a", "b", "b", "b"))
#' p[2] <- 2
#' }
#' 
#' @import methods
#' 
#' @template author
#' @export
Partition <- setClass("Partition", contains = "vector", 
                      validity = function(object) {
                        types <- c("logical", "numeric", "integer", "complex", "character", "raw")
                        if (!(mode(object) %in% types))
                          stop(c("Only \"", paste(types, collapse=", "), "\" are allowed"))
                      })

#' Subsetting \linkS4class{Partition} instances
#' 
#' This method overrides the standard subsetting to prevent
#' alteration (makes partitions, i.e. class labels, immutable).
#' 
#' @param x A \linkS4class{Partition} instance
#' @param i \code{\link[base]{Extract}}
#' @param j \code{\link[base]{Extract}}
#' @param value \code{\link[base]{Extract}}
#' 
#' @template author
#' @aliases [<-,Partition-method
#' @docType methods
#' @export
setMethod("[<-", signature(x="Partition"), 
          function(x, i, j, value) stop("\"Partition\" is immutable"))


#' S4 class to represent coefficients of object pairs for the comparison of two
#' object partitions (say \eqn{P} and \eqn{Q}).
#' 
#' @slot N11 The number of object pairs that are in both partitions together in a cluster
#' @slot N00 The number of object pairs that are in no partition together in a cluster
#' @slot N10 The number of object pairs that are only in partition \eqn{P} together in a cluster
#' @slot N01 The number of object pairs that are only in partition \eqn{Q} together in a cluster
#' 
#' @seealso \code{\link{N11}} \code{\link{N10}} \code{\link{N01}} \code{\link{N00}}
#' 
#' @import methods
#' 
#' @template author
#' @export
PairCoefficients <- setClass("PairCoefficients", 
                             slots = c(N11 = "numeric", N10 = "numeric", 
                                       N01 = "numeric", N00 = "numeric"),
                             validity = function(object) {
                               if (any(object@N11 < 0, object@N10 < 0, 
                                       object@N01 < 0, object@N00 < 0))
                                 stop("Coefficients must be >= 0!")
                             })

#' Method to retrieve the coefficient \eqn{N_{11}}
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N11
#' @export
setGeneric("N11", function(obj) { standardGeneric("N11") })

#' @rdname N11
setMethod("N11", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N11)
          })

#' Method to retrieve the coefficient \eqn{N_{01}}
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N01
#' @export
setGeneric("N01", function(obj) { standardGeneric("N01") })

#' @rdname N01
setMethod("N01", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N01)
          })

#' Method to retrieve the coefficient \eqn{N_{10}}
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N10
#' @export
setGeneric("N10", function(obj) { standardGeneric("N10") })

#' @rdname N10
setMethod("N10", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N10)
          })

#' Method to retrieve the coefficient \eqn{N_{00}}
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N00
#' @export
setGeneric("N00", function(obj) { standardGeneric("N00") })

#' @rdname N00
setMethod("N00", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N00)
          })

#' Method to retrieve the complex coefficient \eqn{N_{21}}
#' 
#' It is defined as \eqn{N_{21} = N_{11} + N_{10}}
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N21
#' @export
setGeneric("N21", function(obj) { standardGeneric("N21") })

#' @rdname N21
setMethod("N21", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N11 + obj@N10)
          })

#' Method to retrieve the complex coefficient \eqn{N_{12}}
#' 
#' It is defined as \eqn{N_{12} = N_{11} + N_{01}}
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N12
#' @export
setGeneric("N12", function(obj) { standardGeneric("N12") })

#' @rdname N12
setMethod("N12", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N11 + obj@N01)
          })

#' Method to retrieve the complex coefficient \eqn{N'_{01}}
#' 
#' It is defined as \eqn{N'_{01} = N_{00} + N_{01}}
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N01p
#' @export
setGeneric("N01p", function(obj) { standardGeneric("N01p") })

#' @rdname N01p
setMethod("N01p", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N00 + obj@N01)
          })

#' Method to retrieve the complex coefficient \eqn{N'_{10}}
#' 
#' It is defined as \eqn{N'_{10} = N_{00} + N_{10}}
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N10p
#' @export
setGeneric("N10p", function(obj) { standardGeneric("N10p") })

#' @rdname N10p
setMethod("N10p", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N00 + obj@N10)
          })

#' Method to retrieve the complex coefficient \eqn{N}
#' 
#' It is defined as \eqn{N = N_{11} + N_{10} + N_{01} + N_{00}} which equals 
#' \eqn{n \choose{2}} with \eqn{n} the number of objects
#' 
#' @param obj Instance of \linkS4class{PairCoefficients}
#' 
#' @template author
#' @name N
#' @export
setGeneric("N", function(obj) { standardGeneric("N") })

#' @rdname N
setMethod("N", signature = "PairCoefficients", 
          definition = function(obj) {
            return(obj@N11 + obj@N01 + obj@N10 + obj@N00)
          })
