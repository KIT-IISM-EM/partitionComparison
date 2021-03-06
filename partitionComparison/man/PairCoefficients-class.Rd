% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa-classes.R
\docType{class}
\name{PairCoefficients-class}
\alias{PairCoefficients-class}
\alias{PairCoefficients}
\title{S4 class to represent coefficients of object pairs for the comparison of two
object partitions (say \eqn{P} and \eqn{Q}).}
\description{
S4 class to represent coefficients of object pairs for the comparison of two
object partitions (say \eqn{P} and \eqn{Q}).
}
\section{Slots}{

\describe{
\item{\code{N11}}{The number of object pairs that are in both partitions together in a cluster}

\item{\code{N00}}{The number of object pairs that are in no partition together in a cluster}

\item{\code{N10}}{The number of object pairs that are only in partition \eqn{P} together in a cluster}

\item{\code{N01}}{The number of object pairs that are only in partition \eqn{Q} together in a cluster}
}}

\seealso{
\code{\link{N11}} \code{\link{N10}} \code{\link{N01}} \code{\link{N00}}
}
\author{
Fabian Ball \email{fabian.ball@kit.edu}
}
