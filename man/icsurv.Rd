\name{icsurv}
\alias{icsurv}
\title{ The class of objects returned by the estimation routines in the
  Icens library. }
\description{
  An object of class \code{icsurv} must contain the following
  components:
  \item{converge}{A boolean indicating whether the iteration producing
    \code{pf} converged.} 
  \item{pf }{The probability vector.}
  

  It can optionally contain any of the following components:
  \item{clmat }{The clique matrix used to obtain \code{pf}.}
  \item{intmap}{The real representations of the support for the
    components of \code{pf}.}
  \item{iter }{A matrix containing every iterative estimate of
    \code{pf}, useful for debugging.}
  \item{lval}{The value of the \bold{log} likelihood at \code{pf}.}
  \item{numiter}{The number of iterations taken.} 
  \item{sigma}{The cumulative sum of \code{pf}.}
  \item{weights}{Weights used in the EMICM algorithm.}
}

\author{ Alain Vandal and Robert Gentleman. }

\seealso{ \code{\link{VEM}}, \code{\link{ISDM}}, \code{\link{EMICM}},
     \code{\link{PGM}}, \code{\link{EM}} }
\keyword{methods}
