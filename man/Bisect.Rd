\name{Bisect}
\alias{Bisect}

\title{An implementation of the bisection algorithm for root finding. }
\description{
  Most of the optimizations in \emph{Icens} have a one dimensional
  root-finding component. Since the quantities involved are generally
  restricted to a subset of [0,1] we use bisection to find the roots.
}
\usage{
Bisect(tA, pvec, ndir, Meps, tolbis=1e-07)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{tA}{ The transpose of the clique matrix. }
  \item{pvec}{ The current estimate of the probability vector. }
  \item{ndir}{ The direction to explore. }
  \item{Meps}{ Machine epsilon, elements of \code{pvec} that are less
    than this are assumed to be zero. }
  \item{tolbis}{ The tolerance used to determine if the algorithm has
    converged.  }
}
\details{
 We search from \code{pvec} in the direction \code{ndir} to obtain the
 new value of \code{pvec} that maximizes the likelihood.
}
\value{
  The new estimate of \code{pvec}.
}
\references{ Any book on optimization. }
\author{ Alain Vandal and Robert Gentleman. }


\keyword{ optimize }
