\name{PMGA}
\alias{PMGA}
\title{ Implement the pool monotone groups algorithm. }
\description{
  For isotonization problems some increase in speed and decrease in
  complexity can be achieved through the use of the pool monotone groups
  algorithm of Y.L. Zhang and M.A. Newton (1997).
  It isotonizes a weighted and ordered set of values.
}
\usage{
PMGA(est, ww=rep(1, length(est)))
}
\arguments{
  \item{est}{ The vector of values, in the appropriate order. }
  \item{ww}{ The weight vector. }
}
\details{
  To be supplied at some later date.
}
\value{
 An object containing the following components:
  \item{est }{ The isotonized estimates. }
  \item{ww }{The weights associated with the isotonized estimates. }
  \item{poolnum }{ The number of values pooled in the current
    estimate. }
  \item{passes }{ The number of passes which were required to isotonize
    the list. } 
}
\references{ Y.L. Zhang and M.A. Newton (1997),
    http://www.stat.wisc.edu/~newton/newton.html) }
\author{ Alain Vandal and Robert Gentleman. }

\seealso{ \code{\link{EMICM}} }

\keyword{optimize}
