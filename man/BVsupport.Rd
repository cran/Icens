\name{BVsupport}
\alias{BVsupport}
\title{ Compute the support for the cliques of a bivariate intersection
  graph.  }
\description{
  Given the regions where the events occurred and the cliques of the
  intersection graph the support of the cliques is computed. For each
  clique it is the intersection of the event time regions for all
  observations in that clique.
}
\usage{
BVsupport(intvlx, intvly, cliques=BVcliques(intvlx, intvly))
}
\arguments{
  \item{intvlx}{ The event time intervals for one dimension. }
  \item{intvly}{ The event time intervals for the other dimension. }
  \item{cliques}{ The list of maximal cliques of the intersection graph,
  optionally. }
}
\value{
  An m by 4 matrix containing the corners of the intervals of support
  for the maximal cliques of the intersection graph corresponding to the
  first two arguments to the function.
}
\references{\emph{Graph--Theoretical Aspects of Bivariate Censored
    Data}, R. Gentleman and A. Vandal, 1999, submitted. }
\author{ A. Vandal and R. Gentleman }

\seealso{ \code{\link{BVcliques}}, \code{\link{BVclmat}} }

\examples{
  data(cmv)
  cmv.cl <- BVcliques(cmv[,1:2], cmv[,3:4])
  boxes <- BVsupport(cmv[,1:2], cmv[,3:4], cmv.cl)
}
\keyword{manip}
