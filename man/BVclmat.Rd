\name{BVclmat}
\alias{BVclmat}
\title{ Comput the clique matrix from the clique list. }
\description{
  Given the clique list, obtained from \code{\link{BVcliques}}, the
  clique matrix is obtained. This is the m (number of cliques) by n
  (number of observations) matrix. A[i,j] is one if individual j is in
  maximal clique i.
}
\usage{
BVclmat(cliques)
}

\arguments{
  \item{cliques}{ The clique list. }
}

\value{
  The m by n clique matrix.
}
\references{  \emph{Graph--Theoretical Aspects of Bivariate Censored
    Data}, R. Gentleman and A. Vandal, 1999, submitted. }
\author{ A. Vandal and R. Gentleman }

\seealso{  \code{\link{BVcliques}}, \code{\link{BVsupport}} }

\examples{
   data(cmv)
   bcl <- BVcliques(cmv[,1:2], cmv[,3:4])
   A <- BVclmat(bcl)
}
\keyword{manip}
