\name{BVcliques}
\alias{BVcliques}
\title{ Find the bivariate cliques from the marginal data. }
\description{
 The maximal cliques of the intersection graph are obtained by first
 finding the cliques for the marginal data and then combining them using
 the algorithm in Gentleman and Vandal (1999).
}
\usage{
BVcliques(intvlx, intvly, Lxopen=TRUE, Rxopen=FALSE,
                    Lyopen=TRUE, Ryopen=FALSE )
}
\arguments{
  \item{intvlx}{ The cliques for one marginal component, alternatively
    the marginal intervals can be supplied. }
  \item{intvly}{ The cliques for the other marginal component, alternatively
    the marginal intervals can be supplied. }
  \item{Lxopen}{ Boolean indicating whether the left end point in the x
    coordinate is open.}
  \item{Rxopen}{ Boolean indicating whether the right end point in the x
      coordinate is open.}
  \item{Lyopen}{ Boolean indicating whether the left end point in the y
      coordinate is open.}
  \item{Ryopen}{ Boolean indicating whether the right end point in the y
    coordinate is open.}
}

\value{
  A list of the maximal cliques of the intersection graph of the data.
}
\references{ \emph{Graph--Theoretical Aspects of Bivariate Censored
    Data}, R. Gentleman and A. Vandal, 1999, submitted. }
\author{ A. Vandal and R. Gentleman }

\seealso{  \code{\link{BVclmat}}, \code{\link{BVsupport}} }

\examples{
   data(cmv)
   cmv.cl <- BVcliques(cmv[,1:2], cmv[,3:4], Lxopen=FALSE, Lyopen=FALSE )
}
\keyword{ manip}


