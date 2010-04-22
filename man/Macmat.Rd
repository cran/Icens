\name{Macmat}
\alias{Macmat}
\title{A function to compute the incidence matrix for an intersection graph. }
\description{
  Returns the Petrie matrix and Petrie pairs of an interval order given its
  list of maximal antichains. These can be obtained from
  \code{\link{Maclist}}. 
}
\usage{
Macmat(ml)
}
\arguments{
  \item{ml}{ A list containing the maximal cliques of the intersection
    graph of the data. }
}
\details{
  Not worth mentioning?
}
\value{
  A list containing two components.  
  \item{pmat }{The Petrie or clique matrix of the underlying interval order.}
  \item{ppairs }{The Petrie pairs for each observation. These indicate
    the first and last maximal clique occupied by the observation.}
}
\references{Computational Methods for Censored Data using Intersection
  Graphs, R. Gentleman and A. Vandal, JCGS, 2000.  }
\author{ Alain Vandal and Robert Gentleman}

\seealso{ \code{\link{Maclist}} }

\examples{
   data(cosmesis)
   csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
   ml1 <- Maclist(csub1)
   mm1 <- Macmat(ml1)
}
\keyword{manip }
