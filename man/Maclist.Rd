\name{Maclist}
\alias{Maclist}
\title{ A function to  }
\description{
  Returns a list of maximal cliques of the intersction graph of the
  real valued intervals supplied in \code{m}. These are one dimensional
  intervals with one interval for each individual. The algorithm is
  coded in interpreted code and should be moved to compiled code for speed.
  How do we handle exact failure times?
  Which algorithm is used?
}
\usage{
Maclist(intvls, Lopen=TRUE, Ropen=FALSE)
}
\arguments{
  \item{intvls}{ A n by 2 matrix, the first column is the left endpoints and
  the second column contains the right endpoints of the failure time
  intervals.}
  \item{Lopen}{A boolean indicating whether the intervals are open on
    the left.}
  \item{Ropen}{ A boolean indicating whether the intervals are open on
    the right.}
}
\value{
  A list of length m. Each element of the list corresponds to one
  maximal antichain. The row numbers (from \code{m}) identify the
  individuals and all row numbers for the individuals in the maximal
  clique. Maximal cliques occur in their natural (left to right) order.
}
\references{Computational Methods for Censored Data using Intersection
  Graphs, R. Gentleman and A. Vandal, JCGS, 2000. }
\author{ Alain Vandal and Robert Gentleman }

\seealso{ \code{\link{Macmat}} }

\examples{
   data(cosmesis)
   csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
   ml1 <- Maclist(csub1)
}

\keyword{manip }
