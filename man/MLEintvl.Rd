\name{MLEintvl}
\alias{MLEintvl}
\title{ Compute the real representation for the maximal cliques. }
\description{
  The intervals on the real line that corresponds to the intersections
  of the maximal cliques are computed and returned.
}
\usage{
MLEintvl(intvls, ml=Maclist(intvls))
}
\arguments{
  \item{intvls}{ The n by 2 matrix containing the event time intervals
    for the individuals under study. }
  \item{ml}{ The \code{\link{Maclist}} computed for the \code{intvls}. }
}
\value{
 An m by 2 matrix, where m is the number of maximal cliques. The first
 column contains the left end point of the real representation for the
 appropriate maximal clique and the second column contains the right end
 point. 
}
\references{ Computational Methods for Censored Data using Intersection
  Graphs, R. Gentleman and A. Vandal, JCGS, 2000. }
\author{Alain Vandal and Robert Gentleman  }
\seealso{ \code{\link{Maclist}} }

\examples{
   data(cosmesis)
   csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
   MLEintvl(csub1)
}
\keyword{manip }
