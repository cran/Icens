\name{ISDM}
\alias{ISDM}
\title{ Estimate the NPMLE of censored data using the ISDM method
  proposed in Lesperance and Kalbfleisch (19 }
\description{
  ISDM is a method for estimating the NPMLE of censored data.
}
\usage{
ISDM(A, pvec, maxiter=500, tol=1e-07, tolbis=1e-08, verbose=FALSE)
}
\arguments{
  \item{A}{ The m by n incidence, or clique, matrix. Or the n by 2
    matrix containing the event intervals.}
  \item{pvec}{ An initial estimate of the probability vector; not required. }
  \item{maxiter}{Maximum number of iterations to be made. }
  \item{tol}{ The tolerance used to determine convergence. }
  \item{tolbis}{ A second tolerance used for the steps. }
  \item{verbose}{ Boolean, should verbose output be printed. }
}
\details{
  Lots of complicated stuff should go here.
}
\value{
A list containing:
  \item{pf }{The estimated NPMLE of the probability vector.}
  \item{numiter }{The number of iterations performed.}
}
\references{ An Algorithm for Computing the Nonparametric MLE of a Mixing
  Distribution, Lesperance, Mary L. and Kalbfleisch, John D., JASA,
  1992 }
\author{Alain Vandal and Robert Gentleman }

\seealso{ \code{\link{VEM}}, \code{\link{EMICM}},
     \code{\link{PGM}} }

\examples{
    data(cosmesis)
    csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
    ISDM(csub1)
#    data(pruitt)
#    ISDM(pruitt)  
}
\keyword{optimize}
