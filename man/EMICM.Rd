\name{EMICM}
\alias{EMICM}
\title{ Compute the NPMLE for censored data using the EMICM. }
\description{
  An implementation of the hybrid EM ICM (Iterative convex minorant)
  estimator of the distribution function proposed by Wellner and Zahn (1997).
}
\usage{
EMICM(A, EMstep=TRUE, ICMstep=TRUE, keepiter=FALSE, tol=1e-07,
maxiter=1000) 
}
\arguments{
  \item{A}{ Either the m by n clique matrix or the n by 2 matrix
    containing the event time intervals. }
  \item{EMstep}{ Boolean, indicating whether to take an EM step in the
    iteration. }
  \item{ICMstep}{ Boolean, indicating whether to take an ICM step. }
  \item{keepiter}{  Boolean determining whether to keep the iteration
    states. }
  \item{tol}{ The maximal L1 distance between successive estimates
    before stopping iteration. }
  \item{maxiter}{ The maximal number of iterations to perform before
    stopping. }
}
\details{
  Lots, and they're complicated too!
}
\value{
  An  object of class \code{\link{icsurv}} containing the following
  components:
  \item{pf }{ The estimated probabilities.}
  \item{sigma }{ The NPMLE of the survival function on the maximal
    antichains. }
  \item{weights }{ The diagonal of the likelihood function's second
    derivative. }
  \item{lastchange }{ A vector of differences between the last two
    iterations. }
  \item{numiter }{ The total number of iterations performed.}
  \item{iter }{ Is only present if \code{keepiter} is true; states of
    sigma  during the iteration.}
  \item{intmap }{ The real representation associated with the
    probabilities reported in \code{pf}.}
}
\references{\emph{A hybrid algorithm for computation of the nonparametric
maximum likelihood estimator from censored data}, J. A. Wellner and
Y. Zhan, 1997, JASA.  }
\author{ Alain Vandal and Robert Gentleman }
\seealso{ \code{\link{EM}},\code{\link{VEM}}, \code{\link{PGM}} }

\examples{
  data(cosmesis)
  csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
  EMICM(csub1)
  data(pruitt)
  EMICM(pruitt)
}
\keyword{optimize}
