\name{EM}
\alias{EM}
\title{ A function to compute the NPMLE of p based on the incidence
  matrix A. }
\description{
 The incidence matrix, \code{A} is the m by n matrix that represents the
 data. There are m probabilities that must be estimated. The EM, or
 expectation maximization, method is applied to these data.
}
\usage{
EM(A, pvec, maxiter=500, tol=1e-12)
}

\arguments{
  \item{A}{ The incidence matrix. }
  \item{pvec}{ The probability vector. }
  \item{maxiter}{ The maximum number of iterations. }
  \item{tol}{ The tolerance used to judge convergence. }
}
\details{
  Lots.
}
\value{
 An object of class \code{\link{icsurv}} containing the following
 components: 
  \item{pf }{The NPMLE of the probability vector.}
  \item{numiter }{The number of iterations used.}
  \item{converge }{ A boolean indicating whether the algorithm
    converged.}
  \item{intmap }{ If present indicates the real representation of the
    support for the values in \code{pf}. }
}
\references{The EM algorithm applied to the maximal cliques of the
  intersection graph of the censored data. \emph{The empirical distribution
    function with arbitrarily grouped, censored and truncated data},
  B. W. Turnbull, 1976, JRSS;B.}
\author{ Alain Vandal and Robert Gentleman. }

\seealso{ \code{\link{VEM}}, \code{\link{ISDM}}, \code{\link{EMICM}},
     \code{\link{PGM}} }

\examples{
    data(cosmesis)
    csub1 <- subset(cosmesis, subset= Trt==0, select=c(L,R))
    EM(csub1)
    data(pruitt)
    EM(pruitt)
}
\keyword{nonparametric }




