\name{VEM}
\alias{VEM}
\title{ Compute the NPMLE of p via the Vertex Exchange Method. }
\description{
  The Vertex Exchange Method is used to obtain the NPMLE of \code{p}.
}
\usage{
VEM(A, pvec, maxiter=500, tol=1e-07, tolbis=1e-07, keepiter=FALSE)
}
\arguments{
  \item{A}{ The m by n incidence matrix or the n by 2 matrix of intervals. }
  \item{pvec}{ The initial estimate for the probability vector. }
  \item{maxiter}{ The maximum number of iterations allowed. }
  \item{tol}{ The tolerance used to determine convergence. }
  \item{tolbis}{ The tolerance used in the bisection stage of the
    algorithm. }
  \item{keepiter}{Should iteration information be retained and returned.}
} 
\details{
  Lots.
}
\value{
  An object of class \code{\link{icsurv}} with the following components.
  \item{pf }{The NPMLE of the probability vector. }
  \item{numiter }{The number of iterations used.}
  \item{lval}{The value of the logarithm of the likelihood at the
    NPMLE. }
  \item{converge}{Boolean stating whether the iteration converged.}
  \item{intmap}{If present it contains the real representations for the
    maximimal cliques. These are the intervals (on the real line) where
    the mass in \code{pf} is placed.}
}
\references{ \emph{A Vertex-exchange-method in $D$-optimal Design Theory
    }, D. Bohning, Metrika, 1986. }
\author{ Robert Gentleman and Alain Vandal }

\seealso{ \code{\link{EM}}, \code{\link{ISDM}}, \code{\link{EMICM}},
     \code{\link{PGM}} }

\examples{
    data(cosmesis)
    csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
    VEM(csub1)
    data(pruitt)
    VEM(pruitt)
}
\keyword{optimize}



