\name{PGM}
\alias{PGM}
\title{ An implementation of the projected gradient methods for finding
  the NPMLE. }
\description{
  An estimate of the NPMLE is obtained by using projected gradient
  methods. This method is a special case of the methods described in Wu
  (1978). 
}
\usage{
PGM(A, pvec, maxiter = 500, tol=1e-07, told=2e-05, tolbis=1e-08,
    keepiter=FALSE)
}
\arguments{
  \item{A}{ \code{A} is either the m by n clique matrix or the n by 2
    matrix containing the left and right end points for each event time. }
  \item{pvec }{ An initial estimate of the probability vector. }
  \item{maxiter }{ The maximum number of iterations to take. }
  \item{tol}{ The tolerance for decreases in likelihood. }
  \item{told}{ \code{told} does not seem to be used. }
  \item{tolbis}{ The tolerance used in the bisection code. }
  \item{keepiter}{ A boolean indicating whether to return the number of
    iterations.  }
}
\details{
  New directions are selected by the projected gradient method. The new
  optimal \code{pvec} is obtained using the bisection algorithm, moving
  in the selected direction. Convergence requires both the \eqn{L_1}
  distance for the improved \code{pvec} and the change in likelihood to
  be below \code{tol}.
}
\value{
  An object of class \code{\link{icsurv}} containing the following
  components: 
  \item{pf }{The NPMLE of \code{pvec}. }
  \item{sigma }{The cummulative sum of \code{pvec}.}
  \item{lval }{ The value of the log likelihood at \code{pvec}.} 
  \item{clmat }{ The clique matrix. }
  \item{method }{ The method used, currently only "MPGM" is possible. }
  \item{lastchange }{ The difference between \code{pf} and the previous
    iterate. }
  \item{numiter}{ The number of iterations carried out. }
  \item{eps}{ The tolerances used. }
  \item{converge }{ A boolean indicating whether convergence occurred
    within \code{maxiter} iterations. }
  \item{iter }{ If \code{keepiter} is true then this is a matrix
    containing all iterations - useful for debugging.}
}
\references{ \emph{Some Algorithmic Aspects of the Theory of Optimal
    Designs}, C.--F. Wu, 1978, Annals. }
\author{ Alain Vandal and Robert Gentleman. }

\seealso{ \code{\link{VEM}}, \code{\link{ISDM}}, \code{\link{EMICM}},
     \code{\link{PGM}}, \code{\link{EM}} }

\examples{
    data(cosmesis)
    csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
    PGM(csub1)
    data(pruitt)
    PGM(pruitt)
}
\keyword{optimize}
