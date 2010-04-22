\name{plot.icsurv}
\alias{plot.icsurv}
\title{ A plot method for the estimates produced by the estimation
  methods in Icens. }
\description{
  Produces nice plots of the estimated NPMLE.
}
\usage{
\S3method{plot}{icsurv}(x, type="eq", surv=FALSE, bounds=FALSE, shade=3, density=30,
angle=45, lty=1, new=TRUE, xlab="Time", ylab="Probability", main="GMLE",
ltybnds=2, ...)
}
\arguments{
  \item{x}{ The estimate of the NPMLE. }
  \item{type}{ Three options, "eq" for equivalence call, "gw" for the
    Groeneboom-Wellner estimate, and "lc" for the left-continuous estimate. }
  \item{surv}{ Logical indicating whether or not to plot the survival curve. }
  \item{bounds}{ Logical indicating whether or not to include bounds
    around the estimate. }
  \item{shade}{ An integer in 1, 2, or 3 denoting what component of the
    plot to shade. }
  \item{density}{ The density of shading lines, in lines per inch. }
  \item{angle}{ The slope of shading lines, given as an angle in degrees
    (counter-clockwise). }
  \item{lty}{ The line type for the estimates. }
  \item{new}{ Logical indicating whether or not to create a new plot. }
  \item{xlab}{ The x-axis label. }
  \item{ylab}{ The y-axis label. }
  \item{main}{ The main title for the plot. }
  \item{ltybnds}{ The line type for the bounds on the estimates. }
  \item{\dots}{ Additional arguments passed to the plot function. }
}

\value{
  No value is returned. A plot of the NPMLE is made on the active
  graphics device.
}

\author{Alain Vandal and Robert Gentleman.  }

\seealso{\code{\link{VEM}}, \code{\link{ISDM}}, \code{\link{EMICM}},
  \code{\link{PGM}}  }

\examples{
    data(cosmesis)
    csub1 <- subset(cosmesis, subset=Trt==0, select=c(L,R))
    e1 <- VEM(csub1)
    par(mfrow=c(2,2))
    plot(e1)
    data(pruitt)
    e2 <- EM(csub1)
    plot(e2)
    e3 <- PGM(csub1)
    plot(e3)
    e4 <- EMICM(csub1)
    plot(e4)
}
\keyword{hplot}
