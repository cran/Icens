\name{Plotboxes}
\alias{Plotboxes}
\title{ Plot the event time regions for bivariate data. }
\description{ Plot rectangles described by the interval given in the
  first two arguments.}
\usage{
Plotboxes(int1, int2, textp=FALSE, showmac=FALSE, showsupp=FALSE, showmp=FALSE,
cliques=NULL, macprod=NULL, density=c(2, 8, 20), col=c(2, 3, 4),
offsetx=0.02, offsety=0.03) 
}

\arguments{
  \item{int1}{ The intervals for the x dimension. }
  \item{int2}{ The intervals for the y dimension. }
  \item{textp}{ Boolen, if true add text. }
  \item{showmac}{ Boolean, if true then the maximal cliques are shown in
  a different color? }
  \item{showsupp}{ Boolean, if true show support boxes. }
  \item{showmp}{ Boolean }
  \item{cliques}{ Maximal cliques. }
  \item{macprod}{ macprod }
  \item{density}{ The density of the polygon shading lines, in lines per inch. }
  \item{col}{ Color for plotting features. }
  \item{offsetx}{ Offset for x-axis. }
  \item{offsety}{ Offset for y-axis. }
}

\value{
  No value is returned. The event rectangles are plotted on the active
  graphics device.
}
\references{ \emph{Graph--Theoretical Aspects of Bivariate Censored
    Data}, R. Gentleman and A. Vandal, 1999, submitted. }
\author{ A. Vandal and R. Gentleman }
\seealso{ \code{\link{BVclmat}}, \code{\link{BVsupport}},
  \code{\link{BVcliques}}  }

\examples{
   data(cmv)
   Plotboxes(cmv[,1:2], cmv[,3:4], showmac=TRUE)
}
\keyword{aplot}



