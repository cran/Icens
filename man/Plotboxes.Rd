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
  \item{textp}{ }
  \item{showmac}{ Boolean, if true then the maximal cliques are shown in
  a different colour? }
  \item{showsupp}{  }
  \item{showmp}{ }
  \item{cliques}{  }
  \item{macprod}{  }
  \item{density}{  }
  \item{col}{  }
  \item{offsetx}{  }
  \item{offsety}{  }
}

\value{
  No value is returned. The event rectangles are plotted on the active
  graphics device.
}
\references{ \emph{Graph--Theoretical Aspects of Bivariate Censored
    Data}, R. Gentleman and A. Vandal, 1999, submitted. }
\author{ A. Vandal and R. Gentleman }
\seealso{ \code{\link{BVclmat}}, \code{\ling{BVsupport}},
  \code{\link{BVcliques}}  }

\examples{
   data(cmv)
   Plotboxes(cmv[,1:2], cmv[,3:4], showmac=TRUE)
}
\keyword{aplot}



