\name{cosmesis}
\alias{cosmesis}
\non_function{}
\title{ The time taken until cosmetic deterioration of breast cosmesis. }
\description{
  The \code{cosmesis} data frame has 95 rows and 3 columns.
 

}
\format{
This data frame contains the following columns:
\describe{
\item{L}{
The left end point of the cosmetic deterioration interval.
}
\item{R}{
The right end point of the cosmetic deterioration interval.
}
\item{Trt}{
The treatment indicator. It is zero for those that received radiotherapy.
}
}
}

\source{
 \emph{A semiparametric model for regression analysis of
   interval-censored failure time data,} D. M. Finkelstein and
 R. A. Wolfe, 1985, Biometrics.
}
\examples{
data(cosmesis)
}
\keyword{datasets}
