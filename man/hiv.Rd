\name{hiv}
\alias{hiv}
\non_function{}
\title{ Intervals for infection time and disease onset for 257
  hemophiliac patients.  }
\description{
The \code{hiv} data frame has 257 rows and 4 columns.
}
\format{
This data frame contains the following columns:
\describe{
\item{yL}{
The left end point of the infection time interval.
}
\item{yR}{
The right end point of the infection time interval.
}
\item{zL}{
The left end point of the disease onset interval.
}
\item{zR}{
The right end point of the disease onset interval.
}
\item{Age}{ Coded as 1 if the estimated age at infection was less than 20 
  and 2 if the estimated age at infection was greater than 20.}
\item{Trt}{ Treatment, Light or Heavy}
}
}
\details{
  The setting is as follows. Individuals were infected with the HIV virus
at some unknown time they subsequently develop AIDS at a second
unknown time. The data consist of two intervals, \eqn{(y_L, y_R)}
and \eqn{(z_L,z_R)}, such that the infection time was in the first
interval and the time of disease onset was in the second interval.
A quantity of interest is the incubation time of the disease which is
\eqn{T=Z-Y}. The authors argue persuasively that
this should be considered as bivariate interval censored data.
They note that simply forming the differences \eqn{(z_L-y_R, z_R-y_L)} and
analysing the
resultant data assumes an incorrect likelihood.
DeGruttola and Lagakos transform the problem slightly to study the joint
distribution of \eqn{Y} and \eqn{T=Z-Y}.
This is equivalent to estimating the joint distribution of \eqn{Z} and
\eqn{Y} then transforming.
The data, as reported, have been discretized into six month intervals.

We use the data as reported in Table 1 of DeGruttola and Lagakos, 1989.
The patients were 257 persons with Type A or B hemophilia treated at
two hospitals in France.
They were then examined intermittently (as they came in for treatment?)
and their HIV and AIDS status was determined.
Kim, De Gruttola and Lagakos report some covariate information and
their paper is concerned with the modeling of that information. In
this paper we concentrate only on the event times and ignore the
covariate information; that topic being worthy of separate
investigation.

}
\source{
DeGruttola, V. and Lagakos, S.W.,
1989,
\emph{Analysis of doubly-censored survival data, with application to {AIDS}},
Biometrics.

Kim, Mimi Y. and De Gruttola, Victor G. and Lagakos, Stephen W., 1993,
\emph{Analyzing Doubly Censored Data With Covariates, With Application
  to  AIDS},
Biometrics.

}
\examples{
data(hiv)
}
\keyword{datasets}
