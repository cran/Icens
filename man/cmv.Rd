\name{cmv}
\alias{cmv}
\non_function{}
\title{ Data on times to shedding of cytomegalovirus and to colonization
  of mycobacterium avium complex.  }
\description{
  The \code{cmv} data frame has 204 rows and 4 columns.
  The intervals should be treated as closed at both ends to replicate
  the analysis in Betensky and Finkelstein.
}
\format{
This data frame contains the following columns:
\describe{
  \item{cmvL}{
    The left end of the CMV shedding interval. 
  }
  \item{cmvR}{
    The right end of the CMV shedding interval.
  }
  \item{macL}{
    The left end of the MAC colonization interval.
  }
  \item{macR}{
    The right end of the MAC colonization interval.
  }
}
}
\details{
  Betensky and Finkelstein, 1999 present data from the AIDS Clinical
Trials Group protocol ACTG 181. This was a natural history substudy of
a comparative trial. Patients were scheduled for clinic visits during
follow--up and data was collected on the time until two events;
shedding of cytomegalovirus (CMV) in the urine and blood and for
colonization of mycobacterium avium complex (MAC) in the sputum or
stool. 
}
\source{
Betensky, R. A. and Finkelstein, D. M., 1999,  
\emph{A nonparametric maximum likelihood estimator for bivariate interval
censored data},
Statistics in Medicine,
}
\examples{
data(cmv)
}
\keyword{datasets}
