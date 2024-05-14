\name{prostate}
\alias{prostate}
\docType{data}
\title{The prostate dataset with a continuous surrogate.
}
\description{
This dataset combines the data that were collected in 17 double-blind randomized clinical trials in advanced prostate cancer.
}
\usage{data("prostate")}
\format{
  A data frame with 412 observations on the following 6 variables.
  \describe{
    \item{\code{TRIAL}}{The ID number of a trial.}
    \item{\code{TREAT}}{The treatment indicator, coded as 0=active control and 1=experimental treatment.}
    \item{\code{PSA}}{Prostate specific antigen (surrogate endpoint)}
    \item{\code{SURVTIME}}{Survival time (the true endpoint).}
    \item{\code{SURVIND}}{Censoring indicator for survival time.}
    \item{\code{PATID}}{The ID number of a patient.}
  }
}

\references{
Alonso A, Bigirumurame T, Burzykowski T, Buyse M, Molenberghs G, Muchene L, Perualila NJ, Shkedy Z, Van der Elst W, et al. (2016). Applied surrogate endpoint evaluation methods with SAS and R. CRC Press New York
}
\examples{
data(prostate)
str(prostate)
head(prostates)
}
\keyword{datasets}
