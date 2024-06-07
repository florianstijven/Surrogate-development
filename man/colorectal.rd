\name{colorectal}
\alias{colorectal}
\docType{data}
\title{The Colorectal dataset with a binary surrogate.}
\description{
This dataset combines the data that were collected in 26 double-blind randomized clinical trials in advanced colorectal cancer.
}
\usage{data("colorectal")}
\format{
  A data frame with 3943 observations on the following 7 variables.
  \describe{
    \item{\code{TRIAL}}{The ID number of a trial.}
    \item{\code{responder}}{Binary tumor response (the candidate surrogate), coded as 2=complete response (CR) or partial response (PR) and 1=stabled disease (SD) or progressive disease (PD).}
     \item{\code{SURVIND}}{Censoring indicator for survival time.}
    \item{\code{TREAT}}{The treatment indicator, coded as 0=active control and 1=experimental treatment.}
    \item{\code{CENTER}}{The center in which a patient was treated. In this dataset, there was only one center per trial, hence TRIAL=CENTER.}
     \item{\code{patientid}}{The ID number of a patient.}
     \item{\code{surv}}{Survival time (the true endpoint).}
  }
}
\references{
Alonso, A., Bigirumurame, T., Burzykowski, T., Buyse, M., Molenberghs, G., Muchene, L., ... & Van der Elst, W. (2016). Applied surrogate endpoint evaluation methods with SAS and R. CRC Press.
}
\examples{
data(colorectal)
str(colorectal)
head(colorectal)
}
\keyword{datasets}
