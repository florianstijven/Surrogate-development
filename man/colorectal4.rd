\name{colorectal4}
\alias{colorectal4}
\docType{data}
\title{The Colorectal dataset with an ordinal surrogate.}
\description{
This dataset combines the data that were collected in 19 double-blind randomized clinical trials in advanced colorectal cancer.
}
\usage{data("colorectal4")}
\format{
  A data frame with 3192 observations on the following 7 variables.
  \describe{
    \item{\code{trialend}}{The ID number of a trial.}
    \item{\code{treatn}}{The treatment indicator, coded as 0=active control and 1=experimental treatment.}
    \item{\code{trueind}}{Censoring indicator for survival time.}
    \item{\code{surrogend}}{Categorical ordered tumor response (the candidate surrogate), coded as 1=complete response (CR), 2=partial response (PR), 3=stabled disease (SD) and 4=progressive disease (PD).}
    \item{\code{patid}}{The ID number of a patient.}
    \item{\code{center}}{The center in which a patient was treated. In this dataset, there was only one center per trial, hence TRIAL=CENTER.}
    \item{\code{truend}}{Survival time (the true endpoint).}
  }
}
\references{
Alonso, A., Bigirumurame, T., Burzykowski, T., Buyse, M., Molenberghs, G., Muchene, L., ... & Van der Elst, W. (2016). Applied surrogate endpoint evaluation methods with SAS and R. CRC Press.
}
\examples{
data(colorectal4)
str(colorectal4)
head(colorectal4)
}
\keyword{datasets}
