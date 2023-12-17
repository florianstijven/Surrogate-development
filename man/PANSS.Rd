\name{PANSS}
\alias{PANSS}
\docType{data}
\title{
PANSS subscales and total score based on the data of five clinical trials in schizophrenia
}
\description{
These are the PANSS subscale and total scale scores of five clinical trial in schizophrenia. A total of \eqn{1941} patients were treated by \eqn{126} investiagators (psychiatrists). There were two treatment conditions (risperidone and control). Patients' schizophrenic symptoms were measured using the PANSS (Kay et al., 1988).}
\usage{data(PANSS)}
\format{
A \code{data.frame} with \eqn{1941} observations on \eqn{9} variables.
  \describe{
    \item{\code{Pat.Id}}{The patient ID.}
        \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = active control and \eqn{1} = Risperidone.}
    \item{\code{Invest}}{The ID of the investigator (psychiatrist) who treated the patient.}
    \item{\code{Neg}}{The Negative symptoms scale score.}
    \item{\code{Exc}}{The Excitement scale score.}
    \item{\code{Cog}}{The Cognition scale score.}
    \item{\code{Pos}}{The Positive symptoms scale score.}
    \item{\code{Dep}}{The Depression scale score.}
    \item{\code{Total}}{The Total PANSS score.}
  }
}
\references{
Kay, S.R., Opler, L.A., & Lindenmayer, J.P. (1988). Reliability and validity of the Positive and Negative Syndrome Scale for schizophrenics. Psychiatric Research, 23, 99-110.
}
\keyword{PANSS}
