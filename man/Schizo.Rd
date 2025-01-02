\name{Schizo}
\alias{Schizo}
\docType{data}
\title{
Data of five clinical trials in schizophrenia
}
\description{
These are the data of five clinical trials in schizophrenia. A total of \eqn{2128} patients were treated by \eqn{198} investiagators (psychiatrists). Patients' schizophrenic symptoms were measured using the PANSS, BPRS, and CGI. There were two treatment conditions (risperidone and control).}
\usage{data(Schizo)}
\format{
A \code{data.frame} with \eqn{2128} observations on \eqn{9} variables.
  \describe{
    \item{\code{Id}}{The patient ID.}
    \item{\code{InvestID}}{The ID of the investigator (psychiatrist) who treated the patient.}
    \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = control and \eqn{1} = Risperidone.}
    \item{\code{CGI}}{The change in the CGI score (= score at the start of the treatment - score at the end of the treatment).}
    \item{\code{PANSS}}{The change in the PANSS score.}
    \item{\code{BPRS}}{The change in the BPRS score.}
    \item{\code{PANSS_Bin}}{The dichotomized PANSS change score, coded as \eqn{1} = a reduction of
20\% or more in the PANSS score (score at the end of the treatment relative to score at the beginning of the treatment), \eqn{0} = otherwise.}
    \item{\code{BPRS_Bin}}{The dichotomized BPRS change score, coded as \eqn{1} = a reduction of 20\% or more in the BPRS score (score at the end of the treatment relative to score at the beginning of the treatment), \eqn{0} = otherwise.}
    \item{\code{CGI_Bin}}{The sichtomized change in the CGI score, coded as \eqn{1} = a change of more than \eqn{3} points on the original scale (score at the end of the treatment relative to score at the beginning of the treatment), \eqn{0} = otherwise.}
}
}
\keyword{Schizo}
