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
A \code{data.frame} with \eqn{2128} observations on \eqn{6} variables.
  \describe{
    \item{\code{Id}}{The patient ID.}
    \item{\code{InvestID}}{The ID of the investigator (psychiatrist) who treated the patient.}
    \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = control and \eqn{1} = Risperidone.}
    \item{\code{CGI}}{The change in the CGI score (= score at the start of the treatment - score at the end of the treatment).}
    \item{\code{PANSS}}{The change in the PANSS score (= score at the start of the treatment - score at the end of the treatment).}
    \item{\code{BPRS}}{The change in the PANSS score (= score at the start of the treatment - score at the end of the treatment).}
}
}
\keyword{Schizo}
