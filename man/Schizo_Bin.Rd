\name{Schizo_Bin}
\alias{Schizo_Bin}
\docType{data}
\title{
Data of five clinical trials in Schizophrenia (with binary outcomes).
}
\description{
These are the data of five clinical trials in Schizophrenia. A total of \eqn{2128} patients were treated by \eqn{198} investiagators (psychiatrists). Patients' Schizophrenia symptoms at baseline and at the end of the study were measured using the PANSS and BPRS. 
The variables BPRS_Bin and PANSS_Bin are binary outcomes that indicate whether clinically meaningful change had occurred (1 = a reduction of \eqn{20}\% or higher in the PANSS/BPRS scores at the last measurement compared to baseline; 0 = no such reduction; Leucht et al., 2005; Kay et al., 1988).
}
\usage{data(Schizo_Bin)}
\format{
A \code{data.frame} with \eqn{2128} observations on \eqn{5} variables.
  \describe{
    \item{\code{Id}}{The patient ID.}
    \item{\code{InvestID}}{The ID of the investigator (psychiatrist) who treated the patient.}
    \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = control and \eqn{1} = Risperidone.}
    \item{\code{PANSS_Bin}}{The dichotomized change in the PANSS score (1 = a
reduction of \eqn{20}\% or more in the PANSS score, 0=otherwise)}
    \item{\code{BPRS_Bin}}{The dichotomized change in the BPRS score (1 = a
reduction of \eqn{20}\% or more in the BPRS score, 0=otherwise)}
}
}
\references{
Kay, S.R., Opler, L.A., & Lindenmayer, J.P. (1988). Reliability and validity of the Positive and Negative Syndrome Scale for schizophrenics. Psychiatric Research, 23, 99-110.

Leucht, S., et al. (2005). Clinical implications of Brief Psychiatric Rating Scale scores. The British Journal of Psychiarty, 187, 366-371.
}
\keyword{Schizo_Bin}
