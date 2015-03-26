\name{Schizo_Bin_ST}
\alias{Schizo_Bin_ST}
\docType{data}
\title{
Data of a clinical trial in Schizophrenia (with binary outcomes, dataset \eqn{1}).
}
\description{
These are the data of a clinical trial in Schizophrenia (a subset of the dataset \code{Schizo_Bin}). A total of \eqn{902} patients were treated by \eqn{131} investigators (psychiatrists). Patients' schizophrenia symptoms at baseline and at the end of the study (after \eqn{8} weeks) were measured using the PANSS and BPRS. 
The variables BPRS_Bin and PANSS_Bin are binary outcomes that indicate whether clinically meaningful change had occurred (1 = a reduction of \eqn{20}\% or higher in the PANSS/BPRS scores at the last measurement compared to baseline; 0 = no such reduction; Leucht et al., 2005; Kay et al., 1988).
}
\usage{data(Schizo_Bin_ST)}
\format{
A \code{data.frame} with \eqn{448} observations on \eqn{5} variables.
  \describe{
    \item{\code{Id}}{The patient ID.}
    \item{\code{InvestID}}{The ID of the investigator (psychiatrist) who treated the patient.}
    \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = control treatment (10 mg. haloperidol) and \eqn{1} = experimental treatment (8, 12, or 16 mg. risperidone).}
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
\keyword{Schizo_Bin_ST}
