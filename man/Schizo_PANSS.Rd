\name{Schizo_PANSS}
\alias{Schizo_PANSS}
\docType{data}
\title{
Longitudinal PANSS data of five clinical trials in schizophrenia
}
\description{
These are the longitudinal PANSS data of five clinical trial in schizophrenia. A total of \eqn{2151} patients were treated by \eqn{198} investiagators (psychiatrists). There were two treatment conditions (risperidone and control). Patients' schizophrenic symptoms were measured using the PANSS at different time moments following start of the treatment. The variables Week1-Week8 express the change scores over time using the raw (semi-continuous) PANSS scores. The variables Week1_bin - Week8_bin are binary indicators of a \eqn{20}\% or higher reduction in PANSS score versus baseline. The latter corresponds to a commonly accepted criterion for defining a clinically meaningful response (Kay et al., 1988).}
\usage{data(Schizo_PANSS)}
\format{
A \code{data.frame} with \eqn{2151} observations on \eqn{6} variables.
  \describe{
    \item{\code{Id}}{The patient ID.}
    \item{\code{InvestID}}{The ID of the investigator (psychiatrist) who treated the patient.}
    \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = placebo and \eqn{1} = Risperidone.}
    \item{\code{Week1}}{The change in the PANSS score \eqn{1} week after starting the treatment (= score at the end of the treatment - score at \eqn{1} week after starting the treatment).}
    \item{\code{Week2}}{The change in the PANSS score \eqn{2} weeks after starting the treatment.}
    \item{\code{Week4}}{The change in the PANSS score \eqn{4} weeks after starting the treatment.}
    \item{\code{Week6}}{The change in the PANSS score \eqn{6} weeks after starting the treatment.}
    \item{\code{Week8}}{The change in the PANSS score \eqn{8} weeks after starting the treatment.}
    \item{\code{Week1_bin}}{The dichotomized change in the PANSS score \eqn{1} week after starting the treatment (\eqn{1}=a \eqn{20}\% or higher reduction in PANSS score versus baseline, \eqn{0}=otherwise).}
   \item{\code{Week2_bin}}{The dichotomized change in the PANSS score \eqn{2} weeks after starting the treatment.} 
   \item{\code{Week4_bin}}{The dichotomized change in the PANSS score \eqn{4} weeks after starting the treatment.} 
   \item{\code{Week6_bin}}{The dichotomized change in the PANSS score \eqn{6} weeks after starting the treatment.} 
   \item{\code{Week8_bin}}{The dichotomized change in the PANSS score \eqn{8} weeks after starting the treatment.} 
}
}
\references{
Kay, S.R., Opler, L.A., & Lindenmayer, J.P. (1988). Reliability and validity of the Positive and Negative Syndrome Scale for schizophrenics. Psychiatric Research, 23, 99-110.
}
\keyword{Schizo_PANSS}
