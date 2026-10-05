#' Longitudinal PANSS and BPRS data from a schizophrenia clinical trial
#'
#' These are longitudinal data from a single clinical trial in schizophrenia,
#' prepared for use in the single-trial setting. A total of 450 patients were
#' included in the data set. Patients' schizophrenia symptoms were measured using
#' the PANSS and BPRS at baseline and at weeks 1, 2, 4, 6, and 8 after the start
#' of treatment. There were two treatment conditions: haloperidol (control) and
#' risperidone (treatment).
#'
#' @format A data frame with 3962 rows and 9 variables:
#' \describe{
#'   \item{Id}{The anonymized patient ID.}
#'   \item{Treat}{The treatment indicator, coded as 0 = control/haloperidol and 1 = treatment/risperidone.}
#'   \item{Time}{The assessment week, with values 1, 2, 4, 6, and 8.}
#'   \item{Endpoint}{The endpoint indicator, either BPRS or PANSS.}
#'   \item{Response}{The observed endpoint value at the corresponding assessment week.}
#'   \item{Response_log}{The log-transformed observed endpoint value.}
#'   \item{Time_cont}{The assessment week represented as a numeric time variable.}
#'   \item{Baseline}{The baseline value of the corresponding endpoint, measured at time 0.}
#'   \item{Baseline_log}{The log-transformed baseline value of the corresponding endpoint, measured at time 0.}
#' }
#'
#' @source Longitudinal clinical trial data in schizophrenia.
#' @docType data
#' @name Schizo_PANSS_BPRS
#' @keywords datasets
#'
"Schizo_PANSS_BPRS"
