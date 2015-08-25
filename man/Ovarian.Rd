\name{Ovarian}
\alias{Ovarian}
\docType{data}
\title{The Ovarian dataset}
\description{
This dataset combines the data that were collected in four double-blind randomized clinical trials in advanced ovarian cancer (Ovarian Cancer Meta-Analysis Project, 1991). In these trials, the objective was to examine the efficacy of cyclophosphamide plus cisplatin (CP) versus  cyclophosphamide plus adriamycin plus cisplatin (CAP) to treat advanced ovarian cancer.}
\usage{data("Ovarian")}
\format{
  A data frame with 1192 observations on the following 7 variables.
  \describe{
    \item{\code{Patient}}{The ID number of a patient.}
    \item{\code{Center}}{The center in which a patient was treated.}
    \item{\code{Treat}}{The treatment indicator, coded as 0=CP (active control) and 1=CAP (experimental treatment).}
    \item{\code{Pfs}}{Progression-free survival (the candidate surrogate).}
    \item{\code{PfsInd}}{Censoring indicator for progression-free survival.}
    \item{\code{Surv}}{Survival time (the true endpoint).}
    \item{\code{SurvInd}}{Censoring indicator for survival time.}
  }
}
\references{
Ovarian Cancer Meta-Analysis Project (1991). Cclophosphamide plus cisplatin plus adriamycin versus cyclophosphamide, doxorubicin, and cisplatin chemotherapy of ovarian carcinoma: a meta-analysis. \emph{Classic papers and current comments, 3,} 237-234.
}
\examples{
data(Ovarian)
str(Ovarian)
head(Ovarian)
}
\keyword{datasets}
