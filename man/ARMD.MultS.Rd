\name{ARMD.MultS}
\alias{ARMD.MultS}
\docType{data}
\title{
Data of the Age-Related Macular Degeneration Study with multiple candidate surrogates
}
\description{
These are the data of a clinical trial involving patients suffering from age-related macular degeneration (ARMD), a condition that involves a progressive loss of vision. A total of \eqn{181} patients participated in the trial. Patients' visual acuity was assessed using standardized vision charts. There were two treatment conditions (placebo and interferon-\eqn{\alpha}). The potential surrogate endpoints are the changes in the visual acuity at \eqn{4}, \eqn{12}, and \eqn{24} weeks after starting treatment. The true endpoint is the change in the visual acuity at \eqn{52} weeks.}
\usage{data(ARMD.MultS)}
\format{
A \code{data.frame} with \eqn{181} observations on \eqn{6} variables.
  \describe{
    \item{\code{Id}}{The Patient ID.}
    \item{\code{Diff4}}{The change in the visual acuity at \eqn{4} weeks after starting treatment. This endpoint is a potential surrogate for \code{Diff52}.}
    \item{\code{Diff12}}{The change in the visual acuity at \eqn{12} weeks after starting treatment. This endpoint is a potential surrogate for \code{Diff52}.}
    \item{\code{Diff24}}{The change in the visual acuity at \eqn{24} weeks after starting treatment. This endpoint is a potential surrogate for \code{Diff52}.}
    \item{\code{Diff52}}{The change in the visual acuity at \eqn{52} weeks after starting treatment. This outcome serves as the true endpoint.}
    \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = placebo and \eqn{1} = interferon-\eqn{\alpha}.}
  }
}
\keyword{ARMD.MultS}
