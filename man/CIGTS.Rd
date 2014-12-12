\name{CIGTS}
\alias{CIGTS}
\docType{data}
\title{Data of the Collaborative Initial Glaucoma Treatment Study}
\description{
These are the data of the Collaborative Initial Glaucoma Treatment Study (CIGTS), a randomized clinical trial designed to compare the effects of surgery and medicine on intraocular pressure (IOP). Elevated IOP is an important risk factor for glaucoma, which is an occular disorder that leads to impaired vision. S and T are IOP at months \eqn{12} and \eqn{96}, respectively. A total of \eqn{N=228} patients had IOP measured at both time points of interest. The continuous S and T were dichotomized in the following way: \eqn{S=1} if IOP at month \eqn{12} is \eqn{<18} mmHg and \eqn{S=0} if IOP at month \eqn{12} is \eqn{>=18} mmHg (and similarly for T).

Acknowledgment

David Musch (Coordinating Center Director) and Brenda Gillespie (Study Statistician) are gratefully acknowledged for providing data from the Collaborative Initial Glaucoma Treatment Study (CIGTS).}

\usage{data(CIGTS)}
\format{A \code{data.frame} with \eqn{228} observations on \eqn{4} variables.
  \describe{
    \item{\code{Id}}{The Patient ID.}
    \item{\code{Treat}}{The treatment indicator, coded as \eqn{-1} = medical treatment and \eqn{1} = surgery.}
    \item{\code{IOP_12}}{The dichotomized IOP level \eqn{12} months following randomization.}
    \item{\code{IOP_96}}{The dichotomized IOP level \eqn{96} months following randomization.}
  }
}

\references{
Musch, D. C., Lichter, P., R., Guire, K. E., Standardi, C. L., & CIGTS Investigators (1999). The Collaborative Initial Glaucoma Treatment Study (CIGTS): Study design, methods, and baseline characteristics of enrolled patients. Ophthalmology, 106, 653-662.
}

\keyword{CIGTS}