\name{Sim.Data.STS}
\alias{Sim.Data.STS}
\title{Simulates a dataset that can be used to assess surrogacy in the single-trial setting}
\description{
The function \code{Sim.Data.STS} simulates a dataset that contains the variables Treat, Surr, True, and Pat.ID. The user can specify the total number of patients, the desired \eqn{R_{indiv}} value (also referred to as the adjusted association (\eqn{\gamma}) in the single-trial meta-analytic setting), and the desired means of the surrogate and the true endpoints in the experimental and control treatment groups.
}

\usage{Sim.Data.STS(N.Total=2000, R.Indiv.Target=.8, Means=c(0, 0, 0, 0), Seed=
sample(1:1000, size=1))}
\arguments{
  \item{N.Total}{The total number of patients in the simulated dataset. Default \eqn{2000}.}
  \item{R.Indiv.Target}{The desired \eqn{R_{indiv}} (or \eqn{\gamma}) value in the simulated dataset. Default \eqn{0.80}.}
  \item{Means}{A vector that specifies the desired mean for the surrogate in the control treatment group, mean for the surrogate in the experimental treatment group, mean for the true endpoint in the control treatment group, and mean for the true endpoint in the experimental treatment group, respectively. Default \code{c(0, 0, 0, 0)}.}
  \item{Seed}{The seed that is used to generate the dataset. Default \code{sample(x=1:1000, size=1)}, i.e., a random number between 1 and 1000.}
}
\details{
The generated object \code{Data.Observed.STS} (of class \code{data.frame}) is placed in the workspace (for easy access).}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{Sim.Data.MTS}}, \code{\link{Single.Trial.RE.AA}}
}
\examples{
# Simulate a dataset: 
Sim.Data.STS(N.Total=2000, R.Indiv.Target=.8, Means=c(1, 5, 20, 37), Seed=1) 
}
\keyword{Simulate data}
\keyword{Single-trial setting}
