\name{MarginalProbs}
\alias{MarginalProbs}
\title{Computes marginal probabilities for a dataset where the surrogate and true endpoints are binary}
\description{
This function computes the marginal probabilities associated with the distribution of the potential outcomes for the true and surrogate endpoint.}
\usage{
MarginalProbs(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat)
}
\arguments{
    \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a binary surrogate value, a binary true endpoint value, and a treatment indicator. }
  \item{Surr}{The name of the variable in \code{Dataset} that contains the binary surrogate endpoint values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the binary true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group. }
}
\value{   
\item{Theta_T0S0}{The odds ratio for \eqn{S} and \eqn{T} in the control group.}
\item{Theta_T1S1}{The odds ratio for \eqn{S} and \eqn{T} in the experimental group.}
\item{Freq.Cont}{The frequencies for \eqn{S} and \eqn{T} in the control group.}
\item{Freq.Exp}{The frequencies for \eqn{S} and \eqn{T} in the experimental group.}
\item{pi1_1_}{The estimated \eqn{\pi_{1 \cdot 1 \cdot}}}
\item{pi0_1_}{The estimated \eqn{\pi_{0 \cdot 1 \cdot}}}
\item{pi1_0_}{The estimated \eqn{\pi_{1 \cdot 0 \cdot}}}
\item{pi0_0_}{The estimated \eqn{\pi_{0 \cdot 0 \cdot}}}
\item{pi_1_1}{The estimated \eqn{\pi_{\cdot 1 \cdot 1}}}
\item{pi_1_0}{The estimated \eqn{\pi_{\cdot 1 \cdot 0}}}
\item{pi_0_1}{The estimated \eqn{\pi_{\cdot 0 \cdot 1}}}
\item{pi_0_0}{The estimated \eqn{\pi_{\cdot 0 \cdot 0}}}
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinBin}}
}
\examples{
# Open the ARMD dataset and recode Diff24 and Diff52 as 1
# when the original value is above 0, and 0 otherwise
data(ARMD)
ARMD$Diff24_Dich <- ifelse(ARMD$Diff24>0, 1, 0)
ARMD$Diff52_Dich <- ifelse(ARMD$Diff52>0, 1, 0)

# Obtain marginal probabilities and ORs
MarginalProbs(Dataset=ARMD, Surr=Diff24_Dich, True=Diff52_Dich, 
Treat=Treat)



}
\keyword{Marginal probabilities}
