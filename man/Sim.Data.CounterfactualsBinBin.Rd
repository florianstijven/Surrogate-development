\name{Sim.Data.CounterfactualsBinBin}
\alias{Sim.Data.CounterfactualsBinBin}
\title{Simulate a dataset that contains counterfactuals for binary endpoints}
\description{
The function \code{Sim.Data.CounterfactualsBinBin} simulates a dataset that contains four (binary) counterfactuals (i.e., potential outcomes) and a (binary) treatment indicator. The counterfactuals \eqn{T_0} and \eqn{T_1} denote the true endpoints of a patient under the control and the experimental treatments, respectively, and the counterfactuals \eqn{S_0} and \eqn{S_1} denote the surrogate endpoints of the patient under the control and the experimental treatments, respectively. The user can specify the number of patients and the desired probabilities of the vector of potential outcomes (i.e., \eqn{\bold{{Y'}_c}}=(T_0, T_1, S_0, S_1)).
}  
\usage{Sim.Data.CounterfactualsBinBin(Pi_s=rep(1/16, 16), 
N.Total=2000, Seed=sample(1:1000, size=1))
}
\arguments{
  \item{Pi_s}{The vector of probabilities of the potential outcomes, i.e., \eqn{pi_{0000}}, \eqn{pi_{0100}}, \eqn{pi_{0010}}, \eqn{pi_{0001}}, \eqn{pi_{0101}}, \eqn{pi_{1000}}, \eqn{pi_{1010}}, \eqn{pi_{1001}}, \eqn{pi_{1110}}, \eqn{pi_{1101}}, \eqn{pi_{1011}}, \eqn{pi_{1111}}, \eqn{pi_{0110}}, \eqn{pi_{0011}}, \eqn{pi_{0111}}, \eqn{pi_{1100}}. Default \code{rep(1/16, 16)}.}
  \item{N.Total}{The desired number of patients in the simulated dataset. Default \eqn{2000}.}
  \item{Seed}{A seed that is used to generate the dataset. Default \code{sample(x=1:1000, size=1)}, i.e., a random number between 1 and 1000.}
  }
\details{
The generated object \code{Data.STSBinBin.Counter} (which contains the counterfactuals) and \code{Data.STSBinBin.Obs} (the "observable data") (of class \code{data.frame}) is placed in the workspace.
}

\value{An object of class \code{Sim.Data.CounterfactualsBinBin} with components,
  \item{Data.STSBinBin.Obs}{The generated dataset that contains the "observed" surrogate endrpoint, true endpoint, and assigned treatment.}
  \item{Data.STSBinBin.Counter}{The generated dataset that contains the counterfactuals.}
    \item{Vector_Pi}{The vector of probabilities of the potential outcomes, i.e., \eqn{pi_{0000}}, \eqn{pi_{0100}}, \eqn{pi_{0010}}, \eqn{pi_{0001}}, \eqn{pi_{0101}}, \eqn{pi_{1000}}, \eqn{pi_{1010}}, \eqn{pi_{1001}}, \eqn{pi_{1110}}, \eqn{pi_{1101}}, \eqn{pi_{1011}}, \eqn{pi_{1111}}, \eqn{pi_{0110}}, \eqn{pi_{0011}}, \eqn{pi_{0111}}, \eqn{pi_{1100}}.}
  \item{Pi_Marginals}{The vector of marginal probabilities \eqn{\pi_{1 \cdot 1 \cdot}}, \eqn{\pi_{0 \cdot 1 \cdot}}, \eqn{\pi_{1 \cdot 0 \cdot}}, \eqn{\pi_{0 \cdot 0 \cdot}}, \eqn{\pi_{\cdot 1 \cdot 1}}, \eqn{\pi_{\cdot 1 \cdot 0}}, \eqn{\pi_{\cdot 0 \cdot 1}}, \eqn{\pi_{\cdot 0 \cdot 0}}.}
 \item{True.C3}{The true \eqn{C_3} value.}
 \item{True.R2_H}{The true \eqn{R_H^2} value.}
 \item{True.Theta_T}{The true odds ratio for \eqn{T}.}
 \item{True.Theta_S}{The true odds ratio for \eqn{S}.}
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs}
\examples{
## Generate a dataset with 2000 patients, and values 1/16
## for all proabilities between the counterfactuals:
Sim.Data.CounterfactualsBinBin(N.Total=2000)
}
\keyword{Causal-Inference framework}
\keyword{Simulate data}
\keyword{Counterfactuals}
\keyword{Binary Binary setting}