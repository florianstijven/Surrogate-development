\name{Sim.Data.STSBinBin}
\alias{Sim.Data.STSBinBin}
\title{Simulates a dataset that can be used to assess surrogacy in the single trial setting when S and T are binary endpoints}
\description{
The function \code{Sim.Data.STSBinBin} simulates a dataset that contains four (binary) counterfactuals (i.e., potential outcomes) and a (binary) treatment indicator. The counterfactuals \eqn{T_0} and \eqn{T_1} denote the true endpoints of a patient under the control and the experimental treatments, respectively, and the counterfactuals \eqn{S_0} and \eqn{S_1} denote the surrogate endpoints of the patient under the control and the experimental treatments, respectively. 
In addition, the function provides the "observable" data based on the dataset of the counterfactuals, i.e., the \eqn{S} and \eqn{T} endpoints given the treatment that was allocated to a patient. 
The user can specify the assumption regarding monotonicity that should be made to generate the data (no monotonicity, monotonicity for \eqn{S} alone, monotonicity for \eqn{T} alone, or monotonicity for both \eqn{S} and \eqn{T}).
}  
\usage{Sim.Data.STSBinBin(Monotonicity=c("No"), N.Total=2000, Seed)}
\arguments{
  \item{Monotonicity}{The assumption regarding monotonicity that should be made when the data are generated, i.e., \code{Monotonicity="No"} (no monotonicity assumed), \code{Monotonicity="True.Endp"} (monotonicity assumed for the true endpoint alone), \code{Monotonicity="Surr.Endp"} (monotonicity assumed for the surrogate endpoint alone), and \code{Monotonicity="Surr.True.Endp"} (monotonicity assumed for both endpoints). Default \code{Monotonicity="No"}.}
  \item{N.Total}{The desired number of patients in the simulated dataset. Default \eqn{2000}.}
  \item{Seed}{A seed that is used to generate the dataset. Default \code{sample(x=1:1000, size=1)}, i.e., a random number between 1 and 1000.}
  }
\details{
The generated objects \code{Data.STSBinBin_Counterfactuals} (which contains the counterfactuals) and  \code{Data.STSBinBin_Obs} (which contains the observable data) of class \code{data.frame} are placed in the workspace. Other relevant output can be accessed based on the fitted object (see \eqn{Value} below) 
}

\value{An object of class \code{Sim.Data.STSBinBin} with components,
  \item{Data.STSBinBin.Obs}{The generated dataset that contains the "observed" surrogate endrpoint, true endpoint, and assigned treatment.}
  \item{Data.STSBinBin.Counter}{The generated dataset that contains the counterfactuals.}
    \item{Vector_Pi}{The vector of probabilities of the potential outcomes, i.e., \eqn{pi_{0000}}, \eqn{pi_{0100}}, \eqn{pi_{0010}}, \eqn{pi_{0001}}, \eqn{pi_{0101}}, \eqn{pi_{1000}}, \eqn{pi_{1010}}, \eqn{pi_{1001}}, \eqn{pi_{1110}}, \eqn{pi_{1101}}, \eqn{pi_{1011}}, \eqn{pi_{1111}}, \eqn{pi_{0110}}, \eqn{pi_{0011}}, \eqn{pi_{0111}}, \eqn{pi_{1100}}.}
  \item{Pi_Marginals}{The vector of marginal probabilities \eqn{\pi_{1 \cdot 1 \cdot}}, \eqn{\pi_{0 \cdot 1 \cdot}}, \eqn{\pi_{1 \cdot 0 \cdot}}, \eqn{\pi_{0 \cdot 0 \cdot}}, \eqn{\pi_{\cdot 1 \cdot 1}}, \eqn{\pi_{\cdot 1 \cdot 0}}, \eqn{\pi_{\cdot 0 \cdot 1}}, \eqn{\pi_{\cdot 0 \cdot 0}}.}
 \item{True.R2_H}{The true \eqn{R_H^2} value.}
 \item{True.Theta_T}{The true odds ratio for \eqn{T}.}
 \item{True.Theta_S}{The true odds ratio for \eqn{S}.}
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs}
\examples{
## Generate a dataset with 2000 patients, 
## assuming no monotonicity:
Sim.Data.STSBinBin(Monotonicity=c("No"), N.Total=200)
}
\keyword{Causal-Inference framework}
\keyword{Simulate data}
\keyword{Counterfactuals}
\keyword{Binary Binary setting}