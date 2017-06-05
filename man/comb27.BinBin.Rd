\name{comb27.BinBin}
\alias{comb27.BinBin}
\title{Assesses the surrogate predictive value of each of the 27 prediction functions in the setting where both \eqn{S} and \eqn{T} are binary endpoints}
\description{
The function \code{comb27.BinBin} assesses a surrogate predictive value of each of the 27 possible prediction functions in the single-trial causal-inference framework when both the surrogate and the true endpoints are binary outcomes. The distribution of frequencies at which each of the 27 possible predicton functions are selected provides additional insights regarding the association between \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T}).  See \bold{Details} below.}
\usage{
comb27.BinBin(pi1_1_, pi1_0_, pi_1_1, pi_1_0, 
pi0_1_, pi_0_1, Monotonicity=c("No"),M=1000, Seed=1)
}
\arguments{
  \item{pi1_1_}{A scalar that contains values for \eqn{P(T=1,S=1|Z=0)}, i.e., the probability that \eqn{S=T=1} when under treatment \eqn{Z=0}.}
  \item{pi1_0_}{A scalar that contains values for \eqn{P(T=1,S=0|Z=0)}.}
  \item{pi_1_1}{A scalar that contains values for \eqn{P(T=1,S=1|Z=1)}.}
  \item{pi_1_0}{A scalar that contains values for \eqn{P(T=1,S=0|Z=1)}.}
  \item{pi0_1_}{A scalar that contains values for \eqn{P(T=0,S=1|Z=0)}.}
  \item{pi_0_1}{A scalar that contains values for \eqn{P(T=0,S=1|Z=1)}.}
  \item{Monotonicity}{Specifies which assumptions regarding monotonicity should be made, only one assumption can be made at the time:  \code{Monotonicity=c("No")}, \code{Monotonicity=c("True.Endp")}, \code{Monotonicity=c("Surr.Endp")}, or \code{Monotonicity=c("Surr.True.Endp")}.  Default \code{Monotonicity=c("No")}.}
  \item{M}{The number of random samples that have to be drawn for the freely varying parameters. Default \code{M=100000}. }
  \item{Seed}{The seed to be used to generate \eqn{\pi_r}. Default \code{Seed=1}.}

}

\details{

In the continuous normal setting, surroagacy can be assessed by studying the association between the individual causal effects on \eqn{S} and \eqn{T} (see \code{\link{ICA.ContCont}}). In that setting, the Pearson correlation is the obvious measure of association.

When \eqn{S} and \eqn{T} are binary endpoints, multiple alternatives exist. Alonso et al. (2016) proposed the individual causal association (ICA; \eqn{R_{H}^{2}}), which captures the association between the individual causal effects of the treatment on \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T}) using information-theoretic principles. 

The function \code{comb27.BinBin} computes \eqn{R_{H}^{2}} using a grid-based approach where all possible combinations of the specified grids for the parameters that are allowed to vary freely are considered. It computes the probability of a prediction error for each of the 27 possible prediction functions.The frequency at which each prediction function is selected provides additional insight about the minimal probability of a prediction error PPE which can be obtained with \code{PPE.BinBin}.  
}


\value{An object of class \code{comb27.BinBin} with components,
\item{index}{count variable}
\item{Monotonicity}{The vector of Monotonicity assumptions}
\item{Pe}{The vector of the prediction error values.}
\item{combo}{The vector containing the codes for the each of the 27 prediction functions.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{H_Delta_T}{The vector of the entropies of \eqn{\Delta_T}.}
\item{H_Delta_S}{The vector of the entropies of \eqn{\Delta_S}.}
\item{I_Delta_T_Delta_S}{The vector of the mutual information of \eqn{\Delta_S} and \eqn{\Delta_T}.}


}
\references{
Alonso A,  Van der Elst W, Molenberghs G, Buyse M and Burzykowski T. (2016). An information-theoretic approach for the evaluation of surrogate endpoints based on causal inference. 

Alonso A, Van der Elst W and Meyvisch P (2016). Assessing a surrogate predictive value: A causal inference approach. 
}
\author{
Paul Meyvisch, Wim Van der Elst, Ariel Alonso, Geert Molenberghs
}
\seealso{
\code{\link{PPE.BinBin}}
}
\examples{
# Conduct the analysis assuming no montonicity
 
\dontrun{ # time consuming code part
comb27.BinBin(pi1_1_ = 0.3412, pi1_0_ = 0.2539, pi0_1_ = 0.119, 
              pi_1_1 = 0.6863, pi_1_0 = 0.0882, pi_0_1 = 0.0784,  
              Seed=1,Monotonicity=c("No"), M=500000) 
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{ICA}
\keyword{PPE}