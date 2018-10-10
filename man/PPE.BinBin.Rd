\name{PPE.BinBin}
\alias{PPE.BinBin}
\title{Evaluate a surrogate predictive value based on the minimum probability of a prediction error in the setting where both \eqn{S} and \eqn{T} are binary endpoints}
\description{
The function \code{PPE.BinBin} assesses a surrogate predictive value using the probability of a prediction error in the single-trial causal-inference framework when both the surrogate and the true endpoints are binary outcomes. It additionally assesses the indivdiual causal association (ICA). See \bold{Details} below.}
\usage{
PPE.BinBin(pi1_1_, pi1_0_, pi_1_1, pi_1_0, 
pi0_1_, pi_0_1, M=10000, Seed=1)
}
\arguments{
  \item{pi1_1_}{A scalar that contains values for \eqn{P(T=1,S=1|Z=0)}, i.e., the probability that \eqn{S=T=1} when under treatment \eqn{Z=0}.}
  \item{pi1_0_}{A scalar that contains values for \eqn{P(T=1,S=0|Z=0)}.}
  \item{pi_1_1}{A scalar that contains values for \eqn{P(T=1,S=1|Z=1)}.}
  \item{pi_1_0}{A scalar that contains values for \eqn{P(T=1,S=0|Z=1)}.}
  \item{pi0_1_}{A scalar that contains values for \eqn{P(T=0,S=1|Z=0)}.}
  \item{pi_0_1}{A scalar that contains values for \eqn{P(T=0,S=1|Z=1)}.}
  \item{M}{The number of valid vectors that have to be obtained. Default \code{M=10000}. }
  \item{Seed}{The seed to be used to generate \eqn{\pi_r}. Default \code{Seed=1}.}

}

\details{

In the continuous normal setting, surroagacy can be assessed by studying the association between the individual causal effects on \eqn{S} and \eqn{T} (see \code{\link{ICA.ContCont}}). In that setting, the Pearson correlation is the obvious measure of association.

When \eqn{S} and \eqn{T} are binary endpoints, multiple alternatives exist. Alonso et al. (2016) proposed the individual causal association (ICA; \eqn{R_{H}^{2}}), which captures the association between the individual causal effects of the treatment on \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T}) using information-theoretic principles. 

The function \code{PPE.BinBin} computes \eqn{R_{H}^{2}} using a grid-based approach where all possible combinations of the specified grids for the parameters that are allowed to vary freely are considered. It additionally computes the minimal probability of a prediction error (PPE) and the reduction on the PPE using information that \eqn{S} conveys on \eqn{T}. Both measures provide complementary information over the \eqn{R_{H}^{2}} and facilitate more straightforward clinical interpretation. No assumption about monotonicity can be made. 
}


\value{An object of class \code{PPE.BinBin} with components,
\item{index}{count variable}
\item{PPE}{The vector of the PPE values.}
\item{RPE}{The vector of the RPE values.}
\item{PPE_T}{The vector of the \eqn{PPE_T} values indicating the probability on a prediction error without using information on \eqn{S}.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{H_Delta_T}{The vector of the entropies of \eqn{\Delta_T}.}
\item{H_Delta_S}{The vector of the entropies of \eqn{\Delta_S}.}
\item{I_Delta_T_Delta_S}{The vector of the mutual information of \eqn{\Delta_S} and \eqn{\Delta_T}.}


}
\references{
Alonso A,  Van der Elst W, Molenberghs G, Buyse M and Burzykowski T. (2016). An information-theoretic approach for the evaluation of surrogate endpoints based on causal inference. 

Meyvisch P., Alonso A.,Van der Elst W, Molenberghs G. (2018). Assessing the predictive value of a binary surrogate for a binary true endpoint, based on the minimum probability of a prediction error.
}
\author{
Paul Meyvisch, Wim Van der Elst, Ariel Alonso, Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinBin.Grid.Sample}}
}
\examples{
# Conduct the analysis 
 
\dontrun{ # time consuming code part
PPE.BinBin(pi1_1_=0.4215, pi0_1_=0.0538, pi1_0_=0.0538,
           pi_1_1=0.5088, pi_1_0=0.0307,pi_0_1=0.0482, 
           Seed=1, M=10000) 
}

}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{ICA}
\keyword{RPE}
\keyword{PPE}