\name{ICA.BinBin.Grid.Sample}
\alias{ICA.BinBin.Grid.Sample}
\title{Assess surrogacy in the causal-inference single-trial setting in the binary-binary case when monotonicity for \eqn{S} and \eqn{T} is assumed using the grid-based sample approach}
\description{
The function \code{ICA.BinBin.Grid.Sample} quantifies surrogacy in the single-trial causal-inference framework (individual causal association and causal concordance) when both the surrogate and the true endpoints are binary outcomes. This method provides an alternative for \code{ICA.BinBin} and \code{ICA.BinBin.Grid.Full}. It  uses an alternative strategy to identify plausible values for \eqn{\pi}.  See \bold{Details} below.}
\usage{
ICA.BinBin.Grid.Sample(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, 
pi_0_1, Monotonicity=c("General"), M=100000, 
Seed=sample(1:100000, size=1))
}
\arguments{
  \item{pi1_1_}{A scalar that contains \eqn{P(T=1,S=1|Z=0)}, i.e., the proability that \eqn{S=T=1} when under treatment \eqn{Z=0}.}
  \item{pi1_0_}{A scalar that contains \eqn{P(T=1,S=0|Z=0)}.}
  \item{pi_1_1}{A scalar that contains \eqn{P(T=1,S=1|Z=1)}.}
  \item{pi_1_0}{A scalar that contains \eqn{P(T=1,S=0|Z=1)}.}
  \item{pi0_1_}{A scalar that contains \eqn{P(T=0,S=1|Z=0)}.}
  \item{pi_0_1}{A scalar that contains \eqn{P(T=0,S=1|Z=1)}.}
  \item{Monotonicity}{Specifies which assumptions regarding monotonicity should be made: \code{Monotonicity=c("General")}, \code{Monotonicity=c("No")}, \code{Monotonicity=c("True.Endp")}, \code{Monotonicity=c("Surr.Endp")}, or \code{Monotonicity=c("Surr.True.Endp")}. See \bold{Details} below. Default \code{Monotonicity=c("General")}.}
  \item{M}{The number of random samples that have to be drawn for the freely varying parameters. Default \code{pi_1100=seq(0, 1, by=.02)}.}
  \item{Seed}{The seed to be used to generate \eqn{\pi_r}. Default \code{M=100000}.}

}

\details{

In the continuous normal setting, surroagacy can be assessed by studying the association between the individual causal effects on \eqn{S} and \eqn{T} (see \code{\link{ICA.ContCont}}). In that setting, the Pearson correlation is the obvious measure of association.

When \eqn{S} and \eqn{T} are binary endpoints, multiple alternatives exist. Alonso et al. (2014) proposed the individual causal association (ICA; \eqn{R_{H}^{2}}), which captures the association between the individual causal effects of the treatment on \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T}) using information-theoretic principles. 

The function \code{ICA.BinBin.Grid.Full} computes \eqn{R_{H}^{2}} using a grid-based approach where all possible combinations of the specified grids for the parameters that are allowed that are allowed to vary freely are considered.  When it is not assumed that monotonicity holds for both \eqn{S} and \eqn{T}, the number of possible combinations become very high. The function \code{ICA.BinBin.Grid.Sample} considers a random sample of all possible combinations. 
}

\value{An object of class \code{ICA.BinBin} with components,
\item{Pi.Vectors}{An object of class \code{data.frame} that contains the valid \eqn{\pi} vectors.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{Theta_T}{The vector of odds ratios for \eqn{T}.}
\item{Theta_S}{The vector of odds ratios for \eqn{S}.}
\item{H_Delta_T}{The vector of the entropies of \eqn{\Delta_T}.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2014). Validation of surrogate endpoints: the binary-binary setting from a causal inference perspective.

Buyse, M., Burzykowski, T., Aloso, A., & Molenberghs, G. (2014). Direct estimation of joint counterfactual probabilities, with application to surrogate marker validation.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.ContCont}}, \code{\link{MICA.ContCont}}, \code{\link{ICA.BinBin}}, \code{\link{ICA.BinBin.Grid.Sample}}
}
\examples{
\dontrun{
# Compute R2_H given the marginals, 
# assuming monotonicity for S and T and grids
# pi_0111=seq(0, 1, by=.001) and 
# pi_1100=seq(0, 1, by=.001)
ICA <- ICA.BinBin.Grid.Sample(pi1_1_=0.261, pi1_0_=0.285, 
pi_1_1=0.637, pi_1_0=0.078, pi0_1_=0.134, pi_0_1=0.127,  
Monotonicity=c("General"), M=10000, Seed=1)

# obtain plot of R2_H
plot(ICA, R2_H=TRUE)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{MarginalProbs}

