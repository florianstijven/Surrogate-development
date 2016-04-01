\name{ICA.BinBin.Grid.Sample.Uncert}
\alias{ICA.BinBin.Grid.Sample.Uncert}
\title{Assess surrogacy in the causal-inference single-trial setting in the binary-binary case when monotonicity for \eqn{S} and \eqn{T} is assumed using the grid-based sample approach, accounting for sampling variability in the marginal \eqn{\pi}.}
\description{
The function \code{ICA.BinBin.Grid.Sample.Uncert} quantifies surrogacy in the single-trial causal-inference framework (individual causal association and causal concordance) when both the surrogate and the true endpoints are binary outcomes. This method provides an alternative for \code{ICA.BinBin} and \code{ICA.BinBin.Grid.Full}. It  uses an alternative strategy to identify plausible values for \eqn{\pi}. The function allows to account for sampling variability in the marginal \eqn{\pi}. See \bold{Details} below.}
\usage{
ICA.BinBin.Grid.Sample.Uncert(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, 
pi_0_1, Monotonicity=c("General"), M=100000, 
Volume.Perc=0, Seed=sample(1:100000, size=1))
}
\arguments{
  \item{pi1_1_}{A vector that contains values for \eqn{P(T=1,S=1|Z=0)}, i.e., the probability that \eqn{S=T=1} when under treatment \eqn{Z=0}. A vector is specified to account for uncertainty, i.e., rather than keeping \eqn{P(T=1,S=1|Z=0)} fixed at one estimated value, a distribution can be specified (see \bold{examples} below) from which a value is drawn in each run.}
  \item{pi1_0_}{A vector that contains values for \eqn{P(T=1,S=0|Z=0)}.}
  \item{pi_1_1}{A vector that contains values for \eqn{P(T=1,S=1|Z=1)}.}
  \item{pi_1_0}{A vector that contains values for \eqn{P(T=1,S=0|Z=1)}.}
  \item{pi0_1_}{A vector that contains values for \eqn{P(T=0,S=1|Z=0)}.}
  \item{pi_0_1}{A vector that contains values for \eqn{P(T=0,S=1|Z=1)}.}
  \item{Monotonicity}{Specifies which assumptions regarding monotonicity should be made: \code{Monotonicity=c("General")}, \code{Monotonicity=c("No")}, \code{Monotonicity=c("True.Endp")}, \code{Monotonicity=c("Surr.Endp")}, or \code{Monotonicity=c("Surr.True.Endp")}. When a general analysis is requested (using \code{Monotonicity=c("General")} in the function call), all settings are considered (no monotonicity, monotonicity for \eqn{S} alone, for \eqn{T} alone, and for both for \eqn{S} and \eqn{T}. Default \code{Monotonicity=c("General")}.}
  \item{M}{The number of random samples that have to be drawn for the freely varying parameters. Default \code{M=100000}. This argument is not used when \code{Volume.Perc=0}. Default \code{M=10000}.}
    \item{Volume.Perc}{Note that the marginals that are observable in the data set a number of restrictions on the unidentified correlations. For example, under montonicity for \eqn{S} and \eqn{T}, it holds that \eqn{\pi_{0111}<=min(\pi_{0\cdot1\cdot}, \pi_{\cdot1\cdot1})} and \eqn{\pi_{1100}<=min(\pi_{1\cdot0\cdot}, \pi_{\cdot1\cdot0})}. For example, when \eqn{min(\pi_{0\cdot1\cdot}, \pi_{\cdot1\cdot1})=0.10} and \eqn{min(\pi_{1\cdot0\cdot}, \pi_{\cdot1\cdot0})=0.08}, then all valid  \eqn{\pi_{0111}<=0.10} and  all valid \eqn{\pi_{1100}<=0.08}. The argument \code{Volume.Perc} specifies the fraction  of the 'volume' of the paramater space that is explored. This volume is computed based on the grids G={0, 0.01, ..., maximum possible value for the counterfactual probability at hand}. E.g., in the previous example, the 'volume' of the parameter space would be \eqn{11*9=99}, and when e.g., the argument \code{Volume.Perc=1} is used a total of \eqn{99} runs will be conducted. Notice that when monotonicity is not assumed, relatively high values of \code{Volume.Perc} will lead to a large number of runs and consequently a long analysis time.}
  \item{Seed}{The seed to be used to generate \eqn{\pi_r}. Default \code{M=100000}.}

}

\details{

In the continuous normal setting, surroagacy can be assessed by studying the association between the individual causal effects on \eqn{S} and \eqn{T} (see \code{\link{ICA.ContCont}}). In that setting, the Pearson correlation is the obvious measure of association.

When \eqn{S} and \eqn{T} are binary endpoints, multiple alternatives exist. Alonso et al. (2014) proposed the individual causal association (ICA; \eqn{R_{H}^{2}}), which captures the association between the individual causal effects of the treatment on \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T}) using information-theoretic principles. 

The function \code{ICA.BinBin.Grid.Full} computes \eqn{R_{H}^{2}} using a grid-based approach where all possible combinations of the specified grids for the parameters that are allowed that are allowed to vary freely are considered.  When it is not assumed that monotonicity holds for both \eqn{S} and \eqn{T}, the number of possible combinations become very high. The function \code{ICA.BinBin.Grid.Sample.Uncert} considers a random sample of all possible combinations. 
}

\value{An object of class \code{ICA.BinBin} with components,
\item{Pi.Vectors}{An object of class \code{data.frame} that contains the valid \eqn{\pi} vectors.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{Theta_T}{The vector of odds ratios for \eqn{T}.}
\item{Theta_S}{The vector of odds ratios for \eqn{S}.}
\item{H_Delta_T}{The vector of the entropies of \eqn{\Delta_T}.}
\item{Volume.No}{The 'volume' of the parameter space when monotonicity is not assumed.}
\item{Volume.T}{The 'volume' of the parameter space when monotonicity for \eqn{T} is assumed.}
\item{Volume.S}{The 'volume' of the parameter space when monotonicity for \eqn{S} is assumed.}
\item{Volume.ST}{The 'volume' of the parameter space when monotonicity for \eqn{S} and \eqn{T} is assumed.}

}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2014). Validation of surrogate endpoints: the binary-binary setting from a causal inference perspective.

Buyse, M., Burzykowski, T., Aloso, A., & Molenberghs, G. (2014). Direct estimation of joint counterfactual probabilities, with application to surrogate marker validation.
}
\author{
Wim Van der Elst, Paul Meyvisch, Ariel Alonso & Geert Molenberghs
}
\seealso{
\code{\link{ICA.ContCont}}, \code{\link{MICA.ContCont}}, \code{\link{ICA.BinBin}}, \code{\link{ICA.BinBin.Grid.Sample.Uncert}}
}
\examples{
# Compute R2_H given the marginals (sample from uniform), 
# assuming no monotonicity 
ICA_No2 <- ICA.BinBin.Grid.Sample.Uncert(pi1_1_=runif(10000, 0.3562, 0.4868),
pi0_1_=runif(10000, 0.0240, 0.0837), pi1_0_=runif(10000, 0.0240, 0.0837),
pi_1_1=runif(10000, 0.4434, 0.5742), pi_1_0=runif(10000, 0.0081, 0.0533),
pi_0_1=runif(10000, 0.0202, 0.0763), Seed=1, Monotonicity=c("No"), M=1000)

summary(ICA_No2)

# obtain plot of R2_H
plot(ICA_No2)
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{MarginalProbs}

