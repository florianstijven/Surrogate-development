\name{ICA.BinBin}
\alias{ICA.BinBin}
\title{Assess surrogacy in the causal-inference single-trial setting in the binary-binary case}
\description{
The function \code{ICA.BinBin} quantifies surrogacy in the single-trial causal-inference framework (individual causal association and causal concordance) when both the surrogate and the true endpoints are binary outcomes. See \bold{Details} below.}
\usage{
ICA.BinBin(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1, 
Monotonicity=c("General"), Sum_Pi_f = seq(from=0.01, to=0.99, by=.01), 
M=10000, Volume.Perc=0, Seed=sample(1:100000, size=1))
}

\arguments{
  \item{pi1_1_}{A scalar that contains \eqn{P(T=1,S=1|Z=0)}, i.e., the proability that \eqn{S=T=1} when under treatment \eqn{Z=0}.}
  \item{pi1_0_}{A scalar that contains \eqn{P(T=1,S=0|Z=0)}.}
  \item{pi_1_1}{A scalar that contains \eqn{P(T=1,S=1|Z=1)}.}
  \item{pi_1_0}{A scalar that contains \eqn{P(T=1,S=0|Z=1)}.}
  \item{pi0_1_}{A scalar that contains \eqn{P(T=0,S=1|Z=0)}.}
  \item{pi_0_1}{A scalar that contains \eqn{P(T=0,S=1|Z=1)}.}
  \item{Monotonicity}{Specifies which assumptions regarding monotonicity should be made: \code{Monotonicity=c("General")}, \code{Monotonicity=c("No")}, \code{Monotonicity=c("True.Endp")}, \code{Monotonicity=c("Surr.Endp")}, or \code{Monotonicity=c("Surr.True.Endp")}. See \bold{Details} below. Default \code{Monotonicity=c("General")}.}
  \item{Sum_Pi_f}{A scalar or vector that specifies the grid of values \eqn{G={g_{1},\: g_{2},\:...,\: g_{k}}} to be considered when the sensitivity analysis is conducted. See \bold{Details} below. Default \code{Sum_Pi_f = seq(from=0.01, to=0.99, by=.01)}.}
  \item{M}{The number of runs that are conducted for a given value of \code{Sum_Pi_f}. This argument is not used when \code{Volume.Perc=0}. Default \code{M=10000}.}
  \item{Volume.Perc}{Note that the marginals that are observable in the data set a number of restrictions on the unidentified correlations. For example, under montonicity for \eqn{S} and \eqn{T}, it holds that \eqn{\pi_{0111}<=min(\pi_{0\cdot1\cdot}, \pi_{\cdot1\cdot1})} and \eqn{\pi_{1100}<=min(\pi_{1\cdot0\cdot}, \pi_{\cdot1\cdot0})}. For example, when \eqn{min(\pi_{0\cdot1\cdot}, \pi_{\cdot1\cdot1})=0.10} and \eqn{min(\pi_{1\cdot0\cdot}, \pi_{\cdot1\cdot0})=0.08}, then all valid  \eqn{\pi_{0111}<=0.10} and  all valid \eqn{\pi_{1100}<=0.08}. The argument \code{Volume.Perc} specifies the fraction  of the 'volume' of the paramater space that is explored. This volume is computed based on the grids G={0, 0.01, ..., maximum possible value for the counterfactual probability at hand}. E.g., in the previous example, the 'volume' of the parameter space would be \eqn{11*9=99}, and when e.g., the argument \code{Volume.Perc=1} is used a total of \eqn{99} runs will be conducted for each given value of \code{Sum_Pi_f}. Notice that when monotonicity is not assumed, relatively high values of \code{Volume.Perc} will lead to a large number of runs and consequently a long analysis time.}
  \item{Seed}{The seed to be used to generate \eqn{\pi_r}. Default \code{Seed=sample(1:100000, size=1)}.}
  }

\details{

In the continuous normal setting, surroagacy can be assessed by studying the association between the individual causal effects on \eqn{S} and \eqn{T} (see \code{\link{ICA.ContCont}}). In that setting, the Pearson correlation is the obvious measure of association.

When \eqn{S} and \eqn{T} are binary endpoints, multiple alternatives exist. Alonso et al. (2014) proposed the individual causal association (ICA; \eqn{R_{H}^{2}}), which captures the association between the individual causal effects of the treatment on \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T}) using information-theoretic principles. 

The function \code{ICA.BinBin} computes \eqn{R_{H}^{2}} based on plausible values of the potential outcomes. Denote by \eqn{\bold{Y}'=(T_0,T_1,S_0,S_1)} the vector of potential outcomes. The vector \eqn{\bold{Y}} can take 16 values and the set of parameters \eqn{\pi_{ijpq}=P(T_0=i,T_1=j,S_0=p,S_1=q)} (with \eqn{i,j,p,q=0/1}) fully characterizes its distribution.

However, the parameters in \eqn{\pi_{ijpq}} are not all functionally independent, e.g., \eqn{1=\pi_{\cdot\cdot\cdot\cdot}}. When no assumptions regarding monotonicity are made, the data impose a total of \eqn{7} restrictions, and thus only \eqn{9} proabilities in \eqn{\pi_{ijpq}} are allowed to vary freely (for details, see Alonso et al., 2014). Based on the data and assuming SUTVA, the marginal probabilites \eqn{\pi_{1 \cdot 1 \cdot}}, \eqn{\pi_{1 \cdot 0 \cdot}}, \eqn{\pi_{\cdot 1 \cdot 1}}, \eqn{\pi_{\cdot 1 \cdot 0}}, \eqn{\pi_{0 \cdot 1 \cdot}}, and \eqn{\pi_{\cdot 0 \cdot 1}} can be computed (by hand or using the function \code{\link{MarginalProbs}}). Define the vector 
\deqn{\bold{b}'=(1, \pi_{1 \cdot 1 \cdot}, \pi_{1 \cdot 0 \cdot}, \pi_{\cdot 1 \cdot 1}, \pi_{\cdot 1 \cdot 0}, \pi_{0 \cdot 1 \cdot}, \pi_{\cdot 0 \cdot 1})} and \eqn{\bold{A}} is a contrast matrix such that the identified restrictions can be written as a system of linear equation
\deqn{\bold{A \pi} = \bold{b}.}

The matrix \eqn{\bold{A}} has rank \eqn{7} and can be partitioned as \eqn{\bold{A=(A_r | A_f)}}, and similarly the vector \eqn{\bold{\pi}} can be partitioned as \eqn{\bold{\pi^{'}=(\pi_r^{'} | \pi_f^{'})}} (where \eqn{f} refers to the submatrix/vector given by the \eqn{9} last columns/components of \eqn{\bold{A/\pi}}). Using these partitions the previous system of linear equations can be rewritten as
\deqn{\bold{A_r \pi_r + A_f \pi_f = b}.}

The following algorithm is used to generate plausible distributions for \eqn{\bold{Y}}. First, select a value of the specified grid of values (specified using \code{Sum_Pi_f} in the function call).  For \eqn{k=1} to \eqn{M} (specified using \code{M} in the function call), generate a vector \eqn{\pi_f} that contains \eqn{9} components that are uniformly sampled from hyperplane subject to the restriction that the sum of the generated components equals \code{Sum_Pi_f} (the function \code{\link{RandVec}}, which uses the \code{randfixedsum} algorithm written by Roger Stafford, is used  to obtain these components). Next, \eqn{\bold{\pi_r=A_r^{-1}(b - A_f \pi_f)}} is computed and the \eqn{\pi_r} vectors where all components are in the \eqn{[0;\:1]} range are retained. This procedure is repeated for each of the \code{Sum_Pi_f} values.  Based on these results, \eqn{R_H^2} is estimated. The obtained values can be used to conduct a sensitivity analysis during the validation exercise.

The previous developments hold when no monotonicity is assumed. When monotonicity for \eqn{S}, \eqn{T}, or for \eqn{S} and \eqn{T} is assumed, some of the probabilities of \eqn{\pi} are zero. For example, when montonicity is assumed for \eqn{T}, then \eqn{P(T_0 <= T_1)=1}, or equivantly, \eqn{\pi_{1000}=\pi_{1010}=\pi_{1001}=\pi_{1011}=0}. When monotonicity is assumed, the procedure described above is modified accordingly (for details, see Alonso et al., 2014). When a general analysis is requested (using \code{Monotonicity=c("General")} in the function call), all settings are considered (no monotonicity, monotonicity for \eqn{S} alone, for \eqn{T} alone, and for both for \eqn{S} and \eqn{T}.)
}

\value{An object of class \code{ICA.BinBin} with components,
\item{Pi.Vectors}{An object of class \code{data.frame} that contains the valid \eqn{\pi} vectors.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{Theta_T}{The vector of odds ratios for \eqn{T}.}
\item{Theta_S}{The vector of odds ratios for \eqn{S}.}
\item{H_Delta_T}{The vector of the entropies of \eqn{\Delta_T}.}
\item{Monotonicity}{The assumption regarding monotonicity that was made.}
\item{Volume.No}{The 'volume' of the parameter space when monotonicity is not assumed. Is only provided when the argument \eqn{Volume.Perc} is used (i.e., when it is not equal to \eqn{0}.}
\item{Volume.T}{The 'volume' of the parameter space when monotonicity for \eqn{T} is assumed. Is only provided when the argument \eqn{Volume.Perc} is used.}
\item{Volume.S}{The 'volume' of the parameter space when monotonicity for \eqn{S} is assumed. Is only provided when the argument \eqn{Volume.Perc} is used.}
\item{Volume.ST}{The 'volume' of the parameter space when monotonicity for \eqn{S} and \eqn{T} is assumed. Is only provided when the argument \eqn{Volume.Perc} is used.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2015). Validation of surrogate endpoints: the binary-binary setting from a causal inference perspective.
}
\author{
Wim Van der Elst, Ariel Alonso, Marc Buyse, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.ContCont}}, \code{\link{MICA.ContCont}}
}
\examples{
\dontrun{# Time consuming code part
# Compute R2_H given the marginals specified as the pi's, making no 
# assumptions regarding monotonicity (general case)
ICA <- ICA.BinBin(pi1_1_=0.2619048, pi1_0_=0.2857143, pi_1_1=0.6372549, 
pi_1_0=0.07843137, pi0_1_=0.1349206, pi_0_1=0.127451, Seed=1, 
Monotonicity=c("General"), Sum_Pi_f = seq(from=0.01, to=.99, by=.01), M=10000)

# obtain plot of the results
plot(ICA, R2_H=TRUE)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{MarginalProbs}

