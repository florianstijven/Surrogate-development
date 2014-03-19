\name{ICA.ContCont}
\alias{ICA.ContCont}
\title{Assess surrogacy in the causal-inference single-trial setting (Individual Causal Association, ICA) in the Continuous-continuous case}
\description{
The function \code{ICA.ContCont} quantifies surrogacy in the single-trial causal-inference framework. See \bold{Details} below.}
\usage{ICA.ContCont(T0S0, T1S1, T0T0=1, T1T1=1, S0S0=1, S1S1=1, T0T1=seq(0, 1, by=.1), 
T0S1=seq(0, 1, by=.1), T1S0=seq(0, 1, by=.1), S0S1=seq(0, 1, by=.1))}
\arguments{
  \item{T0S0}{A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}.}
  \item{T1S1}{A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}.}
  \item{T0T0}{A scalar that specifies the variance of the true endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}. Default 1.}
  \item{T1T1}{A scalar that specifies the variance of the true endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}. Default 1.}
  \item{S0S0}{A scalar that specifies the variance of the surrogate endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}. Default 1.}
  \item{S1S1}{A scalar that specifies the variance of the surrogate endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}. Default 1.}
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1 that should be considered in the computation of \eqn{\rho_{\Delta}}. Default \code{seq(0, 1, by=.1)}, i.e., the values \eqn{0}, \eqn{0.10}, \eqn{0.20}, \dots, \eqn{1}.}
  \item{T0S1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and S1 that should be considered in the computation of \eqn{\rho_{\Delta}}. Default \code{seq(0, 1, by=.1)}.}
  \item{T1S0}{A scalar or vector that contains the correlation(s) between the counterfactuals T1 and S0 that should be considered in the computation of \eqn{\rho_{\Delta}}. Default \code{seq(0, 1, by=.1)}.}
  \item{S0S1}{A scalar or vector that contains the correlation(s) between the counterfactuals S0 and S1 that should be considered in the computation of \eqn{\rho_{\Delta}}. Default \code{seq(0, 1, by=.1)}.}
}
\details{
Based on the causal-inference framework, it is assumed that each subject \emph{j} has four counterfactuals (or potential outcomes), i.e., \eqn{T_{0j}}, \eqn{T_{1j}}, \eqn{S_{0j}}, and \eqn{S_{1j}}. Let \eqn{T_{0j}} and \eqn{T_{1j}} denote the counterfactuals for the true endpoint (\eqn{T}) under the control (\eqn{Z=0}) and the experimental (\eqn{Z=1}) treatments of subject \emph{j}, respectively. Similarly, \eqn{S_{0j}} and \eqn{S_{1j}} denote the corresponding counterfactuals for the surrogate endpoint (\eqn{S}) under the control and experimental treatments, respectively. The individual causal effects of \eqn{Z} on \eqn{T} and \eqn{S} for a given subject \emph{j} are then defined as \eqn{\Delta_{T_{j}}=T_{1j}-T_{0j}} and \eqn{\Delta_{S_{j}}=S_{1j}-S_{0j}}, respectively.

In the single-trial causal-inference framework, surrogacy can be quantified as the correlation between the individual causal effects of \eqn{Z} on \eqn{S} and \eqn{T} (for details, see Alonso et al., submitted):

\deqn{\rho_{\Delta}=\rho(\Delta_{T_{j}},\:\Delta_{S_{j}})=\frac{\sqrt{\sigma_{S_{0}S_{0}}\sigma_{T_{0}T_{0}}}\rho_{S_{0}T_{0}}+\sqrt{\sigma_{S_{1}S_{1}}\sigma_{T_{1}T_{1}}}\rho_{S_{1}T_{1}}-\sqrt{\sigma_{S_{0}S_{0}}\sigma_{T_{1}T_{1}}}\rho_{S_{0}T_{1}}-\sqrt{\sigma_{S_{1}S_{1}}\sigma_{T_{0}T_{0}}}\rho_{S_{1}T_{0}}}{\sqrt{(\sigma_{T_{0}T_{0}}+\sigma_{T_{1}T_{1}}-2\sqrt{\sigma_{T_{0}T_{0}}\sigma_{T_{1}T_{1}}}\rho_{T_{0}T_{1}})(\sigma_{S_{0}S_{0}}+\sigma_{S_{1}S_{1}}-2\sqrt{\sigma_{S_{0}S_{0}}\sigma_{S_{1}S_{1}}}\rho_{S_{0}S_{1}})}},} 

where the correlations \eqn{\rho_{S_{0}T_{1}}}, \eqn{\rho_{S_{1}T_{0}}}, \eqn{\rho_{T_{0}T_{1}}}, and \eqn{\rho_{S_{0}S_{1}}} are not estimable. It is thus warranted to conduct a sensitivity analysis (by considering vectors of possible values for the correlations between the counterfactuals -- rather than point estimates). 

When the user specifies a vector of values that should be considered for one or more of the counterfactual correlations in the above expression, the function \code{ICA.ContCont} constructs all possible matrices that can be formed as based on these values, identifies the matrices that are positive definite (i.e., valid correlation matrices), and computes \eqn{\rho_{\Delta}} for each of these matrices. The obtained vector of \eqn{\rho_{\Delta}} values can subsequently be used to examine (i) the impact of different assumptions regarding the correlations between the counterfactuals on the results (see also \code{\link{plot Causal-Inference}}), and (ii) the extent to which proponents of the causal-inference and meta-analytic frameworks will reach the same conclusion with respect to the appropriateness of the candidate surrogate at hand.   

The function \code{ICA} also generates output that is useful to examine the plausibility of finding a good surrogate endpoint (see \code{GoodSurr} in the \bold{Value} section below). For details, see Alonso et al. (submitted).


\bold{Notes}

A single \eqn{\rho_{\Delta}} value is obtained when all correlations in the function call are scalars.

}
\value{An object of class \code{ICA.ContCont} with components,
\item{Total.Num.Matrices}{An object of class \code{numeric} that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.}
\item{Pos.Def}{A \code{data.frame} that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the \eqn{\rho_{\Delta}} values.}
\item{ICA}{A scalar or vector that contains the individual causal association (ICA; \eqn{\rho_{\Delta}}) value(s).}
\item{GoodSurr}{A \code{data.frame} that contains the ICA (\eqn{\rho_{\Delta}}), \eqn{\sigma_{\Delta_{T}}}, and \eqn{\delta}.}  
}
\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). On the relationship between the causal-inference and meta-analytic paradigms for the validation of surrogate markers. 
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{MICA.ContCont}}, \code{\link{Single.Trial.RE.AA}}, \code{\link{plot Causal-Inference}}
}
\examples{
\dontrun{
# Generate the vector of ICA.ContCont values when rho_T0S0=rho_T1S1=.95, 
# sigma_T0T0=90, sigma_T1T1=100,sigma_ S0S0=10, sigma_S1S1=15, and  
# the grid of values {0, .1, ..., 1} is considered for the correlations
# between the counterfactuals:
SurICA <- ICA.ContCont(T0S0=.95, T1S1=.95, T0T0=90, T1T1=100, S0S0=10, S1S1=15,
T0T1=seq(0, 1, by=.1), T0S1=seq(0, 1, by=.1), T1S0=seq(0, 1, by=.1), 
S0S1=seq(0, 1, by=.1))

# Examine and plot the vector of generated ICA values:
summary(SurICA)
\dontrun{plot(SurICA)}

# Obtain the positive definite matrices than can be formed as based on the 
# specified (vectors) of the correlations (these matrices are used to 
# compute the ICA values)
SurICA$Pos.Def}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}

