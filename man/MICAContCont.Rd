\name{MICA.ContCont}
\alias{MICA.ContCont}
\title{Assess surrogacy in the causal-inference multiple-trial setting (Meta-analytic Individual Causal Association; MICA) in the continuous-continuous case}
\description{
The function \code{MICA.ContCont} quantifies surrogacy in the multiple-trial causal-inference framework. See \bold{Details} below.}
\usage{MICA.ContCont(Trial.R, D.aa, D.bb, T0S0, T1S1, T0T0=1, T1T1=1, S0S0=1, S1S1=1, 
T0T1=seq(-1, 1, by=.1), T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), 
S0S1=seq(-1, 1, by=.1))}
\arguments{
\item{Trial.R}{A scalar that specifies the trial-level correlation coefficient (i.e., \eqn{R_{trial}}) that should be used in the computation of \eqn{\rho_{M}}.}
\item{D.aa}{A scalar that specifies the between-trial variance of the treatment effects on the surrogate endpoint (i.e., \eqn{d_{aa}}) that should be used in the computation of \eqn{\rho_{M}}.}
\item{D.bb}{A scalar that specifies the between-trial variance of the treatment effects on the true endpoint (i.e., \eqn{d_{bb}}) that should be used in the computation of \eqn{\rho_{M}}.}
\item{T0S0}{A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{M}}.}
  \item{T1S1}{A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{M}}.}
  \item{T0T0}{A scalar that specifies the variance of the true endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{M}}. Default 1.}
  \item{T1T1}{A scalar that specifies the variance of the true endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{M}}. Default 1.}
  \item{S0S0}{A scalar that specifies the variance of the surrogate endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{M}}. Default 1.}
  \item{S1S1}{A scalar that specifies the variance of the surrogate endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{M}}. Default 1.}
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1 that should be considered in the computation of \eqn{\rho_{M}}. Default \code{seq(-1, 1, by=.1)}, i.e., the values \eqn{-1}, \eqn{-0.9}, \eqn{-0.8}, \dots, \eqn{1}.}
  \item{T0S1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and S1 that should be considered in the computation of \eqn{\rho_{M}}. Default \code{seq(-1, 1, by=.1)}.}
  \item{T1S0}{A scalar or vector that contains the correlation(s) between the counterfactuals T1 and S0 that should be considered in the computation of \eqn{\rho_{M}}. Default \code{seq(-1, 1, by=.1)}.}
  \item{S0S1}{A scalar or vector that contains the correlation(s) between the counterfactuals S0 and S1 that should be considered in the computation of \eqn{\rho_{M}}. Default \code{seq(-1, 1, by=.1)}.}
}
\details{
Based on the causal-inference framework, it is assumed that each subject \emph{j} in trial \emph{i} has four counterfactuals (or potential outcomes), i.e., \eqn{T_{0ij}}, \eqn{T_{1ij}}, \eqn{S_{0ij}}, and \eqn{S_{1ij}}. Let \eqn{T_{0ij}} and \eqn{T_{1ij}} denote the counterfactuals for the true endpoint (\eqn{T}) under the control (\eqn{Z=0}) and the experimental (\eqn{Z=1}) treatments of subject \emph{j} in trial \emph{i}, respectively. Similarly, \eqn{S_{0ij}} and \eqn{S_{1ij}} denote the corresponding counterfactuals for the surrogate endpoint (\eqn{S}) under the control and experimental treatments of subject \emph{j} in trial \emph{i}, respectively. The individual causal effects of \eqn{Z} on \eqn{T} and \eqn{S} for a given subject \emph{j} in trial \emph{i} are then defined as \eqn{\Delta_{T_{ij}}=T_{1ij}-T_{0ij}} and \eqn{\Delta_{S_{ij}}=S_{1ij}-S_{0ij}}, respectively.

In the multiple-trial causal-inference framework, surrogacy can be quantified as the correlation between the individual causal effects of \eqn{Z} on \eqn{S} and \eqn{T} (for details, see Alonso et al., submitted):

\deqn{\rho_{M}=\rho(\Delta_{Tij},\:\Delta_{Sij})=\frac{\sqrt{d_{bb}d_{aa}}R_{trial}+\sqrt{V(\varepsilon_{\Delta Tij})V(\varepsilon_{\Delta Sij})}\rho_{\Delta}}{\sqrt{V(\Delta_{Tij})V(\Delta_{Sij})}},}

where

\deqn{V(\varepsilon_{\Delta Tij})=\sigma_{T_{0}T_{0}}+\sigma_{T_{1}T_{1}}-2\sqrt{\sigma_{T_{0}T_{0}}\sigma_{T_{1}T_{1}}}\rho_{T_{0}T_{1}},}
\deqn{V(\varepsilon_{\Delta Sij})=\sigma_{S_{0}S_{0}}+\sigma_{S_{1}S_{1}}-2\sqrt{\sigma_{S_{0}S_{0}}\sigma_{S_{1}S_{1}}}\rho_{S_{0}S_{1}},}
\deqn{V(\Delta_{Tij})=d_{bb}+\sigma_{T_{0}T_{0}}+\sigma_{T_{1}T_{1}}-2\sqrt{\sigma_{T_{0}T_{0}}\sigma_{T_{1}T_{1}}}\rho_{T_{0}T_{1}},}
\deqn{V(\Delta_{Sij})=d_{aa}+\sigma_{S_{0}S_{0}}+\sigma_{S_{1}S_{1}}-2\sqrt{\sigma_{S_{0}S_{0}}\sigma_{S_{1}S_{1}}}\rho_{S_{0}S_{1}}.}

The correlations between the counterfactuals (i.e., \eqn{\rho_{S_{0}T_{1}}}, \eqn{\rho_{S_{1}T_{0}}}, \eqn{\rho_{T_{0}T_{1}}}, and \eqn{\rho_{S_{0}S_{1}}}) are not identifiable from the data. It is thus warranted to conduct a sensitivity analysis (by considering vectors of possible values for the correlations between the counterfactuals -- rather than point estimates). 

When the user specifies a vector of values that should be considered for one or more of the correlations that are involved in the computation of \eqn{\rho_{M}}, the function \code{MICA.ContCont} constructs all possible matrices that can be formed as based on the specified values, identifies the matrices that are positive definite (i.e., valid correlation matrices), and computes \eqn{\rho_{M}} for each of these matrices. An examination of the vector of the obtained \eqn{\rho_{M}} values allows for a straightforward examination of the impact of different assumptions regarding the correlations between the counterfactuals on the results (see also \code{\link{plot Causal-Inference ContCont}}), and the extent to which proponents of the causal-inference and meta-analytic frameworks will reach the same conclusion with respect to the appropriateness of the candidate surrogate at hand.   

\bold{Notes}{

A single \eqn{\rho_{M}} value is obtained when all correlations in the function call are scalars.

}}

\value{An object of class \code{MICA.ContCont} with components,
\item{Total.Num.Matrices}{An object of class \code{numeric} which contains the total number of matrices that can be formed as based on the user-specified correlations.}
\item{Pos.Def}{A \code{data.frame} that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the \eqn{\rho_{M}} values.}
\item{ICA}{A scalar or vector of the \eqn{\rho_{\Delta}} values.}
\item{MICA}{A scalar or vector of the \eqn{\rho_{M}} values.}
}

\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). On the relationship between the causal-inference and meta-analytic paradigms for the validation of surrogate markers. 
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\section{Warning}{The theory that relates the causal-inference and the meta-analytic frameworks in the multiple-trial setting (as developped in Alonso et al., submitted) assumes that a reduced or semi-reduced modelling approach is used in the meta-analytic framework. Thus \eqn{R_{trial}},  \eqn{d_{aa}} and \eqn{d_{bb}} should be estimated based on a reduced model (i.e., using the \code{Model=c("Reduced")} argument in the functions \code{\link{UnifixedContCont}}, \code{\link{UnimixedContCont}}, \code{\link{BifixedContCont}}, or \code{\link{BimixedContCont}}) or based on a semi-reduced model (i.e., using the \code{Model=c("SemiReduced")} argument in the functions \code{\link{UnifixedContCont}}, \code{\link{UnimixedContCont}}, or \code{\link{BifixedContCont}}). } 
\seealso{
\code{\link{ICA.ContCont}},  \code{\link{MICA.Sample.ContCont}}, \code{\link{plot Causal-Inference ContCont}}, \code{\link{UnifixedContCont}}, \code{\link{UnimixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{BimixedContCont}}
}
\examples{
# Generate the vector of MICA values when R_trial=.8, rho_T0S0=rho_T1S1=.8,
# sigma_T0T0=90, sigma_T1T1=100,sigma_ S0S0=10, sigma_S1S1=15, D.aa=5, D.bb=10,  
# and when the grid of values {0, .2, ..., 1} is considered for the 
# correlations between the counterfactuals:
SurMICA <- MICA.ContCont(Trial.R=.80, D.aa=5, D.bb=10, T0S0=.8, T1S1=.8,
T0T0=90, T1T1=100, S0S0=10, S1S1=15, T0T1=seq(0, 1, by=.2), 
T0S1=seq(0, 1, by=.2), T1S0=seq(0, 1, by=.2), S0S1=seq(0, 1, by=.2))

# Examine and plot the vector of the generated MICA values:
summary(SurMICA)
plot(SurMICA)


# Same analysis, but now assume that D.aa=.5 and D.bb=.1:
SurMICA <- MICA.ContCont(Trial.R=.80, D.aa=.5, D.bb=.1, T0S0=.8, T1S1=.8,
T0T0=90, T1T1=100, S0S0=10, S1S1=15, T0T1=seq(0, 1, by=.2), 
T0S1=seq(0, 1, by=.2), T1S0=seq(0, 1, by=.2), S0S1=seq(0, 1, by=.2))

# Examine and plot the vector of the generated MICA values:
summary(SurMICA)
plot(SurMICA)
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Multiple-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{MICA}
