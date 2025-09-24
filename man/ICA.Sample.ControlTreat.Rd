\name{ICA.Sample.ControlTreat}
\alias{ICA.Sample.ControlTreat}
\title{Assess surrogacy in the causal-inference single-trial setting (Individual Causal Association, ICA) in the Continuous-continuous case using the grid-based sample approach when data is only avalable for the control treatment}
\description{
The function \code{ICA.Sample.ControlTreat} quantifies surrogacy in the single-trial causal-inference framework when data is only avalable for the control treatment. }
\usage{ICA.Sample.ControlTreat(T0S0, T1S1=seq(-1, 1, by = 0.001), 
T0T0=1, T1T1=1, S0S0=1, S1S1=1, T0T1=seq(-1, 1, by=.001), 
T0S1=seq(-1, 1, by=.001), T1S0=seq(-1, 1, by=.001), S0S1=seq(-1, 1, by=.001), 
M=50000, M.Target=NA)}

\arguments{
  \item{T0S0}{A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}.}
  \item{T1S1}{A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}.}
  \item{T0T0}{A scalar or vector that specifies the variance of the true endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}. Default 1.}
  \item{T1T1}{A scalar or vector that specifies the variance of the true endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}. Default 1.}
  \item{S0S0}{A scalar or vector that specifies the variance of the surrogate endpoint in the control treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}. Default 1.}
  \item{S1S1}{A scalar or vector that specifies the variance of the surrogate endpoint in the experimental treatment condition that should be considered in the computation of \eqn{\rho_{\Delta}}. Default 1.}
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1 that should be considered in the computation of \eqn{\rho_{\Delta}}. Default \code{seq(-1, 1, by=.001)}.}
  \item{T0S1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and S1 that should be considered in the computation of \eqn{\rho_{\Delta}}. Default \code{seq(-1, 1, by=.001)}.}
  \item{T1S0}{A scalar or vector that contains the correlation(s) between the counterfactuals T1 and S0 that should be considered in the computation of \eqn{\rho_{\Delta}}. Default \code{seq(-1, 1, by=.001)}.}
  \item{S0S1}{A scalar or vector that contains the correlation(s) between the counterfactuals S0 and S1 that should be considered in the computation of \eqn{\rho_{\Delta}}. Default \code{seq(-1, 1, by=.001)}.}
  \item{M}{The number of runs that should be conducted. Default \code{50000}.}
    \item{M.Target}{The number of ICA values that should be identified. Only one argument \code{M=} or \code{M.Target=} can be used.}
}
\value{An object of class \code{ICA.ContCont} with components,
\item{Total.Num.Matrices}{An object of class \code{numeric} that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.}
\item{Pos.Def}{A \code{data.frame} that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the \eqn{\rho_{\Delta}} values.}
\item{ICA}{A scalar or vector that contains the individual causal association (ICA; \eqn{\rho_{\Delta}}) value(s).}
\item{GoodSurr}{A \code{data.frame} that contains the ICA (\eqn{\rho_{\Delta}}), \eqn{\sigma_{\Delta_{T}}}, and \eqn{\delta}.}  
\item{Variances}{A \code{data.frame} that contains the variances for S and T in both treatment conditions.}

}
\references{
Van der Elst, W. et al. (submitted). On the Early Identification of Promising Surrogate Endpoints Using Causal Inference. 
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{MICA.ContCont}}, \code{\link{ICA.ContCont}}, \code{\link{Single.Trial.RE.AA}}, 
\code{\link{plot Causal-Inference ContCont}}
}
\examples{
# Generate the vector of ICA values when rho_T0S0=.95, 
# sigma_T0T0=90, sigma_T1T1=100,sigma_ S0S0=10, sigma_S1S1=15, and  
# min=-1 max=1 is considered for the correlations
# between the counterfactuals and rho_T1S1:
SurICA2 <- ICA.Sample.ControlTreat(T0S0=.95, T0T0=90, T1T1=100, S0S0=10, 
S1S1=15, M=5000)

# Examine and plot the vector of generated ICA values:
summary(SurICA2)
plot(SurICA2)
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}

