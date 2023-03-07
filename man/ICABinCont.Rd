\name{ICA.BinCont}
\alias{ICA.BinCont}
\title{Assess surrogacy in the causal-inference single-trial setting in the binary-continuous case}
\description{
The function \code{ICA.BinCont} quantifies surrogacy in the single-trial setting within the causal-inference framework (individual causal association) when the surrogate endpoint is continuous (normally distributed) and the true endpoint is a binary outcome. For details, see Alonso Abad \emph{et al.} (2023). }
\usage{
ICA.BinCont(Dataset, Surr, True, Treat, 
  BS=FALSE,
  G_pi_10=c(0,1), 
  G_rho_01_00=c(-1,1), 
  G_rho_01_01=c(-1,1), 
  G_rho_01_10=c(-1,1), 
  G_rho_01_11=c(-1,1), 
  Theta.S_0, 
  Theta.S_1, 
  M=1000, Seed=123, 
  Monotonicity=FALSE,
  Independence=FALSE,
  HAA=FALSE,
  Cond_ind=FALSE,
  Plots=TRUE, Save.Plots="No", Show.Details=FALSE)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, and a treatment indicator.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate endpoint values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group.}
  \item{BS}{Logical. If \code{BS=TRUE}, the sampling variability is accounted for in the analysis by using a bootstrap procedure. Default \code{BS=FALSE}.}
  \item{G_pi_10}{The lower and upper limits of the uniform distribution from which the probability parameter \eqn{\pi_{10}} is sampled. Default \code{c(0,1)}. When \code{Monotonicity=TRUE} the values of these limits are set as \code{c(0,0)}.}
  \item{G_rho_01_00}{The lower and upper limits of the uniform distribution from which the association parameter \eqn{\rho_{01}^{00}} is sampled. Default \code{c(-1,1)}.}
  \item{G_rho_01_01}{The lower and upper limits of the uniform distribution from which the association parameter \eqn{\rho_{01}^{01}} is sampled. Default \code{c(-1,1)}.}
  \item{G_rho_01_10}{The lower and upper limits of the uniform distribution from which the association parameter \eqn{\rho_{01}^{10}} is sampled. Default \code{c(-1,1)}.}
  \item{G_rho_01_11}{The lower and upper limits of the uniform distribution from which the association parameter \eqn{\rho_{01}^{11}} is sampled. Default \code{c(-1,1)}.}
  \item{Theta.S_0}{The starting values of the means and standard deviations for the mixture distribution of the surrogate endpoint in the control group. The argument should contain eight values, where the first four values represent the starting values for the means and the last four values represent the starting values for the standard deviations. These starting values should be approximated based on the data on hand. Example: \code{Theta.S_0=c(-10,-5,5,10,10,10,10,10)}.}
  \item{Theta.S_1}{The starting values of the means and standard deviations for the mixture distribution of the surrogate endpoint in the treatment group. The argument should contain eight values, where the first four values represent the starting values for the means and the last four values represent the starting values for the standard deviations. These starting values should be approximated based on the data on hand. Example: \code{Theta.S_1=c(-10,-5,5,10,10,10,10,10)}.}
  \item{M}{The number of Monte Carlo iterations. Default \code{M=1000}.}
  \item{Seed}{The random seed to be used in the analysis (for reproducibility). Default \code{Seed=123}.}
  \item{Monotonicity}{Logical. If \code{Monotonicity=TRUE}, the analysis is performed assuming monotonicity, i.e. \eqn{P(T_1 < T_0) = 0}. Default \code{Monotonicity=FALSE}.}
  \item{Independence}{Logical. If \code{Independence=TRUE}, the analysis is performed assuming independence between the treatment effect in both groups, i.e. \eqn{\pi_{ij} = \pi_{i.} \times \pi_{.j}}. Default \code{Independence=FALSE}.}
  \item{HAA}{Logical. If \code{HAA=TRUE}, the analysis is performed assuming homogeneous association, i.e. \eqn{\rho_{01}^{ij} = \rho_{01}}. Default \code{HAA=FALSE}.}
  \item{Cond_ind}{Logical. If \code{Cond_ind=TRUE}, the analysis is performed assuming conditional independence, i.e. \eqn{\rho_{01} = 0}. Default \code{Cond_ind=FALSE}.}
\item{Plots}{Logical. Should histograms of \eqn{S_0} (surrogate endpoint in control group) and \eqn{S_1} (surrogate endpoint in experimental treatment group) be provided together with density of fitted mixtures? Default \code{Plots=TRUE}.}
\item{Save.Plots}{Should the plots (see previous item) be saved? If \code{Save.Plots="No"}, no plots are saved. If plots have to be saved, replace \code{"No"} by the desired location, e.g., \code{Save.Plots="C:/"}. Default \code{Save.Plots="No"}.} 
\item{Show.Details}{Should some details regarding the availability of some output from the function be displayed in the console when the analysis is running? Setting \code{Show.Details=TRUE} could be useful for debugging procedure (if any). Default \code{Show.Details=FALSE}.}
  }
\value{An object of class \code{ICA.BinCont} with components,
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{pi_00}{The vector of \eqn{\pi_{00}^T} values.}
\item{pi_01}{The vector of \eqn{\pi_{01}^T} values.}
\item{pi_10}{The vector of \eqn{\pi_{10}^T} values.}
\item{pi_11}{The vector of \eqn{\pi_{11}^T} values.}
\item{G_rho_01_00}{The vector of the \eqn{\rho_{01}^{00}} values.}
\item{G_rho_01_01}{The vector of the \eqn{\rho_{01}^{01}} values.}
\item{G_rho_01_10}{The vector of the \eqn{\rho_{01}^{10}} values.}
\item{G_rho_01_11}{The vector of the \eqn{\rho_{01}^{11}} values.}
\item{pi_Delta_T_min1}{The vector of the \eqn{\pi_{-1}^{\Delta T}} values.}
\item{pi_Delta_T_0}{The vector of the \eqn{\pi_{0}^{\Delta T}} values.}
\item{pi_Delta_T_1}{The vector of the \eqn{\pi_{1}^{\Delta T}} values.}
\item{pi_0_00}{The vector of \eqn{\pi_{00}} values of \eqn{f(S_0)}.}
\item{pi_0_01}{The vector of \eqn{\pi_{01}} values of \eqn{f(S_0)}.}
\item{pi_0_10}{The vector of \eqn{\pi_{10}} values of \eqn{f(S_0)}.}
\item{pi_0_11}{The vector of \eqn{\pi_{11}} values of \eqn{f(S_0)}.}
\item{mu_0_00}{The vector of mean \eqn{\mu_{0}^{00}} values of \eqn{f(S_0)}.}
\item{mu_0_01}{The vector of mean \eqn{\mu_{0}^{01}} values of \eqn{f(S_0)}.}
\item{mu_0_10}{The vector of mean \eqn{\mu_{0}^{10}} values of \eqn{f(S_0)}.}
\item{mu_0_11}{The vector of mean \eqn{\mu_{0}^{11}} values of \eqn{f(S_0)}.}
\item{sigma2_00_00}{The vector of variance \eqn{\sigma_{00}^{00}} values of \eqn{f(S_0)}.}
\item{sigma2_00_01}{The vector of variance \eqn{\sigma_{00}^{01}} values of \eqn{f(S_0)}.}
\item{sigma2_00_10}{The vector of variance \eqn{\sigma_{00}^{10}} values of \eqn{f(S_0)}.}
\item{sigma2_00_11}{The vector of variance \eqn{\sigma_{00}^{11}} values of \eqn{f(S_0)}.}
\item{pi_1_00}{The vector of \eqn{\pi_{00}} values of \eqn{f(S_1)}.}
\item{pi_1_01}{The vector of \eqn{\pi_{01}} values of \eqn{f(S_1)}.}
\item{pi_1_10}{The vector of \eqn{\pi_{10}} values of \eqn{f(S_1)}.}
\item{pi_1_11}{The vector of \eqn{\pi_{11}} values of \eqn{f(S_1)}.}
\item{mu_1_00}{The vector of mean \eqn{\mu_{1}^{00}} values of \eqn{f(S_1)}.}
\item{mu_1_01}{The vector of mean \eqn{\mu_{1}^{01}} values of \eqn{f(S_1)}.}
\item{mu_1_10}{The vector of mean \eqn{\mu_{1}^{10}} values of \eqn{f(S_1)}.}
\item{mu_1_11}{The vector of mean \eqn{\mu_{1}^{11}} values of \eqn{f(S_1)}.}
\item{sigma2_11_00}{The vector of variance \eqn{\sigma_{11}^{00}} values of \eqn{f(S_1)}.}
\item{sigma2_11_01}{The vector of variance \eqn{\sigma_{11}^{01}} values of \eqn{f(S_1)}.}
\item{sigma2_11_10}{The vector of variance \eqn{\sigma_{11}^{10}} values of \eqn{f(S_1)}.}
\item{sigma2_11_11}{The vector of variance \eqn{\sigma_{11}^{11}} values of \eqn{f(S_1)}.}
\item{mean_Y_S0}{The vector of mean \eqn{\mu_{0}} values of \eqn{f(S_0)}.}
\item{mean_Y_S1}{The vector of mean \eqn{\mu_{1}} values of \eqn{f(S_1)}.}
\item{var_Y_S0}{The vector of variance \eqn{\sigma_{00}} values of \eqn{f(S_0)}.}
\item{var_Y_S1}{The vector of variance \eqn{\sigma_{11}} values of \eqn{f(S_1)}.}
\item{dev_S0}{The vector of deviance values of the normal mixture for \eqn{f(S_0)}.}
\item{dev_S1}{The vector of deviance values of the normal mixture for \eqn{f(S_1)}.}
\item{code_nlm_0}{An integer indicating why the optimization process to estimate the mixture normal parameters of \eqn{f(S_0)} terminated: 1) relative gradient is close to zero, current iterate is probably solution; 2) successive iterates within tolerance, current iterate is probably solution; 3) last global step failed to locate a point lower than the estimate, the estimate might be an approximate local minimum of the function.}
\item{code_nlm_1}{An integer indicating why the optimization process to estimate the mixture normal parameters of \eqn{f(S_1)} terminated: 1) relative gradient is close to zero, current iterate is probably solution; 2) successive iterates within tolerance, current iterate is probably solution; 3) last global step failed to locate a point lower than the estimate, the estimate might be an approximate local minimum of the function.}
\item{mean.S0}{The mean of \eqn{S_0}.}
\item{var.S0}{The variance of \eqn{S_0}.}
\item{mean.S1}{The mean of \eqn{S_1}.}
\item{var.S1}{The variance of \eqn{S_1}.}
}
\references{
Alonso Abad, A., Ong, F., Stijven, F., Van der Elst, W., Molenberghs, G., Van Keilegom, I., Verbeke, G., & Callegaro, A. (2023). An information-theoretic approach for the assessment of a continuous outcome as a surrogate for a binary true endpoint based on causal inference: Application to vaccine evaluation.
}
\author{
Wim Van der Elst, Fenny Ong, Ariel Alonso, and Geert Molenberghs
}
\seealso{
\code{\link{ICA.ContCont}}, \code{\link{MICA.ContCont}}, \code{\link{ICA.BinBin}}
}
\examples{
\dontrun{# Time consuming code part
data(Schizo)
Fit <- ICA.BinCont(Dataset = Schizo, Surr = BPRS, True = PANSS_Bin, 
Theta.S_0=c(-10,-5,5,10,10,10,10,10), Theta.S_1=c(-10,-5,5,10,10,10,10,10), 
Treat=Treat, M=50, Seed=1)

summary(Fit)
plot(Fit)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{Continuous endpoint}
\keyword{Binary endpoint}