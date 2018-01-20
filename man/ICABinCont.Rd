\name{ICA.BinCont}
\alias{ICA.BinCont}
\title{Assess surrogacy in the causal-inference single-trial setting in the binary-continuous case}
\description{
The function \code{ICA.BinCont} quantifies surrogacy in the single-trial causal-inference framework (individual causal association) when the surrogate endpoint is continuous (normally distributed) and the true endpoint is a binary outcome. For details, see Alonso \emph{et al.} (2016). }
\usage{
ICA.BinCont(Dataset, Surr, True, Treat, Diff.Sigma=FALSE, 
  G_pi_00=seq(0, 1, by=.01), G_rho_01_00=seq(-1, 1, by=.01), 
  G_rho_01_01=seq(-1, 1, by=.01), G_rho_01_10=seq(-1, 1, by=.01), 
  G_rho_01_11=seq(-1, 1, by=.01), M=1000, Seed=123, 
  Plots=TRUE, Save.Plots="No", Test.Fit.Mixture=FALSE, 
  Test.Fit.Mixture.Alpha=0.01, Test.Fit.Details=FALSE,
  Keep.All=FALSE)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, and a treatment indicator.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate endpoint values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group.}
  \item{Diff.Sigma}{Logical. If \code{Diff.Sigma=TRUE}, then the mixtures of normal distributions are allowed to have different variances. If \code{Diff.Sigma=FALSE}, then the mixtures of normal distributions are not allowed to have different variances (selecting the latter option speeds up the required calculations). Default \code{Diff.Sigma=FALSE}.}
  \item{G_pi_00}{A grid of values to be considered for \eqn{\pi_{11}}, i.e., the unidentifiable probability \eqn{P(T_1=0,T_0=0)}. Default \code{seq(0, 1, by=.01)}.}
  \item{G_rho_01_00}{A grid of values to be considered for the association parameter \eqn{\rho_{01}^{00}}. Default \code{seq(-1, 1, by=.01)}.}
  \item{G_rho_01_01}{A grid of values to be considered for the association parameter \eqn{\rho_{01}^{01}}. Default \code{seq(-1, 1, by=.01)}.}
  \item{G_rho_01_10}{A grid of values to be considered for the association parameter \eqn{\rho_{01}^{10}}. Default \code{seq(-1, 1, by=.01)}.}
  \item{G_rho_01_11}{A grid of values to be considered for the association parameter \eqn{\rho_{01}^{11}}. Default \code{seq(-1, 1, by=.01)}.}
  \item{M}{The number of valid ICA values to be sampled. Default \code{M=1000}.}
  \item{Seed}{The seed to be used to generate \eqn{\pi_r}. Default \code{Seed=123}.}
\item{Plots}{Logical. Should histograms of \eqn{S_0} (surrogate endpoint in control group)  and \eqn{S_1} (surrogate endpoint in experimental treatment group) be provided together with density of fitted mixtures? Default \code{Plots=TRUE}.}
\item{Save.Plots}{Should the plots (see previous item) be saved? If \code{Save.Plots="No"}, no plots are saved. If plots have to be saved, replace \code{"No"} by the desired location, e.g., \code{Save.Plots="C:/"}. Default \code{Save.Plots="No"}.} 
\item{Test.Fit.Mixture}{Should the fit of the densities of the mixture distributions with the observed densities of the surrogates in the control and experimental treatment groups be conducted? For details on the method used, see Wilcox (1995, 2014). The code used to conduct the analysis is provided by Wilcox, see http://dornsife.usc.edu/labs/rwilcox/software/. Default \code{Test.Fit.Mixture=FALSE}.}
\item{Test.Fit.Mixture.Alpha}{The \eqn{alpha}-level that is used in comparing the observed densities of \eqn{S[0]} and \eqn{S[1]}, see previous points. Default \code{Test.Fit.Mixture.Alpha=0.01}.}
\item{Test.Fit.Details}{Should the details of the Wilcox-testing procedure be saved? Default {Test.Fit.Details=FALSE}}. 
\item{Keep.All}{When \code{Test.Fit.Mixture} is used, the Wilcox-testing procedure is used to evaluate model fit and only models for which the fit is OK (i.e., all p-values above the specified \eqn{\alpha}-level) are retained. To keep all results (irrespective of whether or not model fit is OK), \code{Keep.All=TRUE} can be used. Default \code{Keep.All=FALSE}.}
  }
\value{An object of class \code{ICA.BinCont} with components,
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{pi_00}{The vector of \eqn{\pi_{00}^T} values.}
\item{pi_01}{The vector of \eqn{\pi_{01}^T} values.}
\item{pi_10}{The vector of \eqn{\pi_{10}^T} values.}
\item{pi_11}{The vector of \eqn{\pi_{11}^T} values.}
\item{G_rho_01_00}{The vector of the \eqn{R_{01}^{00}} values.}
\item{G_rho_01_01}{The vector of the \eqn{R_{01}^{01}} values.}
\item{G_rho_01_10}{The vector of the \eqn{R_{01}^{10}} values.}
\item{G_rho_01_11}{The vector of the \eqn{R_{01}^{11}} values.}
\item{pi_Delta_T_min1}{The vector of the \eqn{\pi_{-1}^{\Delta T}} values.}
\item{pi_Delta_T_0}{The vector of the \eqn{\pi_{0}^{\Delta T}} values.}
\item{pi_Delta_T_1}{The vector of the \eqn{\pi_{1}^{\Delta T}} values.}
\item{pi_0_00}{The vector of \eqn{\pi_{00}} values of \eqn{f(S_0)}.}
\item{pi_0_01}{The vector of \eqn{\pi_{01}} values of \eqn{f(S_0)}.}
\item{pi_0_10}{The vector of \eqn{\pi_{10}} values of \eqn{f(S_0)}.}
\item{pi_0_11}{The vector of \eqn{\pi_{11}} values of \eqn{f(S_0)}.}
\item{mu_0_00}{The vector of \eqn{\mu_{0}^{00}} values of \eqn{f(S_0)}.}
\item{mu_0_01}{The vector of \eqn{\mu_{0}^{01}} values of \eqn{f(S_0)}.}
\item{mu_0_10}{The vector of \eqn{\mu_{0}^{10}} values of \eqn{f(S_0)}.}
\item{mu_0_11}{The vector of \eqn{\mu_{0}^{11}} values of \eqn{f(S_0)}.}
\item{sigma2_00_00}{The vector of squared \eqn{\sigma_{00}^{00}} values of \eqn{f(S_0)}.}
\item{sigma2_00_01}{The vector of squared \eqn{\sigma_{00}^{01}} values of \eqn{f(S_0)}.}
\item{sigma2_00_10}{The vector of squared \eqn{\sigma_{00}^{10}} values of \eqn{f(S_0)}.}
\item{sigma2_00_11}{The vector of squared \eqn{\sigma_{00}^{11}} values of \eqn{f(S_0)}.}
\item{pi_1_00}{The vector of \eqn{\pi_{00}} values of \eqn{f(S_1)}.}
\item{pi_1_01}{The vector of \eqn{\pi_{01}} values of \eqn{f(S_1)}.}
\item{pi_1_10}{The vector of \eqn{\pi_{10}} values of \eqn{f(S_1)}.}
\item{pi_1_11}{The vector of \eqn{\pi_{11}} values of \eqn{f(S_1)}.}
\item{mu_1_00}{The vector of \eqn{\mu_{1}^{00}} values of \eqn{f(S_1)}.}
\item{mu_1_01}{The vector of \eqn{\mu_{1}^{01}} values of \eqn{f(S_1)}.}
\item{mu_1_10}{The vector of \eqn{\mu_{1}^{10}} values of \eqn{f(S_1)}.}
\item{mu_1_11}{The vector of \eqn{\mu_{1}^{11}} values of \eqn{f(S_1)}.}
\item{sigma2_11_00}{The vector of squared \eqn{\sigma_{11}^{00}} values of \eqn{f(S_1)}.}
\item{sigma2_11_01}{The vector of squared \eqn{\sigma_{11}^{01}} values of \eqn{f(S_1)}.}
\item{sigma2_11_10}{The vector of squared \eqn{\sigma_{11}^{10}} values of \eqn{f(S_1)}.}
\item{sigma2_11_11}{The vector of squared \eqn{\sigma_{11}^{11}} values of \eqn{f(S_1)}.}
\item{Fit.Mixture_S_0_OK}{Is the fit of the mixture distribution for \eqn{S[0]} OK (i.e., all \eqn{p}-values) of the test procedure above the specified \eqn{\alpha}?}
\item{Fit.Mixture_S_1_OK}{Is the fit of the mixture distribution for \eqn{S[1]} OK (i.e., all \eqn{p}-values) of the test procedure above the specified \eqn{\alpha}?}
\item{Test.Fit.Details}{Details of the Wilcox-testing procedure. This information is provided when the argument \code{Test.Fit.Details=FALSE} was used in the function call.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2016). Surrogate markers validation: the continuous-binary setting from a causal inference perspective.

Wilcox, R.R. (1995) Comparing Two Independent Groups Via Multiple Quantiles. Journal of the Royal Statistical Society. Series D (The Statistician), 44, 91-99.

Wilcox, R. R., Erceg-Hurn, D. M, Clark, F., & Carlson, M. (2014). Comparing two independent groups via the lower and upper quantiles. Journal of Statistical Computation and Simulation.
}
\author{
Wim Van der Elst & Ariel Alonso
}
\seealso{
\code{\link{ICA.ContCont}}, \code{\link{MICA.ContCont}}, \code{\link{ICA.BinBin}}
}
\examples{
\dontrun{# Time consuming code part
data(Schizo)
Fit <- ICA.BinCont(Dataset = Schizo, Surr = BPRS, True = PANSS_Bin, 
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