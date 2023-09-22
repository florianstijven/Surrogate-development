\name{ICA.BinCont.BS}
\alias{ICA.BinCont.BS}
\title{Assess surrogacy in the causal-inference single-trial setting in the binary-continuous case with an additional bootstrap procedure before the assessment}
\description{
The function \code{ICA.BinCont.BS} quantifies surrogacy in the single-trial setting within the causal-inference framework (individual causal association) when the surrogate endpoint is continuous (normally distributed) and the true endpoint is a binary outcome. This function also allows for an additional bootstrap procedure before the assessment to take the imprecision due to finite sample size into account. For details, see Alonso Abad \emph{et al.} (2023). }
\usage{
ICA.BinCont.BS(Dataset, Surr, True, Treat, 
  BS=TRUE,
  nb=300,
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
  \item{BS}{Logical. If \code{BS=TRUE}, the additional bootstrap procedure is performed before the sensitivity analysis to account for the the imprecision due to finite sample size. Default \code{BS=TRUE}.}
  \item{nb}{The number of bootstrap. Default \code{nb=300}.}
  \item{G_pi_10}{The lower and upper limits of the uniform distribution from which the probability parameter \eqn{\pi_{10}} is sampled. Default \code{c(0,1)}. Even though the default is \code{c(0,1)}, due to the restriction that all \eqn{\pi_{ij}} should be between \code{(0,1)}, the value of \eqn{\pi_{10}} will always be between \eqn{(0,min(\pi_{1.}, \pi_{.0}))}. When \code{Monotonicity=TRUE} the values of these limits are set as \code{c(0,0)}.}
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
\item{nboots}{The identification number of bootstrap samples being analyzed in the sensitivity analysis.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{pi_00}{The vector of \eqn{\pi_{00}^T} values.}
\item{pi_01}{The vector of \eqn{\pi_{01}^T} values.}
\item{pi_10}{The vector of \eqn{\pi_{10}^T} values.}
\item{pi_11}{The vector of \eqn{\pi_{11}^T} values.}
\item{G_rho_01_00}{The vector of the \eqn{\rho_{01}^{00}} values.}
\item{G_rho_01_01}{The vector of the \eqn{\rho_{01}^{01}} values.}
\item{G_rho_01_10}{The vector of the \eqn{\rho_{01}^{10}} values.}
\item{G_rho_01_11}{The vector of the \eqn{\rho_{01}^{11}} values.}
\item{mu_0_00}{The vector of mean \eqn{\mu_{0}^{00}} values of \eqn{f(S_0)}.}
\item{mu_0_01}{The vector of mean \eqn{\mu_{0}^{01}} values of \eqn{f(S_0)}.}
\item{mu_0_10}{The vector of mean \eqn{\mu_{0}^{10}} values of \eqn{f(S_0)}.}
\item{mu_0_11}{The vector of mean \eqn{\mu_{0}^{11}} values of \eqn{f(S_0)}.}
\item{mu_1_00}{The vector of mean \eqn{\mu_{1}^{00}} values of \eqn{f(S_1)}.}
\item{mu_1_01}{The vector of mean \eqn{\mu_{1}^{01}} values of \eqn{f(S_1)}.}
\item{mu_1_10}{The vector of mean \eqn{\mu_{1}^{10}} values of \eqn{f(S_1)}.}
\item{mu_1_11}{The vector of mean \eqn{\mu_{1}^{11}} values of \eqn{f(S_1)}.}
\item{sigma_00}{The vector of variance \eqn{\sigma_{00}} values of \eqn{f(S_0)}.}
\item{sigma_11}{The vector of variance \eqn{\sigma_{11}} values of \eqn{f(S_1)}.}
}
\references{
Alonso Abad, A., Ong, F., Stijven, F., Van der Elst, W., Molenberghs, G., Van Keilegom, I., Verbeke, G., & Callegaro, A. (2023). An information-theoretic approach for the assessment of a continuous outcome as a surrogate for a binary true endpoint based on causal inference: Application to vaccine evaluation.
}
\author{
Wim Van der Elst, Fenny Ong, Ariel Alonso, and Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinCont}}
}
\examples{
\dontrun{# Time consuming code part
data(Schizo)
Fit <- ICA.BinCont.BS(Dataset = Schizo, Surr = BPRS, True = PANSS_Bin, nb = 10, 
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