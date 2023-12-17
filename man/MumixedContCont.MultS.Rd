\name{MumixedContCont.MultS}
\alias{MumixedContCont.MultS}
\title{Fits a multivariate mixed-effects model to assess surrogacy in the meta-analytic multiple-trial setting (Continuous-continuous case with multiple surrogates)}
\description{
The function \code{MumixedContCont.MultS} uses the multivariate mixed-effects approach to estimate trial- and individual-level surrogacy when the data of multiple clinical trials are available and multiple surrogates are considered for a single true endpoint. See the \bold{Details} section below. 
}
\usage{
MumixedContCont.MultS(Dataset, Endpoints=True~Surr.1+Surr.2, 
Treat="Treat", Trial.ID="Trial.ID", Pat.ID="Pat.ID", 
Model=c("Full"), Min.Trial.Size=2, Alpha=.05, Opt="nlminb")
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains one or more  surrogate value(s), a true endpoint value, a treatment indicator, a patient ID, and a trial ID.}
  \item{Endpoints}{An equation in the form \code{True~Surr.1+Surr.2} that specifies the true endpoint followed by the surrogate endpoint(s).}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group.}
  \item{Trial.ID}{The name of the variable in \code{Dataset} that contains the trial ID to which the patient belongs.}
  \item{Pat.ID}{The name of the variable in \code{Dataset} that contains the patient's ID.}
  \item{Model}{The type of model that should be fitted, i.e., \code{Model=c("Full")} or \code{Model=c("Reduced")}. For details, see below or Van der Elst \emph{et al.} (2023). 
  Default \code{Model=c("Full")}. }
  \item{Min.Trial.Size}{The minimum number of patients that a trial should contain in order to be included in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded from the analysis. Default \eqn{2}.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}} and \eqn{R^2_{indiv}} (based on the approach of Lee, 1971). Default \eqn{0.05}.}
  \item{Opt}{The optimizer to be used by the \code{lme} function (the fits the mixed-effects model), with options \code{nlminb} or \code{optim}. For  details, see \code{?lmeControl}. Default \code{Opt="nlminb"}.}

}
\details{
When a full model is requested (by using the argument \code{Model=c("Full")} in the function call), the following mixed-effects model is fitted:

\deqn{S1_{ij}=\mu_{S1}+m_{S1i}(\alpha_{S1}+a_{S1i})Z_{ij}+\varepsilon_{S1ij},}
\deqn{S2_{ij}=\mu_{S2}+m_{S2i}(\alpha_{S2}+a_{S2i})Z_{ij}+\varepsilon_{S2ij},}
\deqn{SK_{ij}=\mu_{SK}+m_{SKi}(\alpha_{SK}+a_{SKi})Z_{ij}+\varepsilon_{SKij},}
\deqn{T_{ij}=\mu_{T}+m_{Ti}(\beta_{T}+b_{Ti})Z_{ij}+\varepsilon_{Tij},}

where \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}, 
\eqn{\mu_{S1}}, \eqn{\mu_{S2}}, ... \eqn{\mu_{SK}} and \eqn{\mu_{T}} are the fixed intercepts for \eqn{S1}, \eqn{S2}, ... \eqn{SK} and \eqn{T}, \eqn{m_{S1i}}, \eqn{m_{S2i}}, ... \eqn{m_{SKi}}, and \eqn{m_{Ti}} are the corresponding random intercepts, \eqn{\alpha_{S1}}, \eqn{\alpha_{S2}}, ..., \eqn{\alpha_{SK}} and \eqn{\beta_T} are the fixed treatment effects for \eqn{S1}, \eqn{S2}, ... \eqn{SK} and \eqn{T}, and \eqn{a_{S1i}}, \eqn{a_{S2i}}, ... \eqn{a_{SKi}} and \eqn{b_{Ti}} are the corresponding random treatment effects. The vector of the random effects \eqn{\left(m_{S1i},\:m_{S2i}, \: ... , \: m_{SKi},\: m_{Ti},\: a_{S1i},\: a_{S2i},\: ... , \: a_{SKi},\: b_{Ti}\right)} is assumed to be mean-zero normally distributed with unstructured variance-covariance matrix \eqn{\mathbf{D}}. Similarly, the residuals \eqn{\varepsilon_{S1ij}}, \eqn{\varepsilon_{S2ij}}, ... \eqn{\varepsilon_{SKij}},  \eqn{\varepsilon_{Tij}} are assumed to be mean-zero normally distributed with unstructured variance-covariance matrix \eqn{\mathbf{\Sigma}}.

When a reduced model is requested (by using the argument \code{Model=c("Reduced")} in the function call), the trial-specific intercepts for the surrogate endpoints and the true endpoint in the above model are replaced by common intercepts.

For the full model, \eqn{R^2_{trial}} and \eqn{R^2_{indiv}} are estimated based on \eqn{\mathbf{D}} and \eqn{\mathbf{\Sigma}}, respectively:
 \deqn{R_{trial}^{2}=R^2_{b_{Ti}|m_{S1i},\: m_{S2i},\: ..., \:m_{SKi}, \:  a_{S1i},\: a_{S2i}, \: ... \: a_{SKi}}=
\dfrac{\boldsymbol{D}_{ST}^T \: \boldsymbol{D}^{-1}_{SS} \: \boldsymbol{D}_{ST}}{\boldsymbol{D}_{TT}},}
\deqn{R_{indiv}^{2}=R_{\varepsilon_{Tij}|\varepsilon_{S1ij}, \: \varepsilon_{S2ij}, \: ..., \: \varepsilon_{SKij}}^{2}=
\dfrac{\boldsymbol{\Sigma}_{ST}^T \: \boldsymbol{\Sigma}^{-1}_{SS} \: \boldsymbol{\Sigma}_{ST}}{\boldsymbol{\Sigma}_{TT}}.}

For the reduced model, the reduced \eqn{\mathbf{D}} and \eqn{\mathbf{\Sigma}} are used. 
} 
\value{An object of class \code{MumixedContCont.MultS} with components,
  \item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of patients who have a missing value for the surrogate and/or the true endpoint are excluded. In addition, the data of trials (i) in which only one type of the treatment was administered, and (ii) in which either the surrogate or the true endpoint was a constant are excluded. In addition, the user can specify the minimum number of patients that a trial should contain in order to include the trial in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded. \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
  \item{Obs.Per.Trial}{A \code{data.frame} that contains the total number of patients per trial and the number of patients who were administered the control treatment and the experimental treatment in each of the trials (in \code{Data.Analyze}).}
\item{Fixed.Effects}{A \code{data.frame} that contains the fixed intercepts and treatment effects for the true and the surrogate endpoints.}
\item{Random.Effects}{A \code{data.frame} that contains the random intercepts and treatment effects for the true and the surrogate endpoints.}
  \item{Trial.R2.Lee}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval based on the approach of Lee (1971).}
  \item{Indiv.R2.Lee}{A \code{data.frame} that contains the individual-level coefficient of determination (\eqn{R^2_{indiv}}), its standard error and confidence interval based on the approach of Lee (1971).}
  \item{D}{The variance-covariance matrix of the trial-specific intercepts and treatment effects for the surrogates and true endpoints (when a full model is fitted, i.e., when \code{Model=c("Full")} is used in the function call), or the variance-covariance matrix of the trial-specific treatment effects for the surrogates and true endpoints (when a reduced model is fitted, i.e., when \code{Model=c("Reduced")} is used in the function call).}
  \item{Cond.Number.D.Matrix}{The condition number of the \eqn{\mathbf{D}} matrix.}  
    \item{Cond.Number.Sigma.Matrix}{The condition number of the \eqn{\mathbf{\Sigma}} matrix.}  
    \item{Fitted.Model}{The fitted mixed-effects model.}
}
\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.

Lee, Y. S. (1971). Tables of the upper percentage points of the multiple correlation. \emph{Biometrika, 59,} 175-189.

Van der Elst \emph{et al.} (2024). Multivariate surrogate endpoints for normally distributed continuous endpoints in the meta-analytic setting.
}
\author{
Wim Van der Elst
}
\seealso{\code{\link{MufixedContCont.MultS}}}
\examples{
\dontrun{ # time consuming code part
data(PANSS)

# Do a surrogacy analysis with T=Total PANSS score, 
# S1=Negative symptoms and S2=Positive symptoms
# Fit a full mixed-effects model:
Fit.Neg.Pos <- MumixedContCont.MultS(Dataset = PANSS, 
  Endpoints = Total ~ Neg+Pos, Model = "Full", 
  Treat = "Treat", Trial.ID = "Invest", Pat.ID = "Pat.ID")
  
# Model does not converge, as often happens with the 
# mixed-effects approach. Instead, fit a full multivariate 
# fixed-effects model with weighting according to the  
# number of patients in stage 2 of the two stage approach to assess surrogacy:
Fit.Neg.Pos <- MufixedContCont.MultS(Dataset = PANSS, 
  Endpoints = Total ~ Neg+Pos, Model = "Full", 
  Treat = "Treat", Trial.ID = "Invest", Pat.ID = "Pat.ID")
  
# Obtain a summary of the results
summary(Fit.Neg.Pos)
# 
}
}

\keyword{Multiple-trial setting}
\keyword{Meta-analytic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{XXX-effect models}