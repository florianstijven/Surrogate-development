\name{MufixedContCont.MultS}
\alias{MufixedContCont.MultS}
\title{Fits a multivariate fixed-effects model to assess surrogacy in the meta-analytic multiple-trial setting (Continuous-continuous case with multiple surrogates)}
\description{
The function \code{MufixedContCont.MultS} uses the multivariate fixed-effects approach to estimate trial- and individual-level surrogacy when the data of multiple clinical trials are available and multiple surrogates are considered for a single true endpoint. The user can specify whether a (weighted or unweighted) full or reduced model should be fitted. See the \bold{Details} section below. 
}
\usage{
MufixedContCont.MultS(Dataset, Endpoints=True~Surr.1+Surr.2, 
 Treat="Treat", Trial.ID="Trial.ID", Pat.ID="Pat.ID", 
 Model=c("Full"), Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, 
 Number.Bootstraps=0, Seed=123)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains one or more  surrogate value(s), a true endpoint value, a treatment indicator, a patient ID, and a trial ID.}
  \item{Endpoints}{An equation in the form \code{True~Surr.1+Surr.2} that specifies the true endpoint followed by the surrogate endpoint(s).}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group.}
  \item{Trial.ID}{The name of the variable in \code{Dataset} that contains the trial ID to which the patient belongs.}
  \item{Pat.ID}{The name of the variable in \code{Dataset} that contains the patient's ID.}
  \item{Model}{The type of model that should be fitted, i.e., \code{Model=c("Full")} or \code{Model=c("Reduced")}. For details, see below or Van der Elst \emph{et al.} (2023). 
  Default \code{Model=c("Full")}. }
\item{Weighted}{Logical. If \code{TRUE}, then a weighted regression analysis is conducted at stage 2 of the two-stage approach. If \code{FALSE}, then an unweighted regression analysis is conducted at stage 2 of the two-stage approach. See the \bold{Details} section below. Default \code{TRUE}.}
  \item{Min.Trial.Size}{The minimum number of patients that a trial should contain in order to be included in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded from the analysis. Default \eqn{2}.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}} and \eqn{R^2_{indiv}}. Default \eqn{0.05}.}
\item{Number.Bootstraps}{Lee's (Lee, 1971) approach is done by default to obtain confidence intervals around \eqn{R^2_{trial}} and \eqn{R^2_{indiv}}. Alternatively, a non-parametric bootstrap can be done. By default, \code{Number.Bootstraps=0} and thus no bootstrap is conducted. If a bootstrap is desired, specify the number of bootstrap samples used this argument. For example, \code{Number.Bootstraps=100} conducts a bootstrap with 100 bootstrap samples.}
  \item{Seed}{The seed that is used in the bootstrap. Default \code{Seed=123}.}
}
\details{
When the full multivariate mixed-effects model is fitted to assess surrogacy in the meta-analytic framework (for details, see Van der Elst \emph{et al.}, 2023), computational issues often occur. In that situation, the use of simplified model-fitting strategies may be warranted (for details, see Burzykowski et al., 2005; Tibaldi et al., 2003).   

The function \code{MufixedContCont.MultS} implements one such strategy, i.e., it uses a two-stage multivariate fixed-effects modelling approach to assess surrogacy. 
In the first stage of the analysis, a multivariate linear regression model is fitted. When a full model is requested (by using the argument \code{Model=c("Full")} in the function call), the following model is fitted:

\deqn{S1_{ij}=\mu_{S1i}+\alpha_{S1i}Z_{ij}+\varepsilon_{S1ij},}
\deqn{S2_{ij}=\mu_{S2i}+\alpha_{S2i}Z_{ij}+\varepsilon_{S2ij},}
\deqn{SK_{ij}=\mu_{SKi}+\alpha_{SKi}Z_{ij}+\varepsilon_{SKij},}
\deqn{T_{ij}=\mu_{Ti}+\beta_{Ti}Z_{ij}+\varepsilon_{Tij},} 

where \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}, \eqn{\mu_{S1i}}, \eqn{\mu_{S2i}}, ..., \eqn{\mu_{SKi}} and \eqn{\mu_{Ti}} are the fixed trial-specific intercepts for \eqn{S1}, \eqn{S2}, ... \eqn{SK} and \eqn{T}, and \eqn{\alpha_{S1i}}, \eqn{\alpha_{S2i}}, ..., \eqn{\alpha_{SKi}} and \eqn{\beta_{Ti}} are the trial-specific treatment effects on the surrogates and the true endpoint, respectively. When a reduced model is requested (by using the argument \code{Model=c("Reduced")} in the function call), the following model is fitted:

\deqn{S1_{ij}=\mu_{S1}+\alpha_{S1i}Z_{ij}+\varepsilon_{S1ij},}
\deqn{S2_{ij}=\mu_{S2}+\alpha_{S2i}Z_{ij}+\varepsilon_{S2ij},}
\deqn{SK_{ij}=\mu_{SK}+\alpha_{SKi}Z_{ij}+\varepsilon_{SKij},}
\deqn{T_{ij}=\mu_{Ti}+\beta_{Ti}Z_{ij}+\varepsilon_{Tij},} 

where \eqn{\mu_{S1}}, \eqn{\mu_{S2}}, ..., \eqn{\mu_{SK}} and \eqn{\mu_{T}} are the common intercepts for the surrogates and the true endpoint (i.e., it is assumed that the intercepts for the surrogates and the true endpoints are identical in all trials). The other parameters are the same as defined above. 

In the above models, the error terms \eqn{\varepsilon_{S1ij}}, \eqn{\varepsilon_{S2ij}}, ..., \eqn{\varepsilon_{SKij}} and \eqn{\varepsilon_{Tij}} are assumed to be mean-zero normally distributed with variance-covariance matrix \eqn{\bold{\Sigma}}. 

Next, the second stage of the analysis is conducted. When a full model is requested by the user (by using the argument  \code{Model=c("Full")} in the function call), the following model is fitted:

\deqn{\widehat{\beta}_{Ti}=\lambda_{0}+\lambda_{1}\widehat{\mu}_{S1i}+
\lambda_{2}\widehat{\alpha}_{S1i}+\lambda_{3}\widehat{\mu}_{S2i}+\lambda_{4}\widehat{\alpha}_{S2i}+...+
\lambda_{2K-1}\widehat{\mu}_{SKi}+\lambda_{2K}\widehat{\alpha}_{SKi}+\varepsilon_{i},}

where the parameter estimates are based on the full model that was fitted in stage 1. 

When a reduced model is requested by the user (by using the argument \code{Model=c("Reduced")}), the \eqn{\lambda_{1} \widehat{\mu}_{S1i}}, \eqn{\lambda_{3} \widehat{\mu}_{S2i}}, ... and \eqn{\lambda_{2K} \widehat{\mu}_{SKi}} components are dropped from the above expression. 

When the argument \code{Weighted=FALSE} is used in the function call, the model that is fitted in stage 2 is an unweighted linear regression model. When a weighted model is requested (using the argument \code{Weighted=TRUE} in the function call), the information that is obtained in stage 1 is weighted according to the number of patients in a trial. 

The classical coefficient of determination of the fitted stage 2 model provides an estimate of \eqn{R^2_{trial}}.
}
\value{An object of class \code{MufixedContCont.MultS} with components,
  \item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of patients who have a missing value for the surrogate and/or the true endpoint are excluded. In addition, the data of trials (i) in which only one type of the treatment was administered, and (ii) in which either the surrogate or the true endpoint was a constant are excluded. In addition, the user can specify the minimum number of patients that a trial should contain in order to include the trial in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded. \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
  \item{Obs.Per.Trial}{A \code{data.frame} that contains the total number of patients per trial and the number of patients who were administered the control treatment and the experimental treatment in each of the trials (in \code{Data.Analyze}).}
  \item{Results.Stage.1}{The results of stage 1 of the two-stage model fitting approach: a \code{data.frame} that contains the trial-specific intercepts and treatment effects for the surrogate(s) and the true endpoints (when a full model is requested), or the trial-specific treatment effects for the surrogates and the true endpoints (when a reduced model is requested).}
  \item{Residuals.Stage.1}{A \code{data.frame} that contains the residuals for the surrogate and true endpoints that are obtained in stage 1 of the analysis (\eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}}).}
  \item{Results.Stage.2}{An object of class \code{lm} (linear model) that contains the parameter estimates of the regression model that is fitted in stage 2 of the analysis.}
  \item{Trial.R2.Lee}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval based on the approach of Lee (1971).}
    \item{Trial.R2.Boot}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval based on the non-parametric bootstrap.}
      \item{Trial.R2.Adj.Lee}{A \code{data.frame} that contains the adjusted trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval based on the approach of Lee (1971).}
    \item{Trial.R2.Adj.Boot}{A \code{data.frame} that contains the adjusted trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval based on the non-parametric bootstrap.}
  \item{Indiv.R2.Lee}{A \code{data.frame} that contains the individual-level coefficient of determination (\eqn{R^2_{indiv}}), its standard error and confidence interval based on the approach of Lee (1971).}
    \item{Indiv.R2.Boot}{A \code{data.frame} that contains the individual-level coefficient of determination (\eqn{R^2_{indiv}}), its standard error and confidence interval based on the non-parametric bootstrap.}
  
    \item{Fitted.Model.Stage.1}{The fitted Stage 1 model.}
    \item{Model.R2.Indiv}{A linear model that regresses the residuals of T on the residuals of the different surrogates.}
  \item{D.Equiv}{The variance-covariance matrix of the trial-specific intercept and treatment effects for the surrogates and true endpoints (when a full model is fitted, i.e., when \code{Model=c("Full")} is used in the function call), or the variance-covariance matrix of the trial-specific treatment effects for the surrogates and true endpoints (when a reduced model is fitted, i.e., when \code{Model=c("Reduced")} is used in the function call).  The variance-covariance matrix \code{D.Equiv} is equivalent to the \eqn{\bold{D}} matrix that would be obtained when a (full or reduced) mixed-effect approach is used; see function \code{\link{MumixedContCont.MultS}}).}
}
\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.

Lee, Y. S. (1971). Tables of the upper percentage points of the multiple correlation. \emph{Biometrika, 59,} 175-189.

Tibaldi, F., Abrahantes, J. C., Molenberghs, G., Renard, D., Burzykowski, T., Buyse, M., Parmar, M., et al., (2003). Simplified hierarchical linear models for the evaluation of surrogate endpoints. \emph{Journal of Statistical Computation and Simulation, 73,} 643-658.

Van der Elst \emph{et al.} (2024). Multivariate surrogate endpoints for normally distributed continuous endpoints in the meta-analytic setting.
}
\author{
Wim Van der Elst
}
\seealso{\code{\link{MumixedContCont.MultS}}}
\examples{
\dontrun{ # time consuming code part
data(PANSS)

# Do a surrogacy analysis with T=Total PANSS score, S1=Negative symptoms
# and S2=Positive symptoms
# Fit a full multivariate fixed-effects model with weighting according to the  
# number of patients in stage 2 of the two stage approach to assess surrogacy:
Fit.Neg.Pos <- MufixedContCont.MultS(Dataset = PANSS, 
  Endpoints = Total ~ Neg+Pos, Model = "Full", 
  Treat = "Treat", Trial.ID = "Invest", Pat.ID = "Pat.ID")
  
# Obtain a summary of the results
summary(Fit.Neg.Pos)
}
}

\keyword{Multiple-trial setting}
\keyword{Meta-analytic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Fixed-effect models}