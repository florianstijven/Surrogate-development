\name{UnifixedContCont}
\alias{UnifixedContCont}
\title{Fits univariate fixed-effect models to assess surrogacy in the meta-analytic multiple-trial setting (continuous-continuous case)
}
\description{
The function \code{UnifixedContCont} uses the univariate fixed-effects approach to estimate trial- and individual-level surrogacy when the data of multiple clinical trials are available. The user can specify whether a (weighted or unweighted) full, semi-reduced, or reduced model should be fitted. See the \bold{Details} section below.   
}
\usage{
UnifixedContCont(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, Model=c("Full"), 
Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, Number.Bootstraps=500)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, a treatment indicator, a patient ID, and a trial ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate endpoint values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should either be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group, or as \eqn{1} for the experimental group and \eqn{0} for the control group. The \eqn{-1/1} coding is recommended. }
  \item{Trial.ID}{The name of the variable in \code{Dataset} that contains the trial ID to which the patient belongs.}
  \item{Pat.ID}{The name of the variable in \code{Dataset} that contains the patient's ID.}
  \item{Model}{The type of model that should be fitted, i.e., \code{Model=c("Full")}, \code{Model=c("Reduced")}, or \code{Model=c("SemiReduced")}. See the \bold{Details} section below. 
  Default \code{Model=c("Full")}. }
\item{Weighted}{Logical. If \code{TRUE}, then a weighted regression analysis is conducted at stage 2 of the two-stage approach. If \code{FALSE}, then an unweighted regression analysis is conducted at stage 2 of the two-stage approach. See the \bold{Details} section below. Default \code{TRUE}.}
\item{Min.Trial.Size}{The minimum number of patients that a trial should contain to be included in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded from the analysis. Default \eqn{2}.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}}, \eqn{R_{trial}}, \eqn{R^2_{indiv}}, and \eqn{R_{indiv}}. Default \eqn{0.05}.}
\item{Number.Bootstraps}{The standard errors and confidence intervals for \eqn{R^2_{indiv}} and \eqn{R_{indiv}} are determined as based on a bootstrap procedure. \code{Number.Bootstraps} specifies the number of bootstrap samples that are used. Default \eqn{500}.}
}
\details{
When the full bivariate mixed-effects model is fitted to assess surrogacy in the meta-analytic framework (for details, Buyse & Molenberghs, 2000), computational issues often occur (e.g., the Newton Raphson algorithm may not converge). In that situation, the use of simplified model-fitting strategies may be warranted (for details, see see Burzykowski et al., 2005; Tibaldi et al., 2003). 

The function \code{UnifixedContCont} implements one such strategy, i.e., it uses a two-stage univariate fixed-effects modelling approach to assess surrogacy. In the first stage of the analysis, two univariate linear regression models are fitted to the data of each of the \eqn{i} trials. When a full or semi-reduced model is requested (by using the argument \code{Model=c("Full")} or \code{Model=c("SemiReduced")} in the function call), the following univariate models are fitted:

\deqn{S_{ij}=\mu_{Si}+\alpha_{i}Z_{ij}+\varepsilon_{Sij},}
\deqn{T_{ij}=\mu_{Ti}+\beta_{i}Z_{ij}+\varepsilon_{Tij},} 

where \eqn{i} and \eqn{j} are the trial and subject indicators, \eqn{S_{ij}} and \eqn{T_{ij}} are the surrogate and true endpoint values of subject \eqn{j} in trial \eqn{i}, \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}, \eqn{\mu_{Si}} and \eqn{\mu_{Ti}} are the fixed trial-specific intercepts for S and T, and \eqn{\alpha_{i}} and \eqn{\beta_{i}} are the fixed trial-specific treatment effects on S and T, respectively. The error terms \eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}} are assumed to be independent. 

When a reduced model is requested by the user (by using the argument \code{Model=c("Reduced")} in the function call), the following univariate models are fitted:

\deqn{S_{ij}=\mu_{S}+\alpha_{i}Z_{ij}+\varepsilon_{Sij},}
\deqn{T_{ij}=\mu_{T}+\beta_{i}Z_{ij}+\varepsilon_{Tij},}

where \eqn{\mu_{S}} and \eqn{\mu_{T}} are the common intercepts for S and T (i.e., it is assumed that the intercepts for the surrogate and the true endpoints are identical in each of the trials). The other parameters are the same as defined above, and \eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}} are again assumed to be independent. 

An estimate of \eqn{R^2_{indiv}} is provided by \eqn{r(\varepsilon_{Sij}, \varepsilon_{Tij})^2}. 

Next, the second stage of the analysis is conducted. When a full model is requested (by using the argument \code{Model=c("Full")} in the function call), the following model is fitted:
\deqn{\widehat{\beta}_{i}=\lambda_{0}+\lambda_{1}\widehat{\mu_{Si}}+\lambda_{2}\widehat{\alpha}_{i}+\varepsilon_{i},}
where the parameter estimates for \eqn{\beta_i}, \eqn{\mu_{Si}}, and \eqn{\alpha_i} are based on the full models that were fitted in stage 1. 

When a semi-reduced or reduced model is requested (by using the argument \code{Model=c("SemiReduced")} or \code{Model=c("Reduced")} in the function call), the following model is fitted:
\deqn{\widehat{\beta}_{i}=\lambda_{0}+\lambda_{1}\widehat{\alpha}_{i}+\varepsilon_{i}.}

where the parameter estimates for \eqn{\beta_i} and \eqn{\alpha_i} are based on the semi-reduced or reduced models that were fitted in stage 1.

When the argument \code{Weighted=FALSE} is used in the function call, the model that is fitted in stage 2 is an unweighted linear regression model. When a weighted model is requested (using the argument \code{Weighted=TRUE} in the function call), the information that is obtained in stage 1 is weighted according to the number of patients in a trial. 

The classical coefficient of determination of the fitted stage 2 model provides an estimate of \eqn{R^2_{trial}}. 
}
\value{An object of class \code{UnifixedContCont} with components,
  \item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of patients who have a missing value for the surrogate and/or the true endpoint are excluded. In addition, the data of trials (i) in which only one type of the treatment was administered, and (ii) in which either the surrogate or the true endpoint was a constant (i.e., all patients within a trial had the same surrogate and/or true endpoint value) are excluded. In addition, the user can specify the minimum number of patients that a trial should contain in order to include the trial in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded. \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
  \item{Obs.Per.Trial}{A \code{data.frame} that contains the total number of patients per trial and the number of patients who were administered the control treatment and the experimental treatment in each of the trials (in \code{Data.Analyze}).}
  \item{Results.Stage.1}{The results of stage 1 of the two-stage model fitting approach: a \code{data.frame} that contains the trial-specific intercepts and treatment effects for the surrogate and the true endpoints (when a full or semi-reduced model is requested), or the trial-specific treatment effects for the surrogate and the true endpoints (when a reduced model is requested).}
  \item{Residuals.Stage.1}{A \code{data.frame} that contains the residuals for the surrogate and true endpoints that are obtained in stage 1 of the analysis (\eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}}).}
  \item{Results.Stage.2}{An object of class \code{lm} (linear model) that contains the parameter estimates of the regression model that is fitted in stage 2 of the analysis.}
  \item{Trial.R2}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval.}
  \item{Indiv.R2}{A \code{data.frame} that contains the individual-level coefficient of determination (\eqn{R^2_{indiv}}), its standard error and confidence interval.}
  \item{Trial.R}{A \code{data.frame} that contains the trial-level correlation coefficient (\eqn{R_{trial}}), its standard error and confidence interval.}
  \item{Indiv.R}{A \code{data.frame} that contains the individual-level correlation coefficient (\eqn{R_{indiv}}), its standard error and confidence interval.}
  \item{Cor.Endpoints}{A \code{data.frame} that contains the correlations between the surrogate and the true endpoint in the control treatment group (i.e., \eqn{\rho_{T0S0}}) and in the experimental treatment group (i.e., \eqn{\rho_{T1S1}}), their standard errors and their confidence intervals.}
\item{D.Equiv}{The variance-covariance matrix of the trial-specific intercept and treatment effects for the surrogate and true endpoints (when a full or semi-reduced model is fitted, i.e., when \code{Model=c("Full")} or \code{Model=c("SemiReduced")} is used in the function call), or the variance-covariance matrix of the trial-specific treatment effects for the surrogate and true endpoints (when a reduced model is fitted, i.e., when \code{Model=c("Reduced")} is used in the function call).  The variance-covariance matrix \code{D.Equiv} is equivalent to the \eqn{\bold{D}} matrix that would be obtained when a (full or reduced) bivariate mixed-effect approach is used; see function \code{\link{BimixedContCont}}).}
}
\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.

Tibaldi, F., Abrahantes, J. C., Molenberghs, G., Renard, D., Burzykowski, T., Buyse, M., Parmar, M., et al., (2003). Simplified hierarchical linear models for the evaluation of surrogate endpoints. \emph{Journal of Statistical Computation and Simulation, 73,} 643-658.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{UnimixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{BimixedContCont}}, \code{\link{plot Meta-Analytic}}
}
\examples{
\dontrun{
# Example 1, based on the ARMD data
data(ARMD)

# Fit a full univariate fixed-effects model with weighting according to the  
# number of patients in stage 2 of the two stage approach to assess surrogacy:
Sur <- UnifixedContCont(Dataset=ARMD, Surr=Diff24, True=Diff52, Treat=Treat, Trial.ID=Center, 
Pat.ID=Id, Model="Full", Weighted=TRUE)

# Obtain a summary and plot of the results
summary(Sur)
plot(Sur)

# Example 2
# Conduct an analysis based on a simulated dataset with 2000 patients, 100 trials, 
# and Rindiv=Rtrial=.8
# Simulate the data:
Sim.Data.MTS(N.Total=2000, N.Trial=100, R.Trial.Target=.8, R.Indiv.Target=.8,
Seed=123, Model="Reduced")

# Fit a reduced univariate fixed-effects model without weighting to assess
# surrogacy:
Sur2 <- UnifixedContCont(Dataset=Data.Observed.MTS, Surr=Surr, True=True, Treat=Treat, 
Trial.ID=Trial.ID, Pat.ID=Pat.ID, Model="Reduced", Weighted=FALSE)

# Show a summary and plots of results:
summary(Sur2)
plot(Sur2, Weighted=FALSE)}
}

\keyword{Multiple-trial setting}
\keyword{Meta-analytic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Fixed-effect models}