\name{TwoStageSurvSurv}
\alias{TwoStageSurvSurv}
\title{Assess trial-level surrogacy for two survival endpoints using a two-stage approach}
\description{
The function \code{TwoStageSurvSurv} uses a two-stage approach to estimate \eqn{R^2_{trial}}. In stage 1, trial-specific Cox proportional hazard models are fitted and in stage 2 the trial-specific estimated treatment effects on \eqn{T} are regressed on the trial-specific estimated treatment effects on \eqn{S} (measured on the log hazard ratio scale). The user can specify whether a weighted or unweighted model should be fitted at stage 2. See the \bold{Details} section below.   
}
\usage{
TwoStageSurvSurv(Dataset, Surr, SurrCens, True, TrueCens, Treat,
Trial.ID, Weighted=TRUE, Alpha=.05)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value and censoring indicator, a true endpoint value and censoring indicator, a treatment indicator, and a trial ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate endpoint values.}
  \item{SurrCens}{The name of the variable in \code{Dataset} that contains the censoring indicator for the surrogate endpoint values (1 = event, 0 = censored).}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{TrueCens}{The name of the variable in \code{Dataset} that contains the censoring indicator for the true endpoint values (1 = event, 0 = censored).}
\item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. }
  \item{Trial.ID}{The name of the variable in \code{Dataset} that contains the trial ID to which the patient belongs.}
\item{Weighted}{Logical. If \code{TRUE}, then a weighted regression analysis is conducted at stage 2 of the two-stage approach. If \code{FALSE}, then an unweighted regression analysis is conducted at stage 2 of the two-stage approach. See the \bold{Details} section below. Default \code{TRUE}.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}} and \eqn{R_{trial}}. Default \eqn{0.05}.}
}
\details{
A two-stage approach is used to estimate trial-level surrogacy, following a procedure proposed by Buyse et al. (2011). In stage 1, the following trial-specific Cox proportional hazard models are fitted:

\deqn{S_{ij}(t)=S_{i0}(t) exp(\alpha_{i}Z_{ij}),}
\deqn{T_{ij}(t)=T_{i0}(t) exp(\beta_{i}Z_{ij}),} 

where \eqn{S_{i0}(t)} and \eqn{T_{i0}(t)} are the trial-specific baseline hazard functions, \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}, \eqn{\mu_{Si}}, and \eqn{\alpha_{i}} and \eqn{\beta_{i}} are the trial-specific treatment effects on S and T, respectively. 

Next, the second stage of the analysis is conducted:

\deqn{\widehat{\beta_{i}}=\lambda_{0}+\lambda_{1}\widehat{\alpha_{i}}+\varepsilon_{i},}

where the parameter estimates for \eqn{\beta_i}, \eqn{\mu_{Si}}, and \eqn{\alpha_i} are based on the full model that was fitted in stage 1. 

When the argument \code{Weighted=FALSE} is used in the function call, the model that is fitted in stage 2 is an unweighted linear regression model. When a weighted model is requested (using the argument \code{Weighted=TRUE} in the function call), the information that is obtained in stage 1 is weighted according to the number of patients in a trial. 

The classical coefficient of determination of the fitted stage 2 model provides an estimate of \eqn{R^2_{trial}}.
}

\value{An object of class \code{TwoStageSurvSurv} with components,
  \item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of trials that do not have at least three patients per treatment arm are excluded due to estimation constraints (Burzykowski et al., 2001). \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
  \item{Results.Stage.1}{The results of stage 1 of the two-stage model fitting approach: a \code{data.frame} that contains the trial-specific log hazard ratio estimates of the treatment effects for the surrogate and the true endpoints.}
  \item{Results.Stage.2}{An object of class \code{lm} (linear model) that contains the parameter estimates of the regression model that is fitted in stage 2 of the analysis.}
  \item{Trial.R2}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval.}
  \item{Trial.R}{A \code{data.frame} that contains the trial-level correlation coefficient (\eqn{R_{trial}}), its standard error and confidence interval.}
}
\references{
Burzykowski, T., Molenberghs, G., Buyse, M., Geys, H., & Renard, D. (2001). Validation of surrogate endpoints in multiple randomized clinical trials with failure-time endpoints. \emph{Applied Statistics, 50,} 405-422. 

Buyse, M., Michiels, S., Squifflet, P., Lucchesi, K. J., Hellstrand, K., Brune, M. L., Castaigne, S., Rowe, J. M. (2011). Leukemia-free survival as a surrogate end point for overall survival in the evaluation of maintenance therapy for patients with acute myeloid leukemia in complete remission. \emph{Haematologica, 96,} 1106-1112.

}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\link{plot.TwoStageSurvSurv}
}
\examples{
# Open Ovarian dataset
data(Ovarian)

# Conduct analysis
Results <- TwoStageSurvSurv(Dataset = Ovarian, Surr = Pfs, SurrCens = PfsInd, 
True = Surv, TrueCens = SurvInd, Treat = Treat, Trial.ID = Center)

# Examine results of analysis
summary(Results)
plot(Results)
}

\keyword{Multiple-trial setting}
\keyword{Information-theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Cox proportional hazards model}
\keyword{Survival endpoints}
\keyword{Time-to-event endpoints}