\name{SurvSurv}
\alias{SurvSurv}
\title{Assess surrogacy for two survival endpoints based on information theory and a two-stage approach}
\description{
The function \code{SurvSurv} implements the information-theoretic approach to estimate individual-level surrogacy (i.e., \eqn{R^2_{h.ind}}) and the two-stage approach to estimate trial-level surrogacy (\eqn{R^2_{trial}}, \eqn{R^2_{ht}}) when both endpoints are time-to-event variables (Alonso & Molenberghs, 2008). See the \bold{Details} section below.   
}
\usage{
SurvSurv(Dataset, Surr, SurrCens, True, TrueCens, Treat,
Trial.ID, Weighted=TRUE, Alpha=.05, Number.Bootstraps=500,  
Seed=sample(1:1000, size=1))
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
\item{Number.Bootstraps}{The standard error and confidence interval for \eqn{R^2_{h.ind}} is determined based on a bootstrap procedure. \code{Number.Bootstraps} specifies the number of bootstrap samples that are used. Default \eqn{500}.}
\item{Seed}{The seed to be used in the bootstrap procedure. Default \eqn{sample(1:1000, size=1)}.}
}
\details{
\emph{Individual-level surrogacy}

Alonso & Molenbergs (2008) proposed to redefine the surrogate endpoint \eqn{S} as a time-dependent covariate \eqn{S(t)}, taking value \eqn{0} until the surrogate endpoint occurs and \eqn{1} thereafter. Furthermore, these author considered the models

\deqn{\lambda [t \mid x_{ij}, \beta] = K_{ij}(t) \lambda_{0i}(t) exp(\beta x_{ij}),}
\deqn{\lambda [t \mid x_{ij}, s_{ij}, \beta, \phi] = K_{ij}(t) \lambda_{0i}(t) exp(\beta x_{ij} + \phi S_{ij}),}

where \eqn{K_{ij}(t)} is the risk function for patient \eqn{j} in trial \eqn{i}, \eqn{x_{ij}} is a p-dimensional vector of (possibly) time-dependent covariates, \eqn{\beta} is a p-dimensional vector of unknown coefficients, \eqn{\lambda_{0i}(t)} is a trial-specific baseline hazard function, \eqn{S_{ij}} is a time-dependent covariate version of the surrogate endpoint, and \eqn{\phi} its associated effect.

The mutual information between \eqn{S} and \eqn{T} is estimated as \eqn{I(T,S)=\eqn{1}{n}G^2}, where \eqn{n} is the number of patients and \eqn{G^2} is the log likelihood test comparing the previous two models. Individual-level surrogacy can then be estimated as 

\deqn{R^2_{h.ind} = 1 - exp \left(-\frac{1}{n}G^2 \right).}

O'Quigley and Flandre (2006) pointed out that the previous estimator depends upon the censoring mechanism, even when the censoring mechanism is non-informative. For low levels of censoring this may not be an issue of much concern but for high levels it could lead to biased results. To properly cope with the censoring mechanism in time-to-event outcomes, these authors proposed to estimate the mutual information as \eqn{{I}(T,S)=\frac{1}{k}G^2}, where \eqn{k} is the total number of events experienced. Individual-level surrogacy is then estimated as
\deqn{R^2_{h.ind} = 1 - exp \left(-\frac{1}{k}G^2 \right).}


\emph{Trial-level surrogacy}

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

\value{An object of class \code{SurvSurv} with components,
  \item{Results.Stage.1}{The results of stage 1 of the two-stage model fitting approach: a \code{data.frame} that contains the trial-specific log hazard ratio estimates of the treatment effects for the surrogate and the true endpoints.}
  \item{Results.Stage.2}{An object of class \code{lm} (linear model) that contains the parameter estimates of the regression model that is fitted in stage 2 of the analysis.}
  \item{R2.ht}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{ht}}), its standard error and confidence interval.}
  \item{R2.hind}{A \code{data.frame} that contains the individual-level coefficient of determination (\eqn{R^2_{hind}}), its standard error and confidence interval.}
\item{R2.hInd.By.Trial}{A \code{data.frame} that contains individual-level surrogacy estimates \eqn{R^2_{hInd}} (cluster-based estimates) and their confidence interval for each of the trials seperately.}
}
\references{
Alonso, A. A., & Molenberghs, G. (2008). Evaluating time-to-cancer recurrence as a surrogate marker for survival from an information theory perspective. \emph{Statistical Methods in Medical Research, 17,} 497-504. 

Buyse, M., Michiels, S., Squifflet, P., Lucchesi, K. J., Hellstrand, K., Brune, M. L., Castaigne, S., Rowe, J. M. (2011). Leukemia-free survival as a surrogate end point for overall survival in the evaluation of maintenance therapy for patients with acute myeloid leukemia in complete remission. \emph{Haematologica, 96,} 1106-1112.

O'Quigly, J., & Flandre, P. (2006). Quantification of the Prentice criteria for surrogate endpoints. \emph{Biometrics, 62,} 297-300.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\link{plot.SurvSurv}
}
\examples{
# Open Ovarian dataset
data(Ovarian)

# Conduct analysis
Fit <- SurvSurv(Dataset = Ovarian, Surr = Pfs, SurrCens = PfsInd,
True = Surv, TrueCens = SurvInd, Treat = Treat, 
Trial.ID = Center, Number.Bootstraps=50, Alpha=.05, Seed=1)

# Examine results 
plot(Fit)
summary(Fit)
}
\keyword{Multiple-trial setting}
\keyword{Information-theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Cox proportional hazards model}
\keyword{Survival endpoints}
\keyword{Time-to-event endpoints}