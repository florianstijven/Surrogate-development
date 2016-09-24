\name{MixedContContIT}
\alias{MixedContContIT}
\title{Fits (univariate) mixed-effect models to assess surrogacy in the continuous-continuous case based on the Information-Theoretic framework}
\description{
The function \code{MixedContContIT} uses the information-theoretic approach (Alonso & Molenberghs, 2007) to estimate trial- and individual-level surrogacy based on mixed-effect models when both S and T are continuous endpoints. The user can specify whether a (weighted or unweighted) full, semi-reduced, or reduced model should be fitted. See the \bold{Details} section below.   
}
\usage{
MixedContContIT(Dataset, Surr, True, Treat, Trial.ID, Pat.ID,
Model=c("Full"), Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, ...)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, a treatment indicator, a patient ID, and a trial ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate endpoint values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should either be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group, or as \eqn{1} for the experimental group and \eqn{0} for the control group. }
  \item{Trial.ID}{The name of the variable in \code{Dataset} that contains the trial ID to which the patient belongs.}
  \item{Pat.ID}{The name of the variable in \code{Dataset} that contains the patient's ID.}
  \item{Model}{The type of model that should be fitted, i.e., \code{Model=c("Full")}, \code{Model=c("Reduced")}, or \code{Model=c("SemiReduced")}. See the \bold{Details} section below. Default \code{Model=c("Full")}. }
\item{Weighted}{Logical. In practice it is often the case that different trials (or other clustering units) have different sample sizes. Univariate models are used to assess surrogacy in the information-theoretic approach, so it can be useful to adjust for heterogeneity in information content between the trial-specific contributions (particularly when trial-level surrogacy measures are of primary interest and when the heterogeneity in sample sizes is large). If \code{Weighted=TRUE}, weighted regression models are fitted. If \code{Weighted=FALSE}, unweighted regression analyses are conducted. See the \bold{Details} section below. Default \code{TRUE}.}
\item{Min.Trial.Size}{The minimum number of patients that a trial should contain to be included in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded from the analysis. Default \eqn{2}.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{h}} and \eqn{R^2_{ht}}. Default \eqn{0.05}.}
  \item{\dots}{Other arguments to be passed to the function \code{lmer} (of the R package \code{lme4}) that is used to fit the geralized linear mixed-effect models in the function \code{BimixedContCont}. For details, see \href{https://CRAN.R-project.org/package=lme4}{the lme4 manual}. }
}

\details{

\emph{Individual-level surrogacy}

The following generalised linear mixed-effect models are fitted:

\deqn{g_{T}(E(T_{ij}))=\mu_{T}+m_{Ti}+\beta Z_{ij}+b_{i}Z_{ij},}
\deqn{g_{T}(E(T_{ij}|S_{ij}))=\theta_{0}+c_{Ti}+\theta_{1}Z_{ij}+a_{i}Z_{ij}+\theta_{2i}S_{ij},}

where \eqn{i} and \eqn{j} are the trial and subject indicators, \eqn{g_{T}} is an appropriate link function (i.e., an identity link when a continuous true endpoint is considered), \eqn{S_{ij}} and \eqn{T_{ij}} are the surrogate and true endpoint values of subject \eqn{j} in trial \eqn{i}, and \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}. \eqn{\mu_{T}} and \eqn{\beta} are a fixed intercept and a fixed treatment-effect on the true endpoint, while \eqn{m_{Ti}} and \eqn{b_{i}} are the corresponding random effects. \eqn{\theta_{0}} and \eqn{\theta_{1}} are the fixed intercept and the fixed treatment effect on the true endpoint after accounting for the effect of the surrogate endpoint, and \eqn{c_{Ti}} and \eqn{a_i} are the corresponding random effects. 

The \eqn{-2} log likelihood values of the previous models (i.e., \eqn{L_{1}} and \eqn{L_{2}}, respectively) are subsequently used to compute individual-level surrogacy (based on the so-called Variance Reduction Factor, VFR; for details, see Alonso & Molenberghs, 2007):

\deqn{R^2_{hind}= 1 - exp \left(-\frac{L_{2}-L_{1}}{N} \right),}

where \eqn{N} is the number of trials.



\emph{Trial-level surrogacy}

When a full or semi-reduced model is requested (by using the argument \code{Model=c("Full")} or \code{Model=c("SemiReduced")} in the function call), trial-level surrogacy is assessed by fitting the following mixed models: 

\deqn{S_{ij}=\mu_{S}+m_{Si}+(\alpha+a_{i})Z_{ij}+\varepsilon_{Sij},   (1)}
\deqn{T_{ij}=\mu_{T}+m_{Ti}+(\beta+b_{i})Z_{ij}+\varepsilon_{Tij},   (1)} 

where \eqn{i} and \eqn{j} are the trial and subject indicators, \eqn{S_{ij}} and \eqn{T_{ij}} are the surrogate and true endpoint values of subject \eqn{j} in trial \eqn{i}, \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}, \eqn{\mu_{S}} and \eqn{\mu_{T}} are the fixed intercepts for S and T, \eqn{m_{Si}} and \eqn{m_{Ti}} are the corresponding random intercepts,  \eqn{\alpha} and \eqn{\beta} are the fixed treatment effects on S and T, and \eqn{a_{i}} and \eqn{b_{i}} are the corresponding random effects. The error terms \eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}} are assumed to be independent.

When a reduced model is requested by the user (by using the argument \code{Model=c("Reduced")} in the function call), the following univariate models are fitted:

\deqn{S_{ij}=\mu_{S}+(\alpha+a_{i})Z_{ij}+\varepsilon_{Sij},   (2)}
\deqn{T_{ij}=\mu_{T}+(\beta+b_{i})Z_{ij}+\varepsilon_{Tij},   (2)}

where \eqn{\mu_{S}} and \eqn{\mu_{T}} are the common intercepts for S and T. The other parameters are the same as defined above, and \eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}} are again assumed to be independent.

When the user requested that a full model approach is used (by using the argument \code{Model=c("Full")} in the function call, i.e., when models (1) were fitted), the following model is subsequently fitted:

\deqn{\widehat{\beta}_{i}=\lambda_{0}+\lambda_{1}\widehat{\mu_{Si}}+\lambda_{2}\widehat{\alpha_i}+\varepsilon_{i},   (3)}

where the parameter estimates for \eqn{\beta_i}, \eqn{\mu_{Si}}, and \eqn{\alpha_i} are based on models (1) (see above).  When a weighted model is requested (using the argument \code{Weighted=TRUE} in the function call), model (3) is a weighted regression model (with weights based on the number of observations in trial \eqn{i}). The \eqn{-2} log likelihood value of the (weighted or unweighted) models (3) (\eqn{L_1}) is subsequently compared to the \eqn{-2} log likelihood value of an intercept-only model (\eqn{\widehat{\beta}_{i}=\lambda_{3}}; \eqn{L_0}), and \eqn{R^2_{ht}} is computed based on the Variance Reduction Factor (VFR; for details, see Alonso & Molenberghs, 2007):

\deqn{R^2_{ht}= 1 - exp \left(-\frac{L_1-L_0}{N} \right),}

where \eqn{N} is the number of trials. 

When a semi-reduced or reduced model is requested (by using the argument \code{Model=c("SemiReduced")} or \code{Model=c("Reduced")} in the function call), the following model is fitted:

\deqn{\widehat{\beta_{i}}=\lambda_{0}+\lambda_{1}\widehat{\alpha_i}+\varepsilon_{i},}

where the parameter estimates for \eqn{\beta_i} and \eqn{\alpha_i} are based on models (2). The \eqn{-2} log likelihood value of this (weighted or unweighted) model (\eqn{L_1}) is subsequently compared to the \eqn{-2} log likelihood value of an intercept-only model (\eqn{\widehat{\beta}_{i}=\lambda_{3}}; \eqn{L_0}), and \eqn{R^2_{ht}} is computed based on the reduction in the likelihood (as described above).
 
}

\value{An object of class \code{MixedContContIT} with components,
  \item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of patients who have a missing value for the surrogate and/or the true endpoint are excluded. In addition, the data of trials (i) in which only one type of the treatment was administered, and (ii) in which either the surrogate or the true endpoint was a constant (i.e., all patients within a trial had the same surrogate and/or true endpoint value) are excluded. In addition, the user can specify the minimum number of patients that a trial should contain in order to include the trial in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded. \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
  \item{Obs.Per.Trial}{A \code{data.frame} that contains the total number of patients per trial and the number of patients who were administered the control treatment and the experimental treatment in each of the trials (in \code{Data.Analyze}).}
  \item{Trial.Spec.Results}{A \code{data.frame} that contains the trial-specific intercepts and treatment effects for the surrogate and the true endpoints (when a full or semi-reduced model is requested), or the trial-specific treatment effects for the surrogate and the true endpoints (when a reduced model is requested).}
  \item{R2ht}{A \code{data.frame} that contains the trial-level surrogacy estimate and its confidence interval.}
  \item{R2h.ind}{A \code{data.frame} that contains the individual-level surrogacy estimate and its confidence interval.}
  \item{Cor.Endpoints}{A \code{data.frame} that contains the correlations between the surrogate and the true endpoint in the control treatment group (i.e., \eqn{\rho_{T0S0}}) and in the experimental treatment group (i.e., \eqn{\rho_{T1S1}}), their standard errors and their confidence intervals.}
  \item{Residuals}{A \code{data.frame} that contains the residuals for the surrogate and true endpoints (\eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}}) that are obtained when models (1) or models (2) are fitted (see the \bold{Details} section above).}
}
\references{
Alonso, A, & Molenberghs, G. (2007). Surrogate marker evaluation from an information theory perspective. \emph{Biometrics, 63,} 180-186.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{FixedContContIT}}, \code{\link{plot Information-Theoretic}}
}
\examples{
# Example 1
# Based on the ARMD data:
data(ARMD)
# Assess surrogacy based on a full mixed-effect model
# in the information-theoretic framework:
Sur <- MixedContContIT(Dataset=ARMD, Surr=Diff24, True=Diff52, Treat=Treat, Trial.ID=Center,
Pat.ID=Id, Model="Full")
# Obtain a summary of the results:
summary(Sur)

\dontrun{ # Time consuming (>5sec) code part
# Example 2
# Conduct an analysis based on a simulated dataset with 2000 patients, 200 trials,
# and Rindiv=Rtrial=.8
# Simulate the data:
Sim.Data.MTS(N.Total=2000, N.Trial=200, R.Trial.Target=.8, R.Indiv.Target=.8,
Seed=123, Model="Full")
# Assess surrogacy based on a full mixed-effect model
# in the information-theoretic framework:
Sur2 <- MixedContContIT(Dataset=Data.Observed.MTS, Surr=Surr, True=True, Treat=Treat,
Trial.ID=Trial.ID, Pat.ID=Pat.ID, Model="Full")

# Show a summary of the results:
summary(Sur2)}
}

\keyword{Multiple-trial setting}
\keyword{Information-theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Likelihood Reduction Factor (LRF)}
\keyword{Mixed-effect models}
\keyword{Continuous endpoint}