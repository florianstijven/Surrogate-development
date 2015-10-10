\name{FixedBinBinIT}
\alias{FixedBinBinIT}
\title{Fits (univariate) fixed-effect models to assess surrogacy in the binary-binary case based on the Information-Theoretic framework}
\description{
The function \code{FixedBinBinIT} uses the information-theoretic approach (Alonso & Molenberghs, 2007) to estimate trial- and individual-level surrogacy based on fixed-effect models when both S and T are binary variables. The user can specify whether a (weighted or unweighted) full, semi-reduced, or reduced model should be fitted. See the \bold{Details} section below.   
}
\usage{
FixedBinBinIT(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, 
Model=c("Full"), Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, 
Number.Bootstraps=500, Seed=sample(1:1000, size=1))
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
\item{Number.Bootstraps}{The standard errors and confidence intervals for \eqn{R^2_{h}} and \eqn{R^2_{h.ind}} are determined based on a bootstrap procedure. \code{Number.Bootstraps} specifies the number of bootstrap samples that are used. Default \eqn{500}.}
\item{Seed}{The seed to be used in the bootstrap procedure. Default \eqn{sample(1:1000, size=1)}.}
}
\details{

\emph{Individual-level surrogacy}

The following univariate generalised linear models are fitted:

\deqn{g_{T}(E(T_{ij}))=\mu_{Ti}+\beta_{i}Z_{ij},}
\deqn{g_{T}(E(T_{ij}|S_{ij}))=\gamma_{0i}+\gamma_{1i}Z_{ij}+\gamma_{2i}S_{ij},}

where \eqn{i} and \eqn{j} are the trial and subject indicators, \eqn{g_{T}} is an appropriate link function (i.e., a logit link when binary endpoints are considered), \eqn{S_{ij}} and \eqn{T_{ij}} are the surrogate and true endpoint values of subject \eqn{j} in trial \eqn{i}, and \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}. \eqn{\mu_{Ti}} and \eqn{\beta_{i}} are the trial-specific intercepts and treatment-effects on the true endpoint in trial \eqn{i}. \eqn{\gamma_{0i}} and \eqn{\gamma_{1i}} are the trial-specific intercepts and treatment-effects on the true endpoint in trial \eqn{i} after accounting for the effect of the surrogate endpoint. 

The \eqn{-2} log likelihood values of the previous models in each of the \eqn{i} trials (i.e., \eqn{L_{1i}} and \eqn{L_{2i}}, respectively) are subsequently used to compute individual-level surrogacy based on the  so-called Variance Reduction Factor (VFR; for details, see Alonso & Molenberghs, 2007):

\deqn{R^2_{h}= 1 - \frac{1}{N} \sum_{i} exp \left(-\frac{L_{2i}-L_{1i}}{n_{i}} \right),}

where \eqn{N} is the number of trials and \eqn{n_{i}} is the number of patients within trial \eqn{i}.

When it can be assumed (i) that the treatment-corrected association between the surrogate and the true endpoint is constant across trials, or (ii) when all data come from a single clinical trial (i.e., when \eqn{N=1}), the previous expression simplifies to:

\deqn{R^2_{h.ind}= 1 - exp \left(-\frac{L_{2}-L_{1}}{N} \right).}

The upper bound does not reach to 1 when \eqn{T} is binary, i.e., its maximum is 0.75. Kent (1983) claims that 0.75 is a reasonable upper bound and thus \eqn{R^2_{h.ind}} can usually be interpreted without paying special consideration to the discreteness of \eqn{T}. Alternatively, to address the upper bound problem, a scaled version of the mutual information can be used when both \eqn{S} and \eqn{T} are binary (Joe, 1989):

\deqn{R^2_{b.ind}= \frac{I(T,S)}{min[H(T), H(S)]},}

where the entropy of \eqn{T} and \eqn{S} in the previous expression can be estimated using the log likelihood functions of the GLMs shown above. %%%

\emph{Trial-level surrogacy}

When a full or semi-reduced model is requested (by using the argument \code{Model=c("Full")} or \code{Model=c("SemiReduced")} in the function call), trial-level surrogacy is assessed by fitting the following univariate models: 

\deqn{S_{ij}=\mu_{Si}+\alpha_{i}Z_{ij}+\varepsilon_{Sij},   (1)}
\deqn{T_{ij}=\mu_{Ti}+\beta_{i}Z_{ij}+\varepsilon_{Tij},   (1)} 

where \eqn{i} and \eqn{j} are the trial and subject indicators, \eqn{S_{ij}} and \eqn{T_{ij}} are the surrogate and true endpoint values of subject \eqn{j} in trial \eqn{i}, \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}, \eqn{\mu_{Si}} and \eqn{\mu_{Ti}} are the fixed trial-specific intercepts for S and T, and \eqn{\alpha_{i}} and \eqn{\beta_{i}} are the fixed trial-specific treatment effects on S and T, respectively. The error terms \eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}} are assumed to be independent.

When a reduced model is requested by the user (by using the argument \code{Model=c("Reduced")} in the function call), the following univariate models are fitted:

\deqn{S_{ij}=\mu_{S}+\alpha_{i}Z_{ij}+\varepsilon_{Sij},   (2)}
\deqn{T_{ij}=\mu_{T}+\beta_{i}Z_{ij}+\varepsilon_{Tij},   (2)}

where \eqn{\mu_{S}} and \eqn{\mu_{T}} are the common intercepts for S and T. The other parameters are the same as defined above, and \eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}} are again assumed to be independent.

When the user requested a full model approach (by using the argument \code{Model=c("Full")} in the function call, i.e., when models (1) were fitted), the following model is subsequently fitted:

\deqn{\widehat{\beta}_{i}=\lambda_{0}+\lambda_{1}\widehat{\mu_{Si}}+\lambda_{2}\widehat{\alpha}_{i}+\varepsilon_{i},   (3)}

where the parameter estimates for \eqn{\beta_i}, \eqn{\mu_{Si}}, and \eqn{\alpha_i} are based on models (1) (see above).  When a weighted model is requested (using the argument \code{Weighted=TRUE} in the function call), model (3) is a weighted regression model (with weights based on the number of observations in trial \eqn{i}). The \eqn{-2} log likelihood value of the (weighted or unweighted) model (3) (\eqn{L_1}) is subsequently compared to the \eqn{-2} log likelihood value of an intercept-only model (\eqn{\widehat{\beta}_{i}=\lambda_{3}}; \eqn{L_0}), and \eqn{R^2_{ht}} is computed based based on the Variance Reduction Factor (for details, see Alonso & Molenberghs, 2007):

\deqn{R^2_{ht}= 1 - exp \left(-\frac{L_1-L_0}{N} \right),}

where \eqn{N} is the number of trials. 

When a semi-reduced or reduced model is requested (by using the argument \code{Model=c("SemiReduced")} or \code{Model=c("Reduced")} in the function call), the following model is fitted:

\deqn{\widehat{\beta}_{i}=\lambda_{0}+\lambda_{1}\widehat{\alpha}_{i}+\varepsilon_{i},}

where the parameter estimates for \eqn{\beta_i} and \eqn{\alpha_i} are based on models (1) when a semi-reduced model is fitted or on models (2) when a reduced model is fitted. The \eqn{-2} log likelihood value of this (weighted or unweighted) model (\eqn{L_1}) is subsequently compared to the \eqn{-2} log likelihood value of an intercept-only model (\eqn{\widehat{\beta}_{i}=\lambda_{3}}; \eqn{L_0}), and \eqn{R^2_{ht}} is computed based on the reduction in the likelihood (as described above).
 
}

\value{An object of class \code{FixedBinBinIT} with components,
  \item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of patients who have a missing value for the surrogate and/or the true endpoint are excluded. In addition, the data of trials (i) in which only one type of the treatment was administered, and (ii) in which either the surrogate or the true endpoint was a constant (i.e., all patients within a trial had the same surrogate and/or true endpoint value) are excluded. In addition, the user can specify the minimum number of patients that a trial should contain in order to include the trial in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded. \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
  \item{Obs.Per.Trial}{A \code{data.frame} that contains the total number of patients per trial and the number of patients who were administered the control treatment and the experimental treatment in each of the trials (in \code{Data.Analyze}).}
  \item{Trial.Spec.Results}{A \code{data.frame} that contains the trial-specific intercepts and treatment effects for the surrogate and the true endpoints (when a full or semi-reduced model is requested), or the trial-specific treatment effects for the surrogate and the true endpoints (when a reduced model is requested).}
  \item{R2ht}{A \code{data.frame} that contains the trial-level surrogacy estimate and its confidence interval.}
  \item{R2h.ind}{A \code{data.frame} that contains the individual-level surrogacy estimate \eqn{R^2_{h.ind}} (single-trial based estimate) and its confidence interval.}
  \item{R2h}{A \code{data.frame} that contains the individual-level surrogacy estimate \eqn{R^2_{h}} (cluster-based estimate) and its confidence interval (based on a bootsrtrap).}
  \item{R2b.ind}{A \code{data.frame} that contains the individual-level surrogacy estimate \eqn{R^2_{b.ind}} (single-trial based estimate accounting for upper bound) and its confidence interval (based on a bootstrap).}
  \item{R2h.Ind.By.Trial}{A \code{data.frame} that contains individual-level surrogacy estimates \eqn{R^2_{hInd}} (cluster-based estimates) and their confidence interval for each of the trials seperately.}
}
\references{
Alonso, A, & Molenberghs, G. (2007). Surrogate marker evaluation from an information theory perspective. \emph{Biometrics, 63,} 180-186.

Joe, H. (1989). Relative entropy measures of multivariate dependence. \emph{Journal of the American Statistical Association, 84,} 157-164.

Kent, T. J. (1983). Information gain as a general measure of correlation. \emph{Biometrica, 70,} 163-173.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{FixedBinContIT}}, \code{\link{FixedContBinIT}}, \code{\link{plot Information-Theoretic BinCombn}}
}
\examples{
\dontrun{ # Time consuming (>5sec) code part
# Generate data with continuous Surr and True
Sim.Data.MTS(N.Total=5000, N.Trial=50, R.Trial.Target=.9, R.Indiv.Target=.9,
             Fixed.Effects=c(0, 0, 0, 0), D.aa=10, D.bb=10, Seed=1,
             Model=c("Full"))
# Dichtomize Surr and True
Surr_Bin <- Data.Observed.MTS$Surr
Surr_Bin[Data.Observed.MTS$Surr>.5] <- 1
Surr_Bin[Data.Observed.MTS$Surr<=.5] <- 0
True_Bin <- Data.Observed.MTS$True
True_Bin[Data.Observed.MTS$True>.15] <- 1
True_Bin[Data.Observed.MTS$True<=.15] <- 0
Data.Observed.MTS$Surr <- Surr_Bin
Data.Observed.MTS$True <- True_Bin

# Assess surrogacy using info-theoretic framework
Fit <- FixedBinBinIT(Dataset = Data.Observed.MTS, Surr = Surr, 
True = True, Treat = Treat, Trial.ID = Trial.ID, 
Pat.ID = Pat.ID, Number.Bootstraps=100)

# Examine results
summary(Fit)
plot(Fit, Trial.Level = FALSE, Indiv.Level.By.Trial=TRUE)
plot(Fit, Trial.Level = TRUE, Indiv.Level.By.Trial=FALSE)
}
}
\keyword{plot Information-Theoretic BinBin}
\keyword{Multiple-trial setting}
\keyword{Information-theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Likelihood Reduction Factor (LRF)}
\keyword{Fixed-effect models}
\keyword{Binary endpoint}