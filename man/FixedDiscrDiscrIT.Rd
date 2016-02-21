\name{FixedDiscrDiscrIT}
\alias{FixedDiscrDiscrIT}
\title{Investigates surrogacy for binary or ordinal outcomes using the Information Theoretic framework}
\description{
  The function \code{FixedDiscrDiscrIT} uses the information theoretic
  approach (Alonso and Molenberghs 2007) to estimate trial and individual level
  surrogacy based on fixed-effects models when the surrogate is binary and the 
  true outcome is ordinal, the converse case or when both outcomes are ordinal 
  (the user must specify which form the data is in). The user can specify whether
  a weighted or unweighted analysis is required at the trial level. The penalized
  likelihood approach of Firth (1993) is applied to resolve issues of separation in discrete outcomes for particular trials. Requires packages \code{OrdinalLogisticBiplot} and \code{logistf}.
}
\usage{
FixedDiscrDiscrIT(Dataset, Surr, True, Treat, Trial.ID, 
Weighted = TRUE, Setting = c("binord"))
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true outcome value, a treatment indicator and a trial ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate outcome values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true outcome values.}
  \item{Treat}{The name of the  in \code{Dataset} that contains the treatment group values, 0/1 or -1/+1 are recommended.}
  \item{Trial.ID}{The name of the variable in \code{Dataset} that contains the trial ID to which the patient belongs.}
  \item{Weighted}{Logical. In practice it is often the case that different trials (or other clustering units) have different sample sizes. Univariate models are used to assess surrogacy in the information-theoretic approach, so it can be useful to adjust for heterogeneity in information content between the trial-specific contributions (particularly when trial-level surrogacy measures are of primary interest and when the heterogeneity in sample sizes is large). If \code{Weighted=TRUE}, weighted regression models are fitted. If \code{Weighted=FALSE}, unweighted regression analyses are conducted. See the \bold{Details} section below. Default \code{TRUE}.}
 \item{Setting}{Specifies whether an ordinal or binary surrogate or true outcome are present in \code{Dataset}. \code{Setting=c("binord")} for a binary surrogate and ordinal true outcome, \code{Setting=c("ordbin")} for an ordinal surrogate and binary true outcome and \code{Setting=c("ordord")} where both outcomes are ordinal.}
}

\details{
  \emph{Individual level surrogacy}
  
  The following univariate logistic regression models are fitted when \code{Setting=c("ordbin")}:
    
    \deqn{logit(P(T_{ij}=1))=\mu_{Ti}+\beta_{i}Z_{ij},    (1)}
  \deqn{logit(P(T_{ij}=1|S_{ij}=s))=\gamma_{0i}+\gamma_{1i}Z_{ij}+\gamma_{2i}S_{ij},    (1)}
  where: \eqn{i} and \eqn{j} are the trial and subject indicators; \eqn{S_{ij}} and \eqn{T_{ij}} are the surrogate and true outcome values of subject \eqn{j} in trial \eqn{i}; and \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}; \eqn{\mu_{Ti}} and \eqn{\beta_{i}} are the trial-specific intercepts and treatment-effects on the true endpoint in trial \eqn{i}; and \eqn{\gamma_{0i}} and \eqn{\gamma_{1i}} are the trial-specific intercepts and treatment-effects on the true endpoint in trial \eqn{i} after accounting for the effect of the surrogate endpoint.
  The \eqn{-2} log likelihood values of the previous models in each of the \eqn{i} trials (i.e., \eqn{L_{1i}} and \eqn{L_{2i}}, respectively) are subsequently used to compute individual-level surrogacy based on the  so-called Likelihood Reduction Factor (LRF; for details, see Alonso & Molenberghs, 2006):
    \deqn{R^2_{h}= 1 - \frac{1}{N} \sum_{i} exp \left(-\frac{L_{2i}-L_{1i}}{n_{i}} \right),      }
  where \eqn{N} is the number of trials and \eqn{n_{i}} is the number of patients within trial \eqn{i}.
  
  At the individual level in the discrete case \eqn{R^2_{h}} is bounded above by a number strictly less than one and is re-scaled (see Alonso & Molenberghs (2007)):
    \deqn{\widehat{R^2_{h}}= \frac{R^2_{h}}{1-e^{-2L_{0}}},}
  where \eqn{L_{0}} is the log-likelihood of the intercept only model of the true outcome (\eqn{logit(P(T_{ij}=1)=\gamma_{3}}).
    
    In the case of \code{Setting=c("binord")} or \code{Setting=c("ordord")} proportional odds models in (1) are used to accommodate the ordinal true response outcome, in all other respects the calculation of \eqn{R^2_{h}} would proceed in the same manner.
    
    \emph{Trial-level surrogacy}
    
    When  \code{Setting=c("ordbin")} trial-level surrogacy is assessed by fitting the following univariate logistic regression and proportional odds models for the ordinal surrogate and binary true response variables regressed on treatment for each trial \eqn{i}: 
      \deqn{logit(P(S_{ij} \leq W))=\mu_{S_{wi}}+\alpha_{i}Z_{ij},   (2)} 
    \deqn{logit(P(T_{ij}=1))=\mu_{Ti}+\beta_{i}Z_{ij},   (2)} 
    where: \eqn{i} and \eqn{j} are the trial and subject indicators; \eqn{S_{ij}} and \eqn{T_{ij}} are the surrogate and true outcome values of subject \eqn{j} in trial \eqn{i}; \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}; \eqn{\mu_{S_{wi}}} are the trial-specific intercept values for each cut point \eqn{w}, where \eqn{w=1,..,W-1}, of the ordinal surrogate outcome; \eqn{\mu_{Ti}} are the fixed trial-specific intercepts for T; and \eqn{\alpha_{i}} and \eqn{\beta_{i}} are the fixed trial-specific treatment effects on S and T, respectively. The mean trial-specific intercepts for the surrogate are calculated, \eqn{\overline{\mu}_{S_{wi}}}.The following model is subsequently fitted:
      
      \deqn{\widehat{\beta}_{i}=\lambda_{0}+\lambda_{1}\widehat{\overline{\mu}}_{S_{wi}}+\lambda_{2}\widehat{\alpha}_{i}+\varepsilon_{i},   (3)}
    
    where the parameter estimates for \eqn{\beta_i}, \eqn{\overline{\mu}_{S_{wi}}}, and \eqn{\alpha_i} are based on models (2) (see above).  When a weighted model is requested (using the argument \code{Weighted=TRUE} in the function call), model (2) is a weighted regression model (with weights based on the number of observations in trial \eqn{i}). The \eqn{-2} log likelihood value of the (weighted or unweighted) model (2) (\eqn{L_1}) is subsequently compared to the \eqn{-2} log likelihood value of an intercept-only model (\eqn{\widehat{\beta}_{i}=\lambda_{3}}; \eqn{L_0}), and \eqn{R^2_{ht}} is computed based on the Likelihood Reduction Factor (for details, see Alonso & Molenberghs, 2006):
      
      \deqn{R^2_{ht}= 1 - exp \left(-\frac{L_1-L_0}{N} \right),}
    
    where \eqn{N} is the number of trials. 
    
    When separation (the presence of zero cells) occurs in the cross tabs of treatment and the true or surrogate outcome for a particular trial in models (2) extreme bias can occur in \eqn{R^2_{ht}}. Under separation there are no unique maximum likelihood for parameters \eqn{\beta_i}, \eqn{\overline{\mu}_{S_{wi}}} and \eqn{\alpha_i}, in (2), for the affected trial \eqn{i}. This typically leads to extreme bias in the estimation of these parameters and hence outlying influential points in model (3), bias in \eqn{R^2_{ht}} inevitably follows. 
    
    To resolve the issue of separation the penalized likelihood approach of Firth (1993) is applied. This approach adds an asymptotically negligible component to the score function to allow unbiased estimation of \eqn{\beta_i}, \eqn{\overline{\mu}_{S_{wi}}}, and \eqn{\alpha_i} and in turn \eqn{R^2_{ht}}. The penalized likelihood R function \code{logitf} from the package of the same name is applied in the case of binary separation (Heinze and Schemper, 2002). The function \code{pordlogistf} from the package \code{OrdinalLogisticBioplot} is applied in the case of ordinal separation (Hern'andez, 2013). All instances of separation are reported.

In the case of \code{Setting=c("binord")} or \code{Setting=c("ordord")} the appropriate models (either logistic regression or a proportional odds models) are fitted in (2) to accommodate the form (either binary or ordinal) of the true or surrogate response variable. The rest of the analysis would proceed in a similar manner as that described above.
}
\value{An object of class \code{FixedDiscrDiscrIT} with components, 
\item{Trial.Spec.Results}{A \code{data.frame} that contains the trial-specific intercepts and treatment effects for the surrogate and the true endpoints. Also, the number of observations per trial; whether the trial was able to be included in the analysis for both \eqn{R^2_{h}} and \eqn{R^2_{ht}}; whether separation occurred and hence the penalized likelihood approach used for the surrogate or true outcome.}
\item{R2ht}{A \code{data.frame} that contains the trial-level surrogacy estimate and its confidence interval.}
\item{R2h}{A \code{data.frame} that contains the individual-level surrogacy estimate and its confidence interval.}
}
\references{
Alonso, A, & Molenberghs, G. (2007). Surrogate marker evaluation from an information theory perspective. \emph{Biometrics, 63,} 180-186.

Alonso, A, & Molenberghs, G., Geys, H., Buyse, M. & Vangeneugden, T. (2006). A unifying approach for surrogate marker validation based on Prentice's criteria. \emph{Statistics in medicine, 25,} 205-221.

Firth, D. (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika, 80,} 27-38.

Heinze, G. & Schemper, M. 2002. A solution to the problem of separation in logistic regression.  \emph{Statistics in medicine, 21,} 2409-2419.

Hern'andez, J. C. V.-V. O., J. L. 2013. OrdinalLogisticBiplot: Biplot representations of ordinal variables. R.
}
\author{
Hannah M. Ensor & Christopher J. Weir
}
\seealso{
\code{\link{FixedContContIT}}, \code{\link{plot Information-Theoretic}}, \code{\link{logistf}}
}
\examples{

\dontrun{ # Time consuming (>5sec) code part
# Example 1
# Conduct an analysis based on a simulated dataset with 2000 patients, 100 trials,
# and Rindiv=Rtrial=.8

# Simulate the data:
Sim.Data.MTS(N.Total=2000, N.Trial=100, R.Trial.Target=.8, R.Indiv.Target=.8,
Seed=123, Model="Full")

# create a binary true and ordinal surrogate outcome
Data.Observed.MTS$True<-findInterval(Data.Observed.MTS$True, 
c(quantile(Data.Observed.MTS$True,0.5)))
Data.Observed.MTS$Surr<-findInterval(Data.Observed.MTS$Surr, 
c(quantile(Data.Observed.MTS$Surr,0.333),quantile(Data.Observed.MTS$Surr,0.666)))

# Assess surrogacy based on a full fixed-effect model
# in the information-theoretic framework for a binary surrogate and ordinal true outcome:
SurEval <- FixedDiscrDiscrIT(Dataset=Data.Observed.MTS, Surr=Surr, True=True, Treat=Treat,
Trial.ID=Trial.ID, Setting="ordbin")

# Show a summary of the results:
summary(SurEval)
SurEval$Trial.Spec.Results
SurEval$R2h
SurEval$R2ht
}
}
\keyword{Multiple-trial setting}
\keyword{Information-theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Likelihood Reduction Factor (LRF)}
\keyword{Fixed-effect models}