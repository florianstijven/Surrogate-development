\name{BimixedContCont}
\alias{BimixedContCont}
\title{
Fits a bivariate mixed-effects model to assess surrogacy in the meta-analytic multiple-trial setting 
(Continuous-continuous case)}
\description{
The function \code{BimixedContCont} uses the bivariate mixed-effects approach to estimate trial- and individual-level surrogacy when the data of multiple clinical trials are available. The user can specify whether a full or reduced model should be fitted. See the \bold{Details} section below.   
}
\usage{
BimixedContCont(Dataset, Surr, True, Treat, Trial.ID, Pat.ID, Model=c("Full"), 
Min.Trial.Size=2, Alpha=.05, \dots)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, a treatment indicator, a patient ID, and a trial ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate endpoint values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should either be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group, or as \eqn{1} for the experimental group and \eqn{0} for the control group.}
  \item{Trial.ID}{The name of the variable in \code{Dataset} that contains the trial ID to which the patient belongs.}
  \item{Pat.ID}{The name of the variable in \code{Dataset} that contains the patient's ID.}
  \item{Model}{The type of model that should be fitted, i.e., \code{Model=c("Full")} or \code{Model=c("Reduced")}. See the \bold{Details} section below. 
  Default \code{Model=c("Full")}. }
  \item{Min.Trial.Size}{The minimum number of patients that a trial should contain to be included in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded from the analysis. Default \eqn{2}.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}}, \eqn{R_{trial}}, \eqn{R^2_{indiv}} and \eqn{R_{indiv}}. Default \eqn{0.05}.}  
  \item{\dots}{Other arguments to be passed to the function \code{lmer} (of the R package \code{lme4}) that is used to fit the geralized linear mixed-effect models in the function \code{BimixedContCont}. For details, see \href{https://CRAN.R-project.org/package=lme4}{the lme4 manual}. }
}
\details{
The function \code{BimixedContCont} fits a bivariate mixed-effects model to assess surrogacy (for details, see Buyse et al., 2000). In particular, the following mixed-effects model is fitted: 

\deqn{S_{ij}=\mu_{S}+m_{Si}+(\alpha+a_{i})Z_{ij}+\varepsilon_{Sij},}
\deqn{T_{ij}=\mu_{T}+m_{Ti}+(\beta+b_{i})Z_{ij}+\varepsilon_{Tij},} 

where \eqn{i} and \eqn{j} are the trial and subject indicators, \eqn{S_{ij}} and \eqn{T_{ij}} are the surrogate and true endpoint values of subject \eqn{j} in trial \eqn{i}, \eqn{Z_{ij}} is the treatment indicator for subject \eqn{j} in trial \eqn{i}, \eqn{\mu_{S}} and \eqn{\mu_{T}} are the fixed intercepts for S and T, \eqn{m_{Si}} and \eqn{m_{Ti}} are the corresponding random intercepts, \eqn{\alpha} and \eqn{\beta} are the fixed treatment effects for S and T, and \eqn{a_{i}} and \eqn{b_{i}} are the corresponding random treatment effects, respectively. 

The vector of the random effects (i.e., \eqn{m_{Si}}, \eqn{m_{Ti}}, \eqn{a_{i}} and \eqn{b_{i}}) is assumed to be mean-zero normally distributed with variance-covariance matrix \eqn{\bold{D}}:

\deqn{\bold{D}=\left(\begin{array}{cccc}
d_{SS}\\
d_{ST} & d_{TT}\\
d_{Sa} & d_{Ta} & d_{aa}\\
d_{Sb} & d_{Tb} & d_{ab} & d_{bb}
\end{array}\right).}
 
The trial-level coefficient of determination (i.e., \eqn{R^2_{trial}}) is quantified as:

\deqn{R_{trial}^{2}=\frac{\left(\begin{array}{c}
d_{Sb}\\
d_{ab}
\end{array}\right)^{'}\left(\begin{array}{cc}
d_{SS} & d_{Sa}\\
d_{Sa} & d_{aa}
\end{array}\right)^{-1}\left(\begin{array}{c}
d_{Sb}\\
d_{ab}
\end{array}\right)}{d_{bb}}.}
 
The error terms \eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}} are assumed to be mean-zero normally distributed with variance-covariance matrix \eqn{\bold{\Sigma}}: 

\deqn{\bold{\Sigma}=\left(\begin{array}{cc}\sigma_{SS}\\\sigma_{ST} & \sigma_{TT}\end{array}\right).} 

Based on \eqn{\bold{\Sigma}}, individual-level surrogacy is quantified as: 

\deqn{R_{indiv}^{2}=\frac{\sigma_{ST}^{2}}{\sigma_{SS}\sigma_{TT}}.} 

\bold{\emph{Note}}

When the full bivariate mixed-effects approach is used to assess surrogacy in the meta-analytic framework (for details, see Buyse & Molenberghs, 2000), computational issues often occur. Such problems mainly occur when the number of trials is low, the number of patients in the different trials is low, and/or when the trial-level heterogeneity is small (Burzykowski et al., 2000). 

In that situation, the use of a simplified model-fitting strategy may be warranted (for details, see Burzykowski et al., 2000; Tibaldi et al., 2003).

For example, a reduced bivariate-mixed effect model can be fitted instead of a full model (by using the \code{Model=c("Reduced")} argument in the function call). In the reduced model, the random-effects structure is simplified (i) by assuming that there is no heterogeneity in the random intercepts, or (ii) by assuming that the covariance between the random intercepts and random treatment effects is zero. Note that under this assumption, the computation of the trial-level coefficient of determination (i.e., \eqn{R^2_{trial}}) simplifies to: 

\deqn{R_{trial}^{2}=\frac{d_{ab}^{2}}{d_{aa}d_{bb}}.}
 
Alternatively, the bivariate mixed-effects model may be abandonned and the user may fit a univariate fixed-effects model, a bivariate fixed-effects model, or a univariate mixed-effects model (for details, see Tibaldi et al., 2003). These models are implemented in the functions \code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, and \code{\link{UnimixedContCont}}).
}

\value{An object of class \code{BimixedContCont} with components,
  \item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of patients who have a missing value for the surrogate and/or the true endpoint are excluded. In addition, the data of trials (i) in which only one type of the treatment was administered, and (ii) in which either the surrogate or the true endpoint was a constant (i.e., all patients within a trial had the same surrogate and/or true endpoint value) are excluded. In addition, the user can specify the minimum number of patients that a trial should contain in order to include the trial in the analysis. If the number of patients in a trial is smaller than the value specified by \code{Min.Trial.Size}, the data of the trial are excluded. \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
  \item{Obs.Per.Trial}{A \code{data.frame} that contains the total number of patients per trial and the number of patients who were administered the control treatment and the experimental treatment in each of the trials (in \code{Data.Analyze}).}
  \item{Trial.Spec.Results}{A \code{data.frame} that contains the trial-specific intercepts and treatment effects on the surrogate and the true endpoints when a full model is requested (i.e., \eqn{\mu_{S}+m_{Si}}, \eqn{\mu_{T}+m_{Ti}}, \eqn{\alpha+a_{i}}, and \eqn{\beta+b_{i}}), or the trial-specific treatment effects on the surrogate and the true endpoints  when a reduced model is requested (i.e., \eqn{\alpha+a_{i}}, and \eqn{\beta+b_{i}}). Note that the results that are contained in \code{Trial.Spec.Results} are equivalent to the results in \code{Results.Stage.1} that are obtained when the functions \code{\link{UnifixedContCont}}, \code{\link{UnimixedContCont}}, or \code{\link{BifixedContCont}} are used.}
  \item{Residuals}{A \code{data.frame} that contains the residuals for the surrogate and true endpoints (\eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}}).}
\item{Fixed.Effect.Pars}{A \code{data.frame} that contains the fixed intercept and treatment effects for the surrogate and the true endpoints (i.e., \eqn{\mu_{S}}, \eqn{\mu_{T}}, \eqn{\alpha}, and \eqn{\beta}).}
\item{Random.Effect.Pars}{A \code{data.frame} that contains the random intercept and treatment effects for the surrogate and the true endpoints (i.e., \eqn{m_{Si}}, \eqn{m_{Ti}}, \eqn{a_{i}}, and \eqn{b_{i}}) when a full model is fitted (i.e., when \code{Model=c("Full")} is used in the function call), or that contains the random treatment effects for the surrogate and the true endpoints (i.e., \eqn{a_{i}} and \eqn{b_{i}}) when a reduced model is fitted (i.e., when \code{Model=c("Reduced")} is used in the function call).}  
  \item{Trial.R2}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval.}
  \item{Indiv.R2}{A \code{data.frame} that contains the individual-level coefficient of determination (\eqn{R^2_{indiv}}), its standard error and confidence interval.}
  \item{Trial.R}{A \code{data.frame} that contains the trial-level correlation coefficient (\eqn{R_{trial}}), its standard error and confidence interval.}
  \item{Indiv.R}{A \code{data.frame} that contains the individual-level correlation coefficient (\eqn{R_{indiv}}), its standard error and confidence interval.}
  \item{Cor.Endpoints}{A \code{data.frame} that contains the correlations between the surrogate and the true endpoint in the control treatment group (i.e., \eqn{\rho_{T0S0}}) and in the experimental treatment group (i.e., \eqn{\rho_{T1S1}}), their standard errors and their confidence intervals.}
\item{D}{The variance-covariance matrix of the random effects (the \eqn{\bold{D}} matrix), i.e., a \eqn{4} by \eqn{4} variance-covariance matrix of the random intercept and treatment effects when a full model is fitted (i.e., when \code{Model=c("Full")} is used in the function call), or a \eqn{2} by \eqn{2} variance-covariance matrix of the random treatment effects when a reduced model is fitted (i.e., when \code{Model=c("Reduced")} is used in the function call).}
  \item{Sigma}{The \eqn{2} by \eqn{2} variance-covariance matrix of the residuals (\eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}}).}
}  

\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.

Tibaldi, F., Abrahantes, J. C., Molenberghs, G., Renard, D., Burzykowski, T., Buyse, M., Parmar, M., et al., (2003). Simplified hierarchical linear models for the evaluation of surrogate endpoints. \emph{Journal of Statistical Computation and Simulation, 73,} 643-658.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}}, \code{\link{plot Meta-Analytic}}
}
\examples{
# Open the Schizo dataset (clinial trial in schizophrenic patients)
data(Schizo)

\dontrun{ #Time consuming (>5 sec) code part
# When a reduced bivariate mixed-effect model is used to assess surrogacy, 
# the conditioning number for the D matrix is very high: 
Sur <- BimixedContCont(Dataset=Schizo, Surr=BPRS, True=PANSS, Treat=Treat, Model="Reduced", 
Trial.ID=InvestId, Pat.ID=Id)

# Such problems often occur when the total number of patients, the total number 
# of trials and/or the trial-level heterogeneity
# of the treatment effects is relatively small

# As an alternative approach to assess surrogacy, consider using the functions
# BifixedContCont, UnifixedContCont or UnimixedContCont in the meta-analytic framework,
# or use the information-theoretic approach
}
}

\keyword{Multiple-trial setting}
\keyword{Meta-analytic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Mixed-effect models}

