\name{Single.Trial.ContCont}
\alias{Single.Trial.ContCont}
\alias{plot.Single.Trial.ContCont}
\title{Conducts a surrogacy analysis in the single-trial setting for Continuous S and T}
\description{
The function \code{Single.Trial.ContCont} conducts a surrogacy analysis in the single-trial setting for Continuous S and T. See \bold{Details} below. 
}

\usage{
Single.Trial.ContCont(Dataset, Surr, True, Treat, Pat.ID, Alpha=.05, 
Number.Bootstraps=500, Seed=12345)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, a treatment indicator, and a patient ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should either be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group, or as \eqn{1} for the experimental group and \eqn{0} for the control group. The \eqn{-1/1} coding is recommended. }
  \item{Pat.ID}{The name of the variable in \code{Dataset} that contains the patient's ID.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around Alpha (which is a parameter estimate of a model where the surrogate is regressed on the treatment indicator, see \bold{Details} below), Beta, RE, and \eqn{\gamma}. Default \eqn{0.05}.}  
\item{Number.Bootstraps}{The number of bootstrap samples that are used to obtain the bootstrap-based confidence intervals for RE, PE and the adjusted association (\eqn{\gamma}). Default \eqn{500}. A non-parametric bootstrap is used.}
  \item{Seed}{The seed that is used to generate the bootstrap samples. Default \code{Seed=12345}.}
}
\details{
\bold{Prentice criteria}
The Prentice criteria are examined by fitting the following regression models:
\deqn{S_{j}=\mu_{S}+\alpha Z_{j}+\varepsilon_{Sj}, (1)} 
\deqn{T_{j}=\mu_{T}+\beta Z_{j}+\varepsilon_{Tj}, (2)} 
\deqn{T_{j}=\mu+\gamma Z_{j}+\varepsilon_{j}, (3)} 
\deqn{T_{j}=\tilde{\mu}_{T}+\beta_{S} Z_{j}+\gamma_{Z} S_{j}+\tilde{\varepsilon}_{Tj}, (4)} 

where the error terms of (1) and (2) have a joint zero-mean normal distribution with variance-covariance matrix

\deqn{\boldsymbol{\Sigma}=\left(\begin{array}{cc}
\sigma_{SS}\\
\sigma_{ST} & \sigma_{TT}
\end{array}\right)},

and where \eqn{j} is the subject indicator, \eqn{S_{j}} and \eqn{T_{j}} are the surrogate and true endpoint values of subject \eqn{j}, and \eqn{Z_{j}} is the treatment indicator for subject \eqn{j}. 

To be in line with the Prentice criteria, Z should have a significant effect on S in model 1 (Prentice criterion 1), Z should have a significant effect on T in model 2 (Prentice criterion 2), S should have a significant effect on T in model 3 (Prentice criterion criterion 3), and the effect of Z on T should be fully captured by S in model 4 (Prentice criterion 4).  
 
\bold{Proportion Explained}
In view of the problems with the Prentice criteria (see the Limitations section below), Freedman et al. (1992) proposed to quantify surrogacy as the proportion of the effect of the treatment on T that is explained by S (the proportion explained, \eqn{PE})
\deqn{PE(T,\: S,\: Z)=\dfrac{\beta-\beta_{S}}{\beta}=1-\dfrac{\beta_{S}}{\beta},}
where \eqn{\beta} is the estimate of the effect of the treatment on T without correction for S and \eqn{\beta_{S}} is the estimate of the effect of the treatment on T with correction for S (see above). The intuition behind the PE is that, if all treatment effect is mediated by S (i.e., if \eqn{\beta_{S}=0}), then PE=1. On the other hand, if there is no mediation at all (i.e., if \eqn{\beta=\beta_{S}}), then PE=0. Note that the fourth Prentice criterion is equivalent to the requirement that PE=1.   


\bold{Relative Effect and Adjusted Association}
The Relative Effect (RE) and the adjusted association (\eqn{\gamma}) are based on the bivariate regression model (1) and (2) and its variance-covariance matrix 
\eqn{\boldsymbol{\Sigma}}. The parameter estimates of the fitted regression model and the variance-covariance matrix of the residuals are used to compute RE and the adjusted association (\eqn{\gamma}), respectively:

\deqn{RE=\frac{\beta}{\alpha},}
\deqn{\gamma=\frac{\sigma_{ST}}{\sqrt{\sigma_{SS}\sigma_{TT}}}.}


\bold{Limitations}
The Prentice approach to assess surrogacy has some fundamental limitations. For example, the fourth Prentice criterion requires that the statistical test for the  \eqn{\beta_S} in model 4 is non-significant. This criterion is useful to reject a poor surrogate, but it is not suitable to validate a good surrogate (i.e., a non-significant result may always be attributable to a lack of statistical power). Even when lack of power would not be an issue, the result of the statistical test to evaluate the fourth Prentice criterion cannot prove that the effect of the treatment on the true endpoint is fully captured by the surrogate. 

In general, the single-trial framework is hampered by a number of issues (Burzykowski et al., 2005). For example, a key motivation to validate a surrogate endpoint is to be able to predict the effect of Z on T as based on the effect of Z on S in a new clinical trial where T is not (yet) observed. The RE allows for such a prediction, but this requires the assumption that the relation between \eqn{\alpha} and \eqn{\beta} can be described by a linear regression model that goes through the origin. In other words, it has to be assumed that the RE remains constant across clinical trials. The constant RE assumption is unverifiable in a single-trial setting, but a way out of this problem is to combine the information of multiple clinical trials and generalize the RE concept to a multiple-trial setting (as is done in the multiple-trial meta-analytic approach, see \code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}}, and \code{\link{BimixedContCont}}).  
}

\value{An object of class \code{Single.Trial.ContCont} with components,
\item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of patients who have a missing value for the surrogate and/or the true endpoint are excluded. \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
\item{Prentice.Model.1}{An object of class \code{lm} that contains the fitted model 1 (using the Prentice approach).}
\item{Prentice.Model.2}{An object of class \code{lm} that contains the fitted model 2 (using the Prentice approach).}
\item{Prentice.Model.3}{An object of class \code{lm} that contains the fitted model 3 (using the Prentice approach).}
\item{Prentice.Model.4}{An object of class \code{lm} that contains the fitted model 4 (using the Prentice approach).}
\item{Prentice.Passed}{\code{Logical}. If all four Prentice criteria are fulfilled, \code{Prentice.Passed=TRUE}. If at least one criterion is not fulfilled, \code{Prentice.Passed=FALSE}.}

\item{Alpha}{An object of class \code{data.frame} that contains the parameter estimate for \eqn{\alpha}, its standard error, and its confidence interval. Note that \code{Alpha} is not to be confused with the \code{Alpha} argument in the function call, which specifies the \eqn{\alpha}-level of the confidence intervals of the parameters. }
\item{Beta}{An object of class \code{data.frame} that contains the parameter estimate for \eqn{\beta}, its standard error, and its confidence interval.}
\item{RE.Delta}{An object of class \code{data.frame} that contains the estimated RE, its standard error, and its confidence interval (based on the Delta method).}
%\item{RE.Fieller}{An object of class \code{data.frame} that contains the estimated RE, its standard error, and its confidence interval (based on Fieller's theorem).}
\item{RE.Boot}{An object of class \code{data.frame} that contains the estimated RE, its standard error, and its confidence interval (based on a non-parametric bootstra).} % Note that the occurence of outliers in the sample of bootstrapped RE values may lead to standard errors and/or confidence intervals that are not trustworthy. Such problems mainly occur when the parameter estimate for \eqn{\alpha} is close to 0 (taking its standard error into account). To detect possible outliers, studentized deleted residuals are computed (by fitting an intercept-only model with the bootstrapped RE values as the outcome variable). Bootstrapped RE values with an absolute studentized residual larger than \eqn{t(1-\alpha/2n;n-2)} are marked as outliers (where n = the number of bootstrapped RE values; Kutner et al., 2005). A warning is given when outliers are found, and the position of the outlier(s) in the bootstrap sample is identified. Inspection of the vector of bootstrapped RE values (see \code{RE.Boot.Samples} below) is recommended in this situation, and/or the use of the confidence intervals that are based on the Delta method or Fieller's theorem (rather than the bootstrap-based confidence interval).}
\item{AA}{An object of class \code{data.frame} that contains the adjusted association (i.e., \eqn{\gamma}), its standard error, and its confidence interval (based on the Fisher-Z transformation procedure).}
\item{AA.Boot}{An object of class \code{data.frame} that contains the adjusted association (i.e., \eqn{\gamma}), its standard error, and its confidence interval (based on a bootstrap procedure).}
\item{RE.Boot.Samples}{A vector that contains the RE values that were generated during the bootstrap procedure.}
\item{AA.Boot.Samples}{A vector that contains the adjusted association (i.e., \eqn{\gamma}) values that were generated during the bootstrap procedure.}
\item{Cor.Endpoints}{A \code{data.frame} that contains the correlations between the surrogate and the true endpoint in the control treatment group (i.e., \eqn{\rho_{T0T1}}) and in the experimental treatment group (i.e., \eqn{\rho_{T1S1}}), their standard errors and their confidence intervals.}
\item{Residuals}{A \code{data.frame} that contains the residuals for the surrogate and true endpoints that are obtained when the surrogate and the true endpoint are regressed on the treatment indicator.}
\item{Beta_S}{An object of class \code{data.frame} that contains the parameter estimate for \eqn{\beta_S}, its standard error, and its confidence interval.}
\item{PE.Delta}{An object of class \code{data.frame} that contains the estimated PE, its standard error, and its confidence interval (based on the Delta method).}
%\item{PE.Fieller}{An object of class \code{data.frame} that contains the estimated PE and its confidence interval (based on Fieller's method).}
\item{PE.Boot}{An object of class \code{data.frame} that contains the estimated PE, its standard error, and its confidence interval (based on a bootstrap).}
\item{PE.Boot.Samples}{A vector that contains the PE values that were generated during the bootstrap procedure.}
\item{PE.Fieller}{An object of class \code{data.frame} that contains the estimated PE and its confidence interval (based on the Fieller approach).} 
\item{RE.Fieller}{An object of class \code{data.frame} that contains the estimated RE and its confidence interval (based on the Fieller approach).} 
}

\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Buyse, M., & Molenberghs, G. (1998). The validation of surrogate endpoints in randomized experiments. \emph{Biometrics, 54,} 1014-1029.

Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.

}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}}, \code{\link{BimixedContCont}}, \code{\link{ICA.ContCont}} 
}
\examples{
\dontrun{ # time consuming code part
# Example 1, based on the ARMD data:
data(ARMD)

# Assess surrogacy based on the single-trial meta-analytic approach:
Sur <- Single.Trial.ContCont(Dataset=ARMD, Surr=Diff24, True=Diff52, Treat=Treat, Pat.ID=Id)

# Obtain a summary and plot of the results
summary(Sur)
plot(Sur)


# Example 2
# Conduct an analysis based on a simulated dataset with 2000 patients 
# and Rindiv=.8
# Simulate the data:
Sim.Data.STS(N.Total=2000, R.Indiv.Target=.8, Seed=123)

# Assess surrogacy:
Sur2 <- Single.Trial.ContCont(Dataset=Data.Observed.STS, Surr=Surr, True=True, Treat=Treat, 
Pat.ID=Pat.ID)

# Show a summary and plots of results
summary(Sur2)
plot(Sur2)
}
}
\keyword{Single-trial setting}
\keyword{Meta-analytic framework}
\keyword{Individual-level surrogacy}
\keyword{Trial-level surrogacy}
\keyword{Relative effect}
\keyword{Adjusted Association}