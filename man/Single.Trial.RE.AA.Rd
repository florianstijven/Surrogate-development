\name{Single.Trial.RE.AA}
\alias{Single.Trial.RE.AA}
\alias{plot.Single.Trial.RE.AA}
\title{Conducts a surrogacy analysis based on the single-trial meta-analytic framework}
\description{
The function \code{Single.Trial.RE.AA} conducts a surrogacy analysis based on the single-trial meta-analytic framework of Buyse & Molenberghs (1998). See \bold{Details} below. 
}

\usage{
Single.Trial.RE.AA(Dataset, Surr, True, Treat, Pat.ID, Alpha=.05, 
Number.Bootstraps=500, Seed=sample(1:1000, size=1))
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, a treatment indicator, and a patient ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should either be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group, or as \eqn{1} for the experimental group and \eqn{0} for the control group. The \eqn{-1/1} coding is recommended. }
  \item{Pat.ID}{The name of the variable in \code{Dataset} that contains the patient's ID.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around Alpha (which is a parameter estimate of a model where the surrogate is regressed on the treatment indicator, see \bold{Details} below), Beta, RE, and \eqn{\gamma}. Default \eqn{0.05}.}  
\item{Number.Bootstraps}{The number of bootstrap samples that are used to obtain the bootstrapp-based confidence intervals for RE and the adjusted association (\eqn{\gamma}). Default \eqn{500}.}
  \item{Seed}{The seed that is used to generate the bootstrap samples. Default \code{sample(x=1:1000, size=1)}, i.e., a random number between 1 and 1000.}
}
\details{

The Relative Effect (RE) and the adjusted association (\eqn{\gamma}) are based on the following bivariate regression model (when the surrogate and the true endpoints are continuous variables):
\deqn{S_{j}=\mu_{S}+\alpha Z_{j}+\varepsilon_{Sj},} 
\deqn{T_{j}=\mu_{T}+\beta Z_{j}+\varepsilon_{Tj},} 

where the error terms have a joint zero-mean normal distribution with variance-covariance matrix:

\deqn{\boldsymbol{\Sigma}=\left(\begin{array}{cc}
\sigma_{SS}\\
\sigma_{ST} & \sigma_{TT}
\end{array}\right),}

and where \eqn{j} is the subject indicator, \eqn{S_{j}} and \eqn{T_{j}} are the surrogate and true endpoint values of patient \eqn{j}, and \eqn{Z_{j}} is the treatment indicator for patient \eqn{j}. 

The parameter estimates of the fitted regression model and the variance-covariance matrix of the residuals are used to compute RE and the adjusted association (\eqn{\gamma}), respectively:

\deqn{RE=\frac{\beta}{\alpha},}
\deqn{\gamma=\frac{\sigma_{ST}}{\sqrt{\sigma_{SS}\sigma_{TT}}}.}


\bold{Note}

The single-trial meta-analytic framework is hampered by a number of issues (Burzykowski et al., 2005). For example, a key motivation to validate a surrogate endpoint is to be able to predict the effect of Z on T as based on the effect of Z on S in a new clinical trial where T is not (yet) observed. The RE allows for such a prediction, but this requires the assumption that the relation between \eqn{\alpha} and \eqn{\beta} can be described by a linear regression model that goes through the origin. In other words, it has to be assumed that the RE remains constant across clinical trials. The constant RE assumption is unverifiable in a single-trial setting, but a way out of this problem is to combine the information of multiple clinical trials and generalize the RE concept to a multiple-trial setting (as is done in the multiple-trial meta-analytic approach, see \code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}}, and \code{\link{BimixedContCont}}).  
}

\value{An object of class \code{Single.Trial.RE.AA} with components,
  \item{Data.Analyze}{Prior to conducting the surrogacy analysis, data of patients who have a missing value for the surrogate and/or the true endpoint are excluded. \code{Data.Analyze} is the dataset on which the surrogacy analysis was conducted.}
\item{Alpha}{An object of class \code{data.frame} that contains the parameter estimate for \eqn{\alpha}, its standard error, and its confidence interval. Note that \code{Alpha} is not to be confused with the \code{Alpha} argument in the function call, which specifies the \eqn{\alpha}-level of the confidence intervals of the parameters. }
\item{Beta}{An object of class \code{data.frame} that contains the parameter estimate for \eqn{\beta}, its standard error, and its confidence interval.}
\item{RE.Delta}{An object of class \code{data.frame} that contains the estimated RE, its standard error, and its confidence interval (based on the Delta method).}
\item{RE.Fieller}{An object of class \code{data.frame} that contains the estimated RE, its standard error, and its confidence interval (based on Fieller's theorem).}
\item{RE.Boot}{An object of class \code{data.frame} that contains the estimated RE, its standard error, and its confidence interval (based on bootstrapping). Note that the occurence of outliers in the sample of bootstrapped RE values may lead to standard errors and/or confidence intervals that are not trustworthy. Such problems mainly occur when the parameter estimate for \eqn{\alpha} is close to 0 (taking its standard error into account). To detect possible outliers, studentized deleted residuals are computed (by fitting an intercept-only model with the bootstrapped RE values as the outcome variable). Bootstrapped RE values with an absolute studentized residual larger than \eqn{t(1-\alpha/2n;n-2)} are marked as outliers (where n = the number of bootstrapped RE values; Kutner et al., 2005). A warning is given when outliers are found, and the position of the outlier(s) in the bootstrap sample is identified. Inspection of the vector of bootstrapped RE values (see \code{RE.Boot.Samples} below) is recommended in this situation, and/or the use of the confidence intervals that are based on the Delta method or Fieller's theorem (rather than the bootstrap-based confidence interval).}
\item{AA}{An object of class \code{data.frame} that contains the adjusted association (i.e., \eqn{\gamma}), its standard error, and its confidence interval (based on the Fisher-Z transformation procedure).}
\item{AA.Boot}{An object of class \code{data.frame} that contains the adjusted association (i.e., \eqn{\gamma}), its standard error, and its confidence interval (based on a bootstrap procedure).}
\item{RE.Boot.Samples}{A vector that contains the RE values that were generated during the bootstrap procedure.}
\item{AA.Boot.Samples}{A vector that contains the adjusted association (i.e., \eqn{\gamma}) values that were generated during the bootstrap procedure.}
\item{Cor.Endpoints}{A \code{data.frame} that contains the correlations between the surrogate and the true endpoint in the control treatment group (i.e., \eqn{\rho_{T0T1}}) and in the experimental treatment group (i.e., \eqn{\rho_{T1S1}}), their standard errors and their confidence intervals.}
\item{Residuals}{A \code{data.frame} that contains the residuals for the surrogate and true endpoints that are obtained when the surrogate and the true endpoint are regressed on the treatment indicator.}

}

\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Buyse, M., & Molenberghs, G. (1998). The validation of surrogate endpoints in randomized experiments. \emph{Biometrics, 54,} 1014-1029.

Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.

Kutner, M. H., Nachtsheim, C. J., Neter, J., & Li, W. (2005). \emph{Applied linear statistical models (5th ed.)}. New York: McGraw Hill.
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
Sur <- Single.Trial.RE.AA(Dataset=ARMD, Surr=Diff24, True=Diff52, Treat=Treat, Pat.ID=Id)

# Obtain a summary and plot of the results
summary(Sur)
plot(Sur)


# Example 2
# Conduct an analysis based on a simulated dataset with 2000 patients 
# and Rindiv=.8
# Simulate the data:
Sim.Data.STS(N.Total=2000, R.Indiv.Target=.8, Seed=123)

# Assess surrogacy:
Sur2 <- Single.Trial.RE.AA(Dataset=Data.Observed.STS, Surr=Surr, True=True, Treat=Treat, 
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