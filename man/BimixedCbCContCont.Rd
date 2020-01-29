\name{BimixedCbCContCont}
\alias{BimixedCbCContCont}
\title{
Fits a bivariate mixed-effects model using the cluster-by-cluster (CbC) estimator to assess surrogacy in the meta-analytic multiple-trial setting  (Continuous-continuous case)}
\description{
The function \code{BimixedCbCContCont} uses the cluster-by-cluster (CbC) estimator of the bivariate mixed-effects to estimate trial- and individual-level surrogacy when the data of multiple clinical trials are available. See the \bold{Details} section below.  
}
\usage{
BimixedCbCContCont(Dataset, Surr, True, Treat, Trial.ID,Min.Treat.Size=2,Alpha=0.05)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, a treatment indicator, and a trial ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate endpoint values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should either be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group, or as \eqn{1} for the experimental group and \eqn{0} for the control group.}
  \item{Trial.ID}{The name of the variable in \code{Dataset} that contains the trial ID to which the patient belongs.}
  \item{Min.Treat.Size}{The minimum number of patients in each group (control or experimental) that a trial should contain to be included in the analysis. If the number of patients in a group of a trial is smaller than the value specified by \code{Min.Treat.Size}, the data of the trial are excluded from the analysis. Default \eqn{2}.}
    \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}} and \eqn{R^2_{indiv}}. Default \eqn{0.05}.}
}
\details{
The function \code{BimixedContCont} fits a bivariate mixed-effects model using the CbC estimator (for details, see Florez et al., 2019) to assess surrogacy (for details, see Buyse et al., 2000). In particular, the following mixed-effects model is fitted: 

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
The CbC estimator for the full bivariate mixed-effects model is closed-form (for details, see Florez et al., 2019). Therefore, it is fast. Furthermore, it is recommended when computational issues occur with the full maximum likelihood estimator (implemented in function \code{BimixedContCont}). 

The CbC estimator is performed in two stages: (1) a linear model is fitted in each trial. Evidently, it is require that the design matrix (\eqn{X_i}) is full column rank within each trial, allowing estimation of the fixed effects. When \eqn{X_i} is not full rank, trial i is excluded from the analysis. (2) a global estimator of the fixed effects (\eqn{\beta}) is obtained by weighted averaging the sets of estimates of each trial, and \eqn{\bold{D}} is estimated using a method-of-moments estimator. Optimal weights (for details, see Molenberghs et al., 2018) are used as a weighting scheme. 

The estimator of \eqn{\bold{D}} might lead to a non-positive-definite solution. Therefore, the eigenvalue method (for details, see Rousseeuw and Molenberghs, 1993) is used for non-positive-definiteness adjustment. 
}

\value{An object of class \code{BimixedContCont} with components,
  \item{Obs.Per.Trial}{A \code{data.frame} that contains the total number of patients per trial and the number of patients who were administered the control treatment and the experimental treatment in each of the trials (after excluding clusters). Clusters are excluded for two reasons: (i) the number of patients is smaller than the value especified by \code{Min.Trial.Size}, and (ii) the design matrix (\eqn{X_i}) is not full rank.}
  \item{Trial.removed}{Number of trials excluded from the analysis}
\item{Fixed.Effects}{A \code{data.frame} that contains the fixed intercept and treatment effects for the true and the surrogate endpoints (i.e., \eqn{\mu_{S}}, \eqn{\mu_{T}}, \eqn{\alpha}, and \eqn{\beta}) and their corresponding standard error.}
  \item{Trial.R2}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval.}
  \item{Indiv.R2}{A \code{data.frame} that contains the individual-level coefficient of determination (\eqn{R^2_{indiv}}), its standard error and confidence interval.}
\item{D}{The variance-covariance matrix of the random effects (the \eqn{\bold{D}} matrix), i.e., a \eqn{4} by \eqn{4} variance-covariance matrix of the random intercept and treatment effects.}
\item{DH.pd}{ \code{DH.pd=TRUE} if an adjustment for non-positive definiteness was not needed to estimate \eqn{\bold{D}}. \code{DH.pd=FALSE} if this adjustment was required.}
  \item{Sigma}{The \eqn{2} by \eqn{2} variance-covariance matrix of the residuals (\eqn{\varepsilon_{Sij}} and \eqn{\varepsilon_{Tij}}).}
}  

\references{
Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.

Florez, A. J., Molenberghs G, Verbeke G, Alonso, A. (2019). A closed-form estimator for meta-analysis
and surrogate markers evaluation. \emph{Journal of Biopharmaceutical Statistics, 29(2)} 318-332.

Molenberghs, G., Hermans, L., Nassiri, V., Kenward, M., Van der Elst, W., Aerts, M. and Verbeke, G. (2018). Clusters with random size: maximum likelihood versus weighted estimation. \emph{Statistica Sinica,
28}, 1107-1132.

Rousseeuw, P. J. and Molenberghs, G. (1993) Transformation of non positive semidefinite correlation
matrices. \emph{Communications in Statistics, Theory and Methods, 22}, 965-984.
}
\author{
Alvaro J. Florez, Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\code{\link{BimixedContCont}}, \code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, 
\code{\link{UnimixedContCont}}}
\examples{
# Open the Schizo dataset (clinial trial in schizophrenic patients)
data(Schizo)

# Fit a full bivariate random-effects model by the cluster-by-cluster (CbC) estimator
# a minimum of 2 subjects per group are allowed in each trial
  fit <- BimixedCbCContCont(Dataset=Schizo, Surr=BPRS, True=PANSS, Treat=Treat,Trial.ID=InvestId,
                              Alpha=0.05, Min.Treat.Size = 10)
# Note that an adjustment for non-positive definiteness was requiered and 113 trials were removed.

# Obtain a summary of the results
summary(fit)
                         
}

\keyword{Multiple-trial setting}
\keyword{Meta-analytic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Mixed-effect models}