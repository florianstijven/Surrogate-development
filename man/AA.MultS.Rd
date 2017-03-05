\name{AA.MultS}
\alias{AA.MultS}
\title{Compute the multiple-surrogate adjusted association}
\description{
The function \code{AA.MultS} computes the multiple-surrogate adjusted correlation. This is a generalisation of the adjusted association proposed by Buyse & Molenberghs (1998) (see \code{\link{Single.Trial.RE.AA}}) to the setting where there are multiple endpoints. See \bold{Details} below. 
}
\usage{
AA.MultS(Sigma_gamma, N, Alpha=0.05)
}
\arguments{
  \item{Sigma_gamma}{The variance covariance matrix of the residuals of regression models in which the true endpoint (\eqn{T}) is regressed on the treatment (\eqn{Z}), the first surrogate (\eqn{S1}) is regressed on \eqn{Z}, ..., and the \eqn{k}-th surrogate (\eqn{Sk}) is regressed on \eqn{Z}. See \bold{Details} below.}
\item{N}{The sample size (needed to compute a CI around the multiple adjusted association; \eqn{\gamma_M})}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence interval around \eqn{\gamma_M}. Default \eqn{0.05}.}  
}
\details{
The multiple-surrogate adjusted association (\eqn{\gamma_M}) is obtained by regressing \eqn{T}, \eqn{S1}, \eqn{S2}, ..., \eqn{Sk} on the treatment (\eqn{Z}):
\deqn{T_{j}=\mu_{T}+\beta Z_{j}+\varepsilon_{Tj},}
\deqn{S1_{j}=\mu_{S1}+\alpha_{1}Z_{j}+\varepsilon_{S1j},}
\deqn{\ldots,}
\deqn{Sk_{j}=\mu_{Sk}+\alpha_{k}Z_{j}+\varepsilon_{Skj},}
where the error terms have a joint zero-mean normal distribution with variance-covariance matrix:

\deqn{{\boldsymbol{\Sigma}=\left(\begin{array}{cc}
\sigma_{TT} & \Sigma_{\boldsymbol{S}T}\\
\Sigma^{'}_{\boldsymbol{S}T} & \Sigma_{\boldsymbol{SS}} \\
\end{array}\right).}}

The multiple adjusted association is then computed as
\deqn{\gamma_M = \sqrt(\frac{\left(\Sigma^{'}_{ST} \Sigma^{-1}_{SS} \Sigma_{ST}\right)}{\sigma_{TT}})}
}

\value{An object of class \code{AA.MultS} with components,
\item{Gamma.Delta}{An object of class \code{data.frame} that contains the multiple-surrogate adjusted association (i.e., \eqn{\gamma_M}), its standard error, and its confidence interval (based on the Fisher-Z transformation procedure).} 
\item{Corr.Gamma.Delta}{An object of class \code{data.frame} that contains the bias-corrected multiple-surrogate adjusted association (i.e., corrected \eqn{\gamma_M}), its standard error, and its confidence interval (based on the Fisher-Z transformation procedure).} 
  \item{Sigma_gamma}{The variance covariance matrix of the residuals of regression models in which \eqn{T} is regressed on \eqn{Z}, \eqn{S1} is regressed on \eqn{Z}, ..., and \eqn{Sk} is regressed on \eqn{Z}.}
\item{N}{The sample size (used to compute a CI around the multiple adjusted association; \eqn{\gamma_M})}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence interval around \eqn{\gamma_M}.}  
}
\references{
Buyse, M., & Molenberghs, G. (1998). The validation of surrogate endpoints in randomized experiments. \emph{Biometrics, 54,} 1014-1029.

Van der Elst, W., Alonso, A. A., & Molenberghs, G. (2017). A causal inference-based approach to evaluate surrogacy using multiple surrogates.}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{Single.Trial.RE.AA}} 
}
\examples{
data(ARMD.MultS)

# Regress T on Z, S1 on Z, ..., Sk on Z 
# (to compute the covariance matrix of the residuals)
Res_T <- residuals(lm(Diff52~Treat, data=ARMD.MultS))
Res_S1 <- residuals(lm(Diff4~Treat, data=ARMD.MultS))
Res_S2 <- residuals(lm(Diff12~Treat, data=ARMD.MultS))
Res_S3 <- residuals(lm(Diff24~Treat, data=ARMD.MultS))
Residuals <- cbind(Res_T, Res_S1, Res_S2, Res_S3)

# Make covariance matrix of residuals, Sigma_gamma
Sigma_gamma <- cov(Residuals)

# Conduct analysis
Result <- AA.MultS(Sigma_gamma = Sigma_gamma, N = 188, Alpha = .05)

# Explore results
summary(Result)
}
\keyword{Adjusted Association}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{Multiple surrogates}