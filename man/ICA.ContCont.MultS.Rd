\name{ICA.ContCont.MultS}
\alias{ICA.ContCont.MultS}
\title{Assess surrogacy in the causal-inference single-trial setting (Individual Causal Association, ICA) using a continuous univariate T and multiple continuous S}
\description{
The function \code{ICA.ContCont.MultS} quantifies surrogacy in the single-trial causal-inference framework  where T is continuous and there are multiple continuous S.}
\usage{ICA.ContCont.MultS(M = 500, N, Sigma, 
G = seq(from=-1, to=1, by = .00001),
Seed=c(123), Show.Progress=FALSE)}
\arguments{
  \item{M}{The number of multivariate ICA values (\eqn{R^2_{H}}) that should be sampled. Default \code{M=500}.}
  \item{N}{The sample size of the dataset.}
  \item{Sigma}{A matrix that specifies the variance-covariance matrix between \eqn{T_0}, \eqn{T_1}, \eqn{S_{10}}, \eqn{S_{11}}, \eqn{S_{20}}, \eqn{S_{21}}, ..., \eqn{S_{k0}}, and \eqn{S_{k1}} (in this order, the \eqn{T_0} and \eqn{T_1} data should be in \code{Sigma[c(1,2), c(1,2)]}, the \eqn{S_{10}} and \eqn{S_{11}} data should be in \code{Sigma[c(3,4), c(3,4)]}, and so on). The unidentifiable covariances should be defined as \code{NA} (see example below).}
  \item{G}{A vector of the values that should be considered for the unidentified correlations. Default \code{G=seq(-1, 1, by=.00001)}, i.e., values with range \eqn{-1} to \eqn{1}.}
  \item{Seed}{The seed that is used. Default \code{Seed=123}.}
  \item{Show.Progress}{Should progress of runs be graphically shown? (i.e., 1\% done..., 2\% done..., etc). Mainly useful when a large number of S have to be considered (to follow progress and estimate total run time).}
}
\details{
The multivariate ICA (\eqn{R^2_{H}}) is not identifiable because the individual causal treatment effects on \eqn{T}, \eqn{S_1}, ..., \eqn{S_k} cannot be observed. A simulation-based sensitivity analysis is therefore conducted in which the multivariate ICA (\eqn{R^2_{H}}) is estimated across a set of plausible values for the unidentifiable correlations. To this end, consider the variance covariance matrix of the potential outcomes \eqn{\boldsymbol{\Sigma}} (0 and 1 subscripts refer to the control and experimental treatments, respectively): 
\deqn{\boldsymbol{\Sigma} = \left(\begin{array}{ccccccccc}
\sigma_{T_{0}T_{0}}\\
\sigma_{T_{0}T_{1}} & \sigma_{T_{1}T_{1}}\\
\sigma_{T_{0}S1_{0}} & \sigma_{T_{1}S1_{0}} & \sigma_{S1_{0}S1_{0}}\\
\sigma_{T_{0}S1_{1}} & \sigma_{T_{1}S1_{1}} & \sigma_{S1_{0}S1_{1}} & \sigma_{S1_{1}S1_{1}}\\
\sigma_{T_{0}S2_{0}} & \sigma_{T_{1}S2_{0}} & \sigma_{S1_{0}S2_{0}} & \sigma_{S1_{1}S2_{0}} & \sigma_{S2_{0}S2_{0}}\\
\sigma_{T_{0}S2_{1}} & \sigma_{T_{1}S2_{1}} & \sigma_{S1_{0}S2_{1}} & \sigma_{S1_{1}S2_{1}} & \sigma_{S2_{0}S2_{1}} & \sigma_{S2_{1}S2_{1}}\\
... & ... & ... & ... & ... & ... & \ddots\\
\sigma_{T_{0}Sk_{0}} & \sigma_{T_{1}Sk_{0}} & \sigma_{S1_{0}Sk_{0}} & \sigma_{S1_{1}Sk_{0}} & \sigma_{S2_{0}Sk_{0}} & \sigma_{S2_{1}Sk_{0}} & ... & \sigma_{Sk_{0}Sk_{0}}\\
\sigma_{T_{0}Sk_{1}} & \sigma_{T_{1}Sk_{1}} & \sigma_{S1_{0}Sk_{1}} & \sigma_{S1_{1}Sk_{1}} & \sigma_{S2_{0}Sk_{1}} & \sigma_{S2_{1}Sk_{1}} & ... & \sigma_{Sk_{0}Sk_{1}} & \sigma_{Sk_{1}Sk_{1}}.
\end{array}\right)}

The \code{ICA.ContCont.MultS} function requires the user to specify a distribution \eqn{G} for the unidentified correlations. Next, the identifiable correlations are fixed at their estimated values and the unidentifiable correlations are independently and randomly sampled from \eqn{G}. In the function call, the unidentifiable correlations are marked by specifying \code{NA} in the \code{Sigma} matrix (see example section below). The algorithm generates a large number of 'completed' matrices, and only those that are positive definite are retained (the number of positive definite matrices that should be obtained is specified by the \code{M=} argument in the function call). Based on the identifiable variances, these positive definite correlation matrices are converted to covariance matrices \eqn{\boldsymbol{\Sigma}} and the multiple-surrogate ICA are estimated.

An issue with this approach (i.e., substituting unidentified correlations by random and independent samples from \eqn{G}) is that the probability of obtaining a positive definite matrix is very low when the dimensionality of the matrix increases. One approach to increase the efficiency of the algorithm is to build-up the correlation matrix in a gradual way. In particular, the property that a \eqn{\left(k \times k\right)} matrix is positive definite if and only if all principal minors are positive (i.e., Sylvester's criterion) can be used. In other words, a \eqn{\left(k \times k\right)} matrix is positive definite when the determinants of the upper-left \eqn{\left(2 \times 2\right)}, \eqn{\left(3 \times 3\right)}, ..., \eqn{\left(k \times k\right)} submatrices all have a positive determinant. Thus, when a positive definite \eqn{\left(k \times k\right)} matrix has to be generated, one can start with the upper-left \eqn{\left(2 \times 2\right)} submatrix and randomly sample a value from the unidentified correlation (here: \eqn{\rho_{T_0T_0}}) from \eqn{G}. When the determinant is positive (which will always be the case for a \eqn{\left(2 \times 2\right)} matrix), the same procedure is used for the upper-left \eqn{\left(3 \times 3\right)} submatrix, and so on. When a particular draw from \eqn{G} for a particular submatrix does not give a positive determinant, new values are sampled for the unidentified correlations until a positive determinant is obtained. In this way, it can be guaranteed that the final \eqn{\left(k \times k\right)} submatrix will be positive definite. The latter approach is used in the current function. This procedure is used to generate many positive definite matrices. Based on these matrices, \eqn{\boldsymbol{\Sigma_{\Delta}}} is generated and the multivariate ICA (\eqn{R^2_{H}}) is computed (for details, see Van der Elst et al., 2017). 
}
\value{An object of class \code{ICA.ContCont.MultS} with components,
\item{R2_H}{The multiple-surrogate individual causal association value(s).}
\item{Corr.R2_H}{The corrected multiple-surrogate individual causal association value(s).}
}
\references{
Van der Elst, W., Alonso, A. A., & Molenberghs, G. (2017). Univariate versus multivariate surrogate endpoints.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{MICA.ContCont}}, \code{\link{ICA.ContCont}}, \code{\link{Single.Trial.RE.AA}}, 
\code{\link{plot Causal-Inference ContCont}}, \code{\link{ICA.ContCont.MultS_alt}}
}
\examples{
\dontrun{ #time-consuming code parts
# Specify matrix Sigma (var-cavar matrix T_0, T_1, S1_0, S1_1, ...)
# here for 1 true endpoint and 3 surrogates

s<-matrix(rep(NA, times=64),8)
s[1,1] <- 450; s[2,2] <- 413.5; s[3,3] <- 174.2; s[4,4] <- 157.5; 
s[5,5] <- 244.0; s[6,6] <- 229.99; s[7,7] <- 294.2; s[8,8] <- 302.5
s[3,1] <- 160.8; s[5,1] <- 208.5; s[7,1] <- 268.4 
s[4,2] <- 124.6; s[6,2] <- 212.3; s[8,2] <- 287.1
s[5,3] <- 160.3; s[7,3] <- 142.8 
s[6,4] <- 134.3; s[8,4] <- 130.4 
s[7,5] <- 209.3; 
s[8,6] <- 214.7 
s[upper.tri(s)] = t(s)[upper.tri(s)]

# Marix looks like (NA indicates unidentified covariances):
#            T_0    T_1  S1_0  S1_1  S2_0   S2_1  S2_0  S2_1
#            [,1]  [,2]  [,3]  [,4]  [,5]   [,6]  [,7]  [,8]
# T_0  [1,] 450.0    NA 160.8    NA 208.5     NA 268.4    NA
# T_1  [2,]    NA 413.5    NA 124.6    NA 212.30    NA 287.1
# S1_0 [3,] 160.8    NA 174.2    NA 160.3     NA 142.8    NA
# S1_1 [4,]    NA 124.6    NA 157.5    NA 134.30    NA 130.4
# S2_0 [5,] 208.5    NA 160.3    NA 244.0     NA 209.3    NA
# S2_1 [6,]    NA 212.3    NA 134.3    NA 229.99    NA 214.7
# S3_0 [7,] 268.4    NA 142.8    NA 209.3     NA 294.2    NA
# S3_1 [8,]    NA 287.1    NA 130.4    NA 214.70    NA 302.5

# Conduct analysis
ICA <- ICA.ContCont.MultS(M=100, N=200, Show.Progress = TRUE,
  Sigma=s, G = seq(from=-1, to=1, by = .00001), Seed=c(123))

# Explore results
summary(ICA)
plot(ICA)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{Multiple surrogates}
\keyword{Multivariate}
