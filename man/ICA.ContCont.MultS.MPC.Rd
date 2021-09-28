\name{ICA.ContCont.MultS.MPC}
\alias{ICA.ContCont.MultS.MPC}
\title{Assess surrogacy in the causal-inference single-trial setting (Individual Causal Association, ICA) using a continuous univariate T and multiple continuous S, by simulating correlation matrices using a modified algorithm based on partial correlations}
\description{
The function \code{ICA.ContCont.MultS.MPC} quantifies surragacy in the single-trial causal-inference framework in which the true endpoint (T) and multiple surrogates (S) are continuous. This function is a modification of the \code{\link{ICA.ContCont.MultS.PC}} algorithm based on partial correlations. it mitigates the effect of non-informative surrogates and effectively explores the PD space to capture the ICA range (Florez, et al. 2021).}

\usage{ICA.ContCont.MultS.MPC(M=1000,N,Sigma,prob = NULL,Seed=123,
Save.Corr=F, Show.Progress=FALSE)}
\arguments{
  \item{M}{The number of multivariate ICA values (\eqn{R^2_{H}}) that should be sampled. Default \code{M=1000}.}
  \item{N}{The sample size of the dataset.}
  \item{Sigma}{A matrix that specifies the variance-covariance matrix between \eqn{T_0}, \eqn{T_1}, \eqn{S_{10}}, \eqn{S_{11}}, \eqn{S_{20}}, \eqn{S_{21}}, ..., \eqn{S_{k0}}, and \eqn{S_{k1}} (in this order, the \eqn{T_0} and \eqn{T_1} data should be in \code{Sigma[c(1,2), c(1,2)]}, the \eqn{S_{10}} and \eqn{S_{11}} data should be in \code{Sigma[c(3,4), c(3,4)]}, and so on). The unidentifiable covariances should be defined as \code{NA} (see example below).}
  \item{prob}{vector of probabilities to choose the number of surrogates (r) with their non-identifiable correlations equal to zero. The default (\code{prob=NULL}) vector of probabilities is:
\deqn{
\pi_{r} = \frac{{p \choose r}}{\sum_{i=1}^{p}{p \choose i}}, \mbox{ for }r=0,\ldots,p.
}
  In this way, each possible combination of $r$ surrogates has the same probability of being selected.}
  \item{Save.Corr}{If true, the lower diagonal elements of the correlation matrix (identifiable and unidientifiable elements) are stored. If false, these results are not saved.}
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
  
  The identifiable correlations are fixed at their estimated values and the unidentifiable correlations are independently and randomly sampled using a modification of an algorithm based on partial correlations (PC), called modified partial correlation (MPC) algorithm. In the function call, the unidentifiable correlations are marked by specifying \code{NA} in the \code{Sigma} matrix (see example section below). 
  
  
  The PC algorithm generate each correlation matrix progressively based on parameterization of terms of the correlations \eqn{\rho_{i,i+1}}, for \eqn{i=1,\ldots,d-1}, and the partial correlations \eqn{\rho_{i,j|i+1,\ldots,j-1}}, for \eqn{j-i>2} (for details, see Joe, 2006 and Florez et al., 2018). The MPC algorithm randomly fixed some of the unidentifiable correlations to zero in order to explore the PD, which is coherent with the estimable entries of the correlation matrix, to capture the ICA range more efficiently.
  
    Based on the identifiable variances, these correlation matrices are converted to covariance matrices \eqn{\boldsymbol{\Sigma}} and the multiple-surrogate ICA are estimated (for details, see Van der Elst et al., 2017).

This approach to simulate the unidentifiable parameters of \eqn{\boldsymbol{\Sigma}} is computationally more efficient than the one used in the function \code{\link{ICA.ContCont.MultS}}. 
}
\value{An object of class \code{ICA.ContCont.MultS.PC} with components,
\item{R2_H}{The multiple-surrogate individual causal association value(s).}
\item{Corr.R2_H}{The corrected multiple-surrogate individual causal association value(s).}
\item{Lower.Dig.Corrs.All}{A \code{data.frame} that contains the matrix that contains the identifiable and unidentifiable correlations (lower diagonal elements) that were used to compute (\eqn{R^2_{H}}) in the run.}
\item{surr.eval.r}{Matrix indicating the surrogates of which their unidentifiable correlations are fixed to zero in each simulation.}
}
\references{
Florez, A., Molenberghs, G., Van der Elst, W., Alonso, A. A. (2021). An efficient algorithm for causally assessing surrogacy in a multivariate setting.

Florez, A., Alonso, A. A., Molenberghs, G. & Van der Elst, W. (2020). Generating random correlation matrices with fixed values: An application to the evaluation of multivariate surrogate endpoints. \emph{Computational Statistics & Data Analysis 142}.

Joe, H. (2006). Generating random correlation matrices based on partial correlations. \emph{Journal of Multivariate Analysis, 97(10)}:2177-2189.

Van der Elst, W., Alonso, A. A., & Molenberghs, G. (2017). Univariate versus multivariate surrogate endpoints.

}
\author{
Wim Van der Elst, Ariel Alonso, Geert Molenberghs & Alvaro Florez
}
\seealso{
\code{\link{MICA.ContCont}}, \code{\link{ICA.ContCont}}, \code{\link{Single.Trial.RE.AA}}, 
\code{\link{plot Causal-Inference ContCont}}, \code{\link{ICA.ContCont.MultS}}, \code{\link{ICA.ContCont.MultS_alt}}
}
\examples{
\dontrun{ 
# Specify matrix Sigma (var-cavar matrix T_0, T_1, S1_0, S1_1, ...)
# here we have 1 true endpoint and 10 surrogates (8 of these are non-informative)

Sigma = ks::invvech(
  c(25, NA, 17.8, NA, -10.6, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 
    4, NA, -0.32, NA, -1.32, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, 16, 
    NA, -4, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 1, NA, 0.48, NA, 
    0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, 16, NA, 0, NA, 0, NA, 0, NA, 0, 
    NA, 0, NA, 0, NA, 0, NA, 0, NA, 1, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, NA, 0, 
    NA, 0, 16, NA, 8, NA, 8, NA, 8, NA, 8, NA, 8, NA, 8, NA, 8, NA, 1, NA, 0.5, NA, 0.5, 
    NA, 0.5, NA, 0.5, NA, 0.5, NA, 0.5, NA, 0.5, 16, NA, 8, NA, 8, NA, 8, NA, 8, NA, 8, 
    NA, 8, NA, 1, NA, 0.5, NA, 0.5, NA, 0.5, NA, 0.5, NA, 0.5, NA, 0.5, 16, NA, 8, NA, 
    8, NA, 8, NA, 8, NA, 8, NA, 1,NA,0.5,NA,0.5,NA,0.5,NA,0.5,NA,0.5, 16, NA, 8, NA, 8, 
    NA, 8, NA, 8, NA, 1, NA, 0.5, NA, 0.5, NA, 0.5, NA, 0.5, 16, NA, 8, NA, 8, NA, 8, NA,
    1, NA, 0.5, NA, 0.5, NA, 0.5, 16, NA, 8, NA, 8, NA, 1, NA, 0.5, NA, 0.5, 16, NA, 8, NA,
    1, NA, 0.5, 16, NA, 1)) 

# Conduct analysis using the PC and MPC algorithm 
## first evaluating two surrogates
ICA.PC.2 = ICA.ContCont.MultS.PC(M = 30000, N=200, Sigma[1:6,1:6], Seed = 123) 
ICA.MPC.2 = ICA.ContCont.MultS.MPC(M = 30000, N=200, Sigma[1:6,1:6],prob=NULL, 
Seed = 123, Save.Corr=T, Show.Progress = TRUE) 

## later evaluating two surrogates
ICA.PC.10 = ICA.ContCont.MultS.PC(M = 150000, N=200, Sigma, Seed = 123) 
ICA.MPC.10 = ICA.ContCont.MultS.MPC(M = 150000, N=200, Sigma,prob=NULL, 
Seed = 123, Save.Corr=T, Show.Progress = TRUE) 

# Explore results
range(ICA.PC.2$R2_H)
range(ICA.PC.10$R2_H)

range(ICA.MPC.2$R2_H)
range(ICA.MPC.10$R2_H)
## as we observe, the MPC algorithm displays a wider interval of possible values for the ICA
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{Multiple surrogates}
\keyword{Multivariate}
