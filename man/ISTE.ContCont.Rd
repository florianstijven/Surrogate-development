\name{ISTE.ContCont}
\alias{ISTE.ContCont}
\title{Individual-level surrogate threshold effect for continuous normally distributed surrogate and true endpoints.}
\description{Computes the individual-level surrogate threshold effect in the causal-inference single-trial setting where both the surrogate and the true endpoint are continuous normally distributed variables. For details, see paper in the references section.}
\usage{ISTE.ContCont(Mean_T1, Mean_T0, Mean_S1, Mean_S0, N, Delta_S=c(-10, 0, 10), 
zeta.PI=0.05, PI.Bound=0, PI.Lower=TRUE, Show.Prediction.Plots=TRUE, Save.Plots="No", 
T0S0, T1S1, T0T0=1, T1T1=1, S0S0=1, S1S1=1, T0T1=seq(-1, 1, by=.001), 
T0S1=seq(-1, 1, by=.001), T1S0=seq(-1, 1, by=.001),
S0S1=seq(-1, 1, by=.001), M.PosDef=500, Seed=123)}
\arguments{
\item{Mean_T1}{A scalar or vector that specifies the mean of the true endpoint in the experimental treatment condition (a vector is used to account for estimation uncertainty).}
\item{Mean_T0}{A scalar or vector that specifies the mean of the true endpoint in the control condition (a vector is used to account for estimation uncertainty).}
\item{Mean_S1}{A scalar or vector that specifies the mean of the surrogate endpoint in the experimental treatment condition (a vector is used to account for estimation uncertainty).}
\item{Mean_S0}{A scalar or vector that specifies the mean of the surrogate endpoint in the control condition (a vector is used to account for estimation uncertainty).}
\item{N}{The sample size of the clinical trial.}
\item{Delta_S}{The vector or scalar of \eqn{\Delta S} values for which the expected \eqn{\Delta T} and its prediction error has to be computed.}
\item{zeta.PI}{The alpha-level to be used in the computation of the prediction interval around \eqn{E(\Delta T)}. Default \code{zeta.PI=0.05}, i.e., the \eqn{95\%} prediction interval.}
\item{PI.Bound}{The ISTE is defined as the value of \eqn{\Delta S} for which the lower (or upper) bound of the \eqn{(1-\alpha)\%} prediction interval around \eqn{E(\Delta T)} is 0. If another threshold value than 0 is desired, this can be requested by using the \code{PI.Bound} argument. For example, the argument \code{PI.Bound=5} can be used in the function call to obtain the values of \eqn{\Delta S} for which the lower (or upper) bound of the \eqn{(1-\alpha)\%} prediction intervals (in the different runs of the algorithm)around \eqn{\Delta T} equal 5.}
\item{PI.Lower}{Logical. Should a lower (\code{PI.Lower=TRUE}) or upper (\code{PI.Lower=FALSE}) prediction interval be used in the computation of ISTE? Default \code{PI.Lower=TRUE}.}
\item{Show.Prediction.Plots}{Logical. Should plots that depict \eqn{E(\Delta T)} against \eqn{\Delta S} (prediction function), the prediction interval, and the ISTE for the different runs of the algorithm be shown? Default \code{Show.Prediction.Plots=TRUE}.}
\item{Save.Plots}{Should the prediction plots (see previous item) be saved? If \code{Save.Plots="No"} is used (the default argument), the plots are not saved. If the plots have to be saved, replace \code{"No"} by the desired location, e.g., \code{Save.Plots="C:/Analysis directory/"} on a windows computer or \code{Save.Plots="/Users/wim/Desktop/Analysis directory/"} on macOS or Linux.} 
\item{T0S0}{A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the control treatment condition that should be considered in the computation of ISTE.}
\item{T1S1}{A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the experimental treatment condition that should be considered in the computation of ISTE.}
\item{T0T0}{A scalar that specifies the variance of the true endpoint in the control treatment condition that should be considered in the computation of ISTE. Default 1.}
\item{T1T1}{A scalar that specifies the variance of the true endpoint in the experimental treatment condition that should be considered in the computation of ISTE. Default 1.}
\item{S0S0}{A scalar that specifies the variance of the surrogate endpoint in the control treatment condition that should be considered in the computation of ISTE. Default 1.}
\item{S1S1}{A scalar that specifies the variance of the surrogate endpoint in the experimental treatment condition that should be considered in the computation of ISTE. Default 1.}
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1 that should be considered in the computation of ISTE. Default \code{seq(-1, 1, by=.001)}.}
\item{T0S1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and S1 that should be considered in the computation of ISTE. Default \code{seq(-1, 1, by=.001)}.}
\item{T1S0}{A scalar or vector that contains the correlation(s) between the counterfactuals T1 and S0 that should be considered in the computation of ISTE. Default \code{seq(-1, 1, by=.001)}.}
\item{S0S1}{A scalar or vector that contains the correlation(s) between the counterfactuals S0 and S1 that should be considered in the computation of ISTE. Default \code{seq(-1, 1, by=.001)}.}
\item{M.PosDef}{The number of positive definite \eqn{\Sigma} matrices that should be identified. This will also determine the amount of ISTE values that are identified. Default \code{M.PosDef=500}.}
\item{Seed}{The seed to be used in the analysis (for reproducibility). Default \code{Seed=123}.}
}
\details{
See paper in the references section.
}
\value{An object of class \code{ICA.ContCont} with components,
%\item{ISTE_Low_PI_0}{The vector of individual surrogate threshold effect (ISTE) values, i.e., the values of \eqn{\Delta S} for which the lower bound of the \eqn{(1-\alpha)\%} prediction interval around \eqn{\Delta T} is 0.}  
%\item{ISTE_Up_PI_0}{The vector of individual surrogate threshold effect (ISTE) values, i.e., the values of \eqn{\Delta S} for which the upper bound of the \eqn{(1-\alpha)\%} prediction interval around \eqn{\Delta T} is 0.}  
\item{ISTE_Low_PI}{The vector of individual surrogate threshold effect (ISTE) values, i.e., the values of \eqn{\Delta S} for which the lower bound of the \eqn{(1-\alpha)\%} prediction interval around \eqn{\Delta T} is 0 (or another threshold value, which can be requested by using the \code{PI.Bound} argument in the function call).} 
\item{ISTE_Up_PI}{Same as \code{ISTE_Low_PI}, but using the upper bound of the \eqn{(1-\alpha)\%} prediction interval.} 
\item{MSE}{The vector of mean squared error values that are obtained in the prediction of \eqn{\Delta T} based on \eqn{\Delta S}.} 
\item{gamma0}{The vector of intercepts that are obtained in the prediction of \eqn{\Delta T} based on \eqn{\Delta S}.}
\item{gamma1}{The vector of slope that are obtained in the prediction of \eqn{\Delta T} based on \eqn{\Delta S}.}
\item{Delta_S_For_Which_Delta_T_equal_0}{The vector of \eqn{\Delta S} values for which \eqn{E(\Delta T = 0)}.}
\item{S_squared_pred}{The vector of variances of the prediction errors for \eqn{\Delta T}.}
\item{Predicted_Delta_T}{The vector/matrix of predicted values of \eqn{\Delta T} for the \eqn{\Delta S} values that were requested in the function call (argument \code{Delta_S}).}
\item{PI_Interval_Low}{The vector/matrix of lower bound values of the \eqn{(1-\alpha)\%} prediction interval around \eqn{\Delta T} for the \eqn{\Delta S} values that were requested in the function call (argument \code{Delta_S}).}
\item{PI_Interval_Up}{The vector/matrix of upper bound values of the \eqn{(1-\alpha)\%} prediction interval around \eqn{\Delta T} for the \eqn{\Delta S} values that were requested in the function call (argument \code{Delta_S}).}
\item{T0T0}{The vector of variances of T0 (true endpoint in the control treatment) that are used in the computation (this is a constant if the variance is fixed in the function call).}
\item{T1T1}{The vector of variances of T1 (true endpoint in the experimental treatment) that are used in the computations (this is a constant if the variance is fixed in the function call).}
\item{S0S0}{The vector of variances of S0 (surrogate endpoint in the control treatment) that are used in the computations (this is a constant if the variance is fixed in the function call).}
\item{S1S1}{The vector of variances of S1 (surrogate endpoint in the experimental treatment) that are used in the computations (this is a constant if the variance is fixed in the function call).}
\item{Mean_DeltaT}{The vector of treatment effect values on the true endpoint that are used in the computations (this is a constant if the means of T0 and T1 are fixed in the function call).} 
\item{Mean_DeltaS}{The vector of treatment effect values on the surrogate endpoint that are used in the computations (this is a constant if the means of S0 and S1 are fixed in the function call).} 
\item{Total.Num.Matrices}{An object of class \code{numeric} that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.}
\item{Pos.Def}{A \code{data.frame} that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the ISTE values.}
\item{ICA}{Apart from ISTE, ICA is also computed (the individual causal association). For details, see \code{\link{ICA.ContCont}}.}
\item{zeta.PI}{The \code{zeta.PI} value specified in the function call.}
\item{PI.Bound}{The \code{PI.Bound} value specified in the function call.}
\item{PI.Lower}{The \code{PI.Lower} value specified in the function call.}
\item{Delta_S}{The \code{Delta_S} value(s) specified in the function call.}
}
\references{
Van der Elst, W., Alonso, A. A., and Molenberghs, G. (submitted). The individual-level surrogate threshold effect in a causal-inference setting.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.ContCont}}
}
\examples{
# Define input for analysis using the Schizo dataset, 
# with S=BPRS and T = PANSS. 
# For each of the identifiable quantities,
# uncertainty is accounted for by specifying a uniform
# distribution with min, max values corresponding to
# the 95\% confidence interval of the quantity.
T0S0 <- runif(min = 0.9524, max = 0.9659, n = 1000)
T1S1 <- runif(min = 0.9608, max = 0.9677, n = 1000)

S0S0 <- runif(min=160.811, max=204.5009, n=1000)
S1S1 <- runif(min=168.989, max = 194.219, n=1000)
T0T0 <- runif(min=484.462, max = 616.082, n=1000)
T1T1 <- runif(min=514.279, max = 591.062, n=1000)

Mean_T0 <- runif(min=-13.455, max=-9.489, n=1000)
Mean_T1 <- runif(min=-17.17, max=-14.86, n=1000)
Mean_S0 <- runif(min=-7.789, max=-5.503, n=1000)
Mean_S1 <- runif(min=-9.600, max=-8.276, n=1000)

# Do the ISTE analysis
\dontrun{
ISTE <- ISTE.ContCont(Mean_T1=Mean_T1, Mean_T0=Mean_T0, 
 Mean_S1=Mean_S1, Mean_S0=Mean_S0, N=2128, Delta_S=c(-50:50), 
 zeta.PI=0.05, PI.Bound=0, Show.Prediction.Plots=TRUE,
 Save.Plots="No", T0S0=T0S0, T1S1=T1S1, T0T0=T0T0, T1T1=T1T1, 
 S0S0=S0S0, S1S1=S1S1)

# Examine results:
summary(ISTE)

# Plots of results. 
  # Plot ISTE
plot(ISTE)
  # Other plots, see plot.ISTE.ContCont for details
plot(ISTE, Outcome="MSE")
plot(ISTE, Outcome="gamma0")
plot(ISTE, Outcome="gamma1")
plot(ISTE, Outcome="Exp.DeltaT")
plot(ISTE, Outcome="Exp.DeltaT.Low.PI")
plot(ISTE, Outcome="Exp.DeltaT.Up.PI")
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{Individual-level surrogate threshold effect}
\keyword{ISTE}

