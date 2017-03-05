\name{plot ICA.ContCont.MultS}
\alias{plot.ICA.ContCont.MultS}
\alias{plot.ICA.ContCont.MultS_alt}
\title{
Plots the Individual Causal Association in the setting where there are multiple continuous S and a continuous T}
\description{
This function provides a plot that displays the frequencies, percentages, or cumulative percentages of the multivariate individual causal association (\eqn{R^2_{H}}). These figures are useful to examine the sensitivity of the obtained results with respect to the assumptions regarding the correlations between the counterfactuals.}
\usage{
  \method{plot}{ICA.ContCont.MultS}(x, R2_H=FALSE, Corr.R2_H=TRUE, 
   Type="Percent", Labels=FALSE,  
   Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), col, 
   Prediction.Error.Reduction=FALSE, \dots)
}
\arguments{
  \item{x}{An object of class \code{ICA.ContCont.MultS}. See \code{\link{ICA.ContCont.MultS}} or \code{\link{ICA.ContCont.MultS_alt}}.}
\item{R2_H}{Should a plot of the \eqn{R^2_{H}} be provided? Default \code{FALSE}.}
\item{Corr.R2_H}{Should a plot of the corrected \eqn{R^2_{H}} be provided? Default \code{TRUE}.}
  \item{Type}{The type of plot that is produced. When \code{Type=Freq} or \code{Type=Percent}, the Y-axis shows frequencies or percentages of \eqn{R^2_{H}}. When \code{Type=CumPerc}, the Y-axis shows cumulative percentages of \eqn{R^2_{H}}.  Default "Percent".}  
  \item{Labels}{Logical. When \code{Labels=TRUE}, the percentage of \eqn{R^2_{H}} values that are equal to or larger than the midpoint value of each of the bins are displayed (on top of each bin). Default \code{FALSE}.}  
  \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
  \item{col}{The color of the bins. Default \code{col <- c(8)}.}
  \item{Prediction.Error.Reduction}{Should a plot be shown that shows the prediction error (reisdual error) in predicting \eqn{Delta}T using an intercept only model, and that shows the prediction error (reisdual error) in predicting \eqn{Delta}T using \eqn{Delta S_1}, \eqn{Delta S_2}, ...? Default \code{Prediction.Error.Reduction=FALSE}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{hist()}.}
}
\references{
Van der Elst, W., Alonso, A. A., & Molenberghs, G. (2017). Univariate versus multivariate surrogate endpoints.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{ICA.ContCont}, \link{ICA.ContCont.MultS}, \link{ICA.ContCont.MultS_alt}, \link{MICA.ContCont}, \link{plot MinSurrContCont}}
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

# Marix looks like:
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
  Sigma=s, G = seq(from=-1, to=1, by = .00001), Seed=c(123), 
  Model = "Delta_T ~ Delta_S1 + Delta_S2 + Delta_S3")

# Explore results
summary(ICA)
plot(ICA)
}
}
\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Multiple-trial setting}
\keyword{Sensitivity}
\keyword{Multiple surrogates}
\keyword{Multivariate setting}