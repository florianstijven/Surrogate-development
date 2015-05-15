\name{plot Causal-Inference ContCont}
\alias{plot Causal-Inference ContCont}
\alias{plot.ICA.ContCont}
\alias{plot.MICA.ContCont}
\title{
Plots the (Meta-Analytic) Individual Causal Association when S and T are continuous outcomes
}
\description{
This function provides a plot that displays the frequencies, percentages, or cumulative percentages of the individual causal association (ICA; \eqn{\rho_{\Delta}}) and/or the meta-analytic individual causal association (MICA; \eqn{\rho_{M}}) values. These figures are useful to examine the sensitivity of the obtained results with respect to the assumptions regarding the correlations between the counterfactuals (for details, see Alonso et al., submitted; Van der Elst et al., submitted). Optionally, it is also possible to obtain plots that are useful in the examination of the plausibility of finding a good surrogate endpoint when an object of class \code{ICA.ContCont} is considered.}
\usage{
\method{plot}{ICA.ContCont}(x, Xlab.ICA, Main.ICA, Type="Percent", 
Labels=FALSE, ICA=TRUE, Good.Surr=FALSE, Main.Good.Surr, 
Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), col, \dots)

\method{plot}{MICA.ContCont}(x, ICA=TRUE, MICA=TRUE, Type="Percent", 
Labels=FALSE, Xlab.ICA, Main.ICA, Xlab.MICA, Main.MICA,
Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), col, \dots)
}
\arguments{
  \item{x}{An object of class \code{ICA.ContCont} or \code{MICA.ContCont}. See \code{\link{ICA.ContCont}} or \code{\link{MICA.ContCont}}.}
  \item{ICA}{Logical. When \code{ICA=TRUE}, a plot of the ICA is provided. Default \code{TRUE}.}
  \item{MICA}{Logical. This argument only has effect when the \code{plot()} function is applied to an object of class \code{MICA.ContCont}. When \code{MICA=TRUE}, a plot of the MICA is provided. Default \code{TRUE}.}
 \item{Type}{The type of plot that is produced. When \code{Type=Freq} or \code{Type=Percent}, the Y-axis shows frequencies or percentages of \eqn{\rho_{\Delta}}, \eqn{\rho_{M}}, and/or \eqn{\delta}. When \code{Type=CumPerc}, the Y-axis shows cumulative percentages of \eqn{\rho_{\Delta}}, \eqn{\rho_{M}}, and/or \eqn{\delta}.  Default "Percent".}  
 \item{Labels}{Logical. When \code{Labels=TRUE}, the percentage of \eqn{\rho_{\Delta}}, \eqn{\rho_{M}}, and/or \eqn{\delta} values that are equal to or larger than the midpoint value of each of the bins are displayed (on top of each bin). Default \code{FALSE}.}  
\item{Xlab.ICA}{The legend of the X-axis of the ICA plot. Default "\eqn{\rho_{\Delta}}".}
 \item{Main.ICA}{The title of the ICA plot. Default "ICA".}  
\item{Xlab.MICA}{The legend of the X-axis of the MICA plot. Default "\eqn{\rho_{M}}".}
 \item{Main.MICA}{The title of the MICA plot. Default "MICA".}    
\item{Good.Surr}{Logical. When \code{Good.Surr=TRUE}, a plot of \eqn{\delta} is provided. This plot is useful in the context of examinating the plausibility of finding a good surrogate endpoint. Only applies when an object of class \code{ICA.ContCont} is considered. For details, see Alonso et al. (submitted). Default \code{FALSE}.}
\item{Main.Good.Surr}{The title of the plot of \eqn{\delta}. Only applies when an object of class \code{ICA.ContCont} is considered. For details, see Alonso et al. (submitted).}
  \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
  \item{col}{The color of the bins. Default \code{col <- c(8)}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{hist()}.}
}
\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). On the relationship between the causal inference and meta-analytic paradigms for the validation of surrogate markers. 

Van der Elst, W., Alonso, A., & Molenberghs, G. (submitted). An exploration of the relationship between causal inference and meta-analytic measures of surrogacy.}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{ICA.ContCont}, \link{MICA.ContCont}, \link{plot MinSurrContCont}}
\examples{
# Plot of ICA

# Generate the vector of ICA values when rho_T0S0=rho_T1S1=.95, and when the
# grid of values {0, .2, ..., 1} is considered for the correlations
# between the counterfactuals:
SurICA <- ICA.ContCont(T0S0=.95, T1S1=.95, T0T1=seq(0, 1, by=.2), T0S1=seq(0, 1, by=.2), 
T1S0=seq(0, 1, by=.2), S0S1=seq(0, 1, by=.2))

# Plot the results:
plot(SurICA)

# Same plot but add the percentages of ICA values that are equal to or larger 
# than the midpoint values of the bins
plot(SurICA, Labels=TRUE)

# Plot of both ICA and MICA

# Generate the vector of ICA and MICA values when R_trial=.8, rho_T0S0=rho_T1S1=.8, 
# D.aa=5, D.bb=10, and when the grid of values {0, .2, ..., 1} is considered 
# for the correlations between the counterfactuals:
SurMICA <- MICA.ContCont(Trial.R=.80, D.aa=5, D.bb=10, T0S0=.8, T1S1=.8, 
T0T1=seq(0, 1, by=.2), T0S1=seq(0, 1, by=.2), T1S0=seq(0, 1, by=.2), 
S0S1=seq(0, 1, by=.2))

# Plot the vector of generated ICA and MICA values
plot(SurMICA, ICA=TRUE, MICA=TRUE)
}

\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Multiple-trial setting}
\keyword{Sensitivity}
\keyword{Plausibility of a surrogate}