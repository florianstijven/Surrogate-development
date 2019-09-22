\name{plot ISTE.ContCont}
\alias{plot ISTE.ContCont}
\alias{plot.ISTE.ContCont}
\title{
Plots the individual-level surrogate threshold effect (STE) values and related metrics
}
\description{
This function plots the individual-level surrogate threshold effect (STE) values and related metrics, e.g., the expected \eqn{\Delta T} values for a vector of \eqn{\Delta S} values.}
\usage{
\method{plot}{ISTE.ContCont}(x, Outcome="ISTE", breaks=50, \dots)
}
\arguments{
\item{x}{An object of class \code{ISTE.ContCont}. See \code{\link{ISTE.ContCont}}.}
\item{Outcome}{The outcome for which a histogram has to be produced. When \code{Outcome="ISTE"}, 
a histogram of the ISTE is produced. When \code{Outcome="MSE"}, a histogram of the MSE values (of regression models in which \eqn{\Delta T} is regressed on \eqn{\Delta S}) is given. When \code{Outcome="gamma0"}, a histogram of \eqn{\gamma[0]} values (of regression models in which \eqn{\Delta T} is regressed on \eqn{\Delta S}) is given. When \code{Outcome="gamma1"}, a histogram of \eqn{\gamma[1]} values (of regression models in which \eqn{\Delta T} is regressed on \eqn{\Delta S}) is given. When \code{Outcome="Exp.DeltaT"}, a histogram of the expected \eqn{\Delta T} values for a vector of \eqn{\Delta S} values (specified in the call of the \code{ISTE.ContCont} function) values is given. When \code{Outcome="Exp.DeltaT.Low.PI"}, a histogram of the lower prediction intervals of the expected \eqn{\Delta T} values for a vector of \eqn{\Delta S} values (specified in the call of the \code{ISTE.ContCont} function) values is given. When \code{Outcome="Exp.DeltaT.Up.PI"}, a histogram of the upper prediction intervals of the expected \eqn{\Delta T} values for a vector of \eqn{\Delta S} values (specified in the call of the \code{ISTE.ContCont} function) values is given. Dafault \code{Outcome="ISTE"}. When \code{Outcome="Delta_S_For_Which_Delta_T_equal_0"}, a histogram of \eqn{omega} is shown with \eqn{E(Delta T | Delta S > omega)>0}. 
}
  \item{breaks}{The number of breaks used in the histogram(s). Default \code{breaks=50}.} 
 \item{\dots}{Extra graphical parameters to be passed to \code{hist()}.}
}
\references{
Van der Elst, W., Alonso, A. A., and Molenberghs, G. (submitted). The individual-level surrogate threshold effect in a causal-inference setting.}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{ISTE.ContCont}}
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
 alpha.PI=0.05, PI.Bound=0, Show.Prediction.Plots=TRUE,
 Save.Plots="No", T0S0=T0S0, T1S1=T1S1, T0T0=T0T0, T1T1=T1T1, 
 S0S0=S0S0, S1S1=S1S1)

# Examine results:
summary(ISTE)

# Plots of results. 
  # Plot main ISTE results
plot(ISTE)
  # Other plots
plot(ISTE, Outcome="MSE")
plot(ISTE, Outcome="gamma0")
plot(ISTE, Outcome="gamma1")
plot(ISTE, Outcome="Exp.DeltaT")
plot(ISTE, Outcome="Exp.DeltaT.Low.PI")
plot(ISTE, Outcome="Exp.DeltaT.Up.PI")
}
}
\keyword{Plot ISTE}
\keyword{Plot individual-level surrogate threshold effect}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Sensitivity}