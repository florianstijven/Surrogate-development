\name{plot TrialLevelIT}
\alias{plot.TrialLevelIT}
\title{Provides a plots of trial-level surrogacy in the information-theoretic framework based on the output of the \code{TrialLevelIT()} function
}
\description{
Produces a plot that provides a graphical representation of trial-level surrogacy based on the output of the \code{TrialLevelIT()} function (information-theoretic framework).   
}
\usage{
\method{plot}{TrialLevelIT}(x, Xlab.Trial, 
Ylab.Trial, Main.Trial, Par=par(oma=c(0, 0, 0, 0), 
mar=c(5.1, 4.1, 4.1, 2.1)), \dots)
}
\arguments{
  \item{x}{An object of class \code{TrialLevelIT}.}
  \item{Xlab.Trial}{The legend of the X-axis of the plot that depicts trial-level surrogacy. Default "Treatment effect on the surrogate endpoint (\eqn{\alpha_{i}})".}
  \item{Ylab.Trial}{The legend of the Y-axis of the plot that depicts trial-level surrogacy. Default "Treatment effect on the true endpoint (\eqn{\beta_{i}})".}
  \item{Main.Trial}{The title of the plot that depicts trial-level surrogacy. Default "Trial-level surrogacy".}
  \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{plot()}.}
}
\references{
Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\link{UnifixedContCont}, \link{BifixedContCont}, \link{UnifixedContCont}, \link{BimixedContCont}, \link{TrialLevelIT}
}
\examples{
# Generate vector treatment effects on S
set.seed(seed = 1)
Alpha.Vector <- seq(from = 5, to = 10, by=.1) + runif(min = -.5, max = .5, n = 51)

# Generate vector treatment effects on T
set.seed(seed=2)
Beta.Vector <- (Alpha.Vector * 3) + runif(min = -5, max = 5, n = 51)

# Apply the function to estimate R^2_{h.t}
Fit <- TrialLevelIT(Alpha.Vector=Alpha.Vector,
Beta.Vector=Beta.Vector, N.Trial=50, Model="Reduced")

# Plot the results
plot(Fit)
}
\keyword{Plot surrogacy}
\keyword{Information-theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Multiple-trial setting}