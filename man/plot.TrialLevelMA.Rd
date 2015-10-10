\name{plot TrialLevelMA}
\alias{plot.TrialLevelMA}
\title{Provides a plots of trial-level surrogacy in the meta-analytic framework based on the output of the \code{TrialLevelMA()} function
}
\description{
Produces a plot that provides a graphical representation of trial-level surrogacy based on the output of the \code{TrialLevel()} function (meta-analytic framework).   
}
\usage{
\method{plot}{TrialLevelMA}(x, Weighted=TRUE, Xlab.Trial, 
Ylab.Trial, Main.Trial, Par=par(oma=c(0, 0, 0, 0), 
mar=c(5.1, 4.1, 4.1, 2.1)), \dots)
}
\arguments{
  \item{x}{An object of class \code{TrialLevelMA}.}
  \item{Weighted}{Logical. If \code{Weighted=TRUE}, the circles that depict the trial-specific treatment effects on the true endpoint against the surrogate endpoint are proportional to the number of patients in the trial. If \code{Weighted=FALSE}, all circles have the same size. Default \code{TRUE}.}
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
\link{UnifixedContCont}, \link{BifixedContCont}, \link{UnifixedContCont}, \link{BimixedContCont}, \link{TrialLevelMA}
}
\examples{
# Generate vector treatment effects on S
set.seed(seed = 1)
Alpha.Vector <- seq(from = 5, to = 10, by=.1) + runif(min = -.5, max = .5, n = 51)
# Generate vector treatment effects on T
set.seed(seed=2)
Beta.Vector <- (Alpha.Vector * 3) + runif(min = -5, max = 5, n = 51)
# Vector of sample sizes of the trials (here, all n_i=10)
N.Vector <- rep(10, times=51)

# Apply the function to estimate R^2_{trial}
Fit <- TrialLevelMA(Alpha.Vector=Alpha.Vector,
Beta.Vector=Beta.Vector, N.Vector=N.Vector)

# Plot the results and obtain summary
plot(Fit)
summary(Fit)
}

\keyword{Plot surrogacy}
\keyword{Meta-analytic framework}
\keyword{Trial-level surrogacy}
\keyword{Multiple-trial setting}