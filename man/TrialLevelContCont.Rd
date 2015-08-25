\name{TrialLevelContCont}
\alias{TrialLevelContCont}
\title{Estimates trial-level surrogacy based on treatment effects on S and T (continuous-continuous case)
}
\description{
The function \code{TrialLevelContCont} estimates trial-level surrogacy based on the vectors of treatment effects on \eqn{S} (i.e., \eqn{\alpha_{i}}) and \eqn{T} (i.e., \eqn{\beta_{i}}) in the different trials. In particular, \eqn{\beta_{i}} is regressed on \eqn{\alpha_{i}} and the classical coefficient of determination of the fitted model provides an estimate of \eqn{R^2_{trial}}. In addition, the standard error and CI are provided.   
}
\usage{
TrialLevelContCont(Alpha.Vector, Beta.Vector, 
N.Vector, Weighted=TRUE, Alpha=.05)
}
\arguments{
\item{Alpha.Vector}{The vector of treatment effects on \eqn{S} in the different trials, i.e., \eqn{\alpha_{i}}.}
\item{Beta.Vector}{The vector of treatment effects on \eqn{T} in the different trials, i.e., \eqn{\beta_{i}}.}
\item{N.Vector}{The vector of trial sizes \eqn{N_{i}}.}
\item{Weighted}{Logical. If \code{TRUE}, then a weighted regression analysis is conducted. If \code{FALSE}, then an unweighted regression analysis is conducted. Default \code{TRUE}.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}} and \eqn{R_{trial}}. Default \eqn{0.05}.}
}
\value{An object of class \code{TrialLevelContCont} with components,
\item{Alpha.Vector}{The vector of treatment effects on \eqn{S} in the different trials.}
\item{Beta.Vector}{The vector of treatment effects on \eqn{T} in the different trials.}
\item{N.Vector}{The vector of trial sizes \eqn{N_{i}}.}
\item{Trial.R2}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{trial}}), its standard error and confidence interval.}
\item{Trial.R}{A \code{data.frame} that contains the trial-level correlation coefficient (\eqn{R_{trial}}), its standard error and confidence interval.}
}
\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{UnimixedContCont}}, \code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{BimixedContCont}}, \code{\link{plot Meta-Analytic}}
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
Fit <- TrialLevelContCont(Alpha.Vector=Alpha.Vector,
Beta.Vector=Beta.Vector, N.Vector=N.Vector)

# Plot the results and obtain summary
plot(Fit)
summary(Fit)
}

\keyword{Multiple-trial setting}
\keyword{Meta-analytic framework}
\keyword{Trial-level surrogacy}