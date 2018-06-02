\name{Fano.BinBin}
\alias{Fano.BinBin}
\title{Evaluate the possibility of finding a good surrogate in the setting where both \eqn{S} and \eqn{T} are binary endpoints}
\description{
The function \code{Fano.BinBin} evaluates the existence of a good surrogate in the single-trial causal-inference framework when both the surrogate and the true endpoints are binary outcomes. See \bold{Details} below.}
\usage{
Fano.BinBin(pi1_,  pi_1, rangepi10=c(0,min(pi1_,1-pi_1)), 
fano_delta=c(0.1), M=100, Seed=1)
}
\arguments{
  \item{pi1_}{A scalar or a vector of plausibel values that represents the proportion of responders under treatment.}
  \item{pi_1}{A scalar or a vector of plausibel values that represents the proportion of responders under control.}
  \item{rangepi10}{Represents the range from which \eqn{\pi_{10}} is sampled. By default, Monte Carlo simulation will be constrained to the interval \eqn{[0,\min(\pi_{1\cdot},\pi_{\cdot0})]} but this allows the user to specify a more narrow range. \code{rangepi10=c(0,0)} is equivalent to the assumption of monotonicity for the true endpoint.}
  \item{fano_delta}{A scalar or a vector that specifies the values for the upper bound of the prediction error \eqn{\delta}. Default \code{fano_delta=c(0.2)}.}
  \item{M}{The number of random samples that have to be drawn for the freely varying parameter \eqn{\pi_{10}}. Default \code{M=1000}. The number of random samples should be sufficiently large in relation to the length of the interval \code{rangepi10}. Typically \code{M=1000} yields a sufficiently fine grid. In case \code{rangepi10} is a single value: \code{M=1} }
  \item{Seed}{The seed to be used to sample the freely varying parameter \eqn{\pi_{10}}. Default \code{Seed=1}.}

}

\details{
Values for \eqn{\pi_{10}} have to be uniformly sampled from the interval \eqn{[0,\min(\pi_{1\cdot},\pi_{\cdot0})]}. Any sampled value for \eqn{\pi_{10}} will fully determine the bivariate distribution of potential outcomes for the true endpoint. The treatment effect should be positive.

The vector \eqn{\bold{\pi_{km}}} fully determines \eqn{R^2_{HL}}.
}

\value{An object of class \code{Fano.BinBin} with components,
\item{R2_HL}{The sampled values for \eqn{R^2_{HL}}.}
\item{H_Delta_T}{The sampled values for \eqn{H{\Delta T}}.}
\item{PPE_T}{The sampled values for \eqn{PPE_{T}}.}
\item{minpi10}{The minimum value for \eqn{\pi_{10}}.}
\item{maxpi10}{The maximum value for \eqn{\pi_{10}}.}
\item{samplepi10}{The sampled value for \eqn{\pi_{10}}.}
\item{delta}{The specified vector of upper bounds for the prediction errors.}
\item{uncertainty}{Indexes the sampling of \eqn{pi1\_}.}
\item{pi_00}{The sampled values for \eqn{\pi_{00}}.}
\item{pi_11}{The sampled values for \eqn{\pi_{11}}.}
\item{pi_01}{The sampled values for \eqn{\pi_{01}}.}
\item{pi_10}{The sampled values for \eqn{\pi_{10}}.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2014). Validation of surrogate endpoints: the binary-binary setting from a causal inference perspective.
}
\author{
Paul Meyvisch, Wim Van der Elst, Ariel Alonso
}
\seealso{
\code{\link{plot.Fano.BinBin}}
}
\examples{
# Conduct the analysis assuming no montonicity
# for the true endpoint, using a range of
# upper bounds for prediction errors 
Fano.BinBin(pi1_ = 0.5951 ,  pi_1 = 0.7745, 
fano_delta=c(0.05, 0.1, 0.2), M=1000)


# Conduct the same analysis now sampling from
# a range of values to allow for uncertainty

Fano.BinBin(pi1_ = runif(n=20,min=0.504,max=0.681), 
pi_1 = runif(n=20,min=0.679,max=0.849), 
fano_delta=c(0.05, 0.1, 0.2), M=10, Seed=2)
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Fano}
\keyword{ICA}
\keyword{MarginalProbs}