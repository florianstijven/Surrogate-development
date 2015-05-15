\name{Test.Mono}
\alias{Test.Mono}
\title{
Test whether the data are compatible with monotonicity for S and/or T (binary endpoints)
}
\description{
For some situations, the observable marginal probabilities contain sufficient information to exclude a particular monotonicity scenario. For example, under monotonicity for \eqn{S} and \eqn{T}, one of the restrictions that the data impose is \eqn{\pi_{0111}<min(\pi_{0 \cdot 1 \cdot}, \pi_{\cdot 1 \cdot 1})}. If the latter condition does not hold in the dataset at hand, monotonicity for \eqn{S} and \eqn{T} can be excluded. 
}
\usage{
Test.Mono(pi1_1_, pi0_1_, pi1_0_, pi_1_1, pi_1_0, pi_0_1)
}
\arguments{
  \item{pi1_1_}{A scalar that contains \eqn{P(T=1,S=1|Z=0)}.}
  \item{pi0_1_}{A scalar that contains \eqn{P(T=0,S=1|Z=0)}.}
  \item{pi1_0_}{A scalar that contains \eqn{P(T=1,S=0|Z=0)}.}
  \item{pi_1_1}{A scalar that contains \eqn{P(T=1,S=1|Z=1)}.}
  \item{pi_1_0}{A scalar that contains \eqn{P(T=1,S=0|Z=1)}.}
  \item{pi_0_1}{A scalar that contains \eqn{P(T=0,S=1|Z=1)}.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2015). Validation of surrogate endpoints: the binary-binary setting from a causal inference perspective.
}
\author{
Wim Van der Elst, Ariel Alonso, Marc Buyse, & Geert Molenberghs
}
\examples{
Test.Mono(pi1_1_=0.2619048, pi1_0_=0.2857143, pi_1_1=0.6372549, 
pi_1_0=0.07843137, pi0_1_=0.1349206, pi_0_1=0.127451)
}
\keyword{Monotonicity}
\keyword{Test Monotonicity}