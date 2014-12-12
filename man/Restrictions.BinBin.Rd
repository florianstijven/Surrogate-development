\name{Restrictions.BinBin}
\alias{Restrictions.BinBin}
\title{Examine restrictions in \eqn{\bold{\pi}_{f}} under different montonicity assumptions for binary \eqn{S} and \eqn{T}}
\description{
The function \code{Restrictions.BinBin} gives an overview of the restrictions in \eqn{\bold{\pi}_{f}}
under different assumptions regarding montonicity when both \eqn{S} and \eqn{T} are binary.}

\usage{
Restrictions.BinBin(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1)
}

\arguments{
  \item{pi1_1_}{A scalar that contains \eqn{P(T=1,S=1|Z=0)}, i.e., the proability that \eqn{S=T=1} when under treatment \eqn{Z=0}.}
  \item{pi1_0_}{A scalar that contains \eqn{P(T=1,S=0|Z=0)}.}
  \item{pi_1_1}{A scalar that contains \eqn{P(T=1,S=1|Z=1)}.}
  \item{pi_1_0}{A scalar that contains \eqn{P(T=1,S=0|Z=1)}.}
  \item{pi0_1_}{A scalar that contains \eqn{P(T=0,S=1|Z=0)}.}
  \item{pi_0_1}{A scalar that contains \eqn{P(T=0,S=1|Z=1)}.}
  }

\value{An overview of the restrictions for the freely varying parameters imposed by the data is provided}

\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2014). Validation of surrogate endpoints: the binary-binary setting from a causal inference perspective.
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{
\code{\link{MarginalProbs}}
}

\examples{
Restrictions.BinBin(pi1_1_=0.262, pi0_1_=0.135, pi1_0_=0.286, 
pi_1_1=0.637, pi_1_0=0.078, pi_0_1=0.127)
}

\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}