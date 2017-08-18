\name{ECT}
\alias{ECT}
\title{Apply the Entropy Concentration Theorem}
\description{
The Entropy Concentration Theorem (ECT; Edwin, 1982) states that if \eqn{N} is large enough, then \eqn{100(1-F)}\% of all \eqn{\bold{p*}} and \eqn{\Delta H} is determined by the upper tail are \eqn{1-F} of a \eqn{\chi^2} distribution, with \eqn{DF = q - m - 1} (which equals \eqn{8} in a surrogate evaluation context).
}
\usage{
ECT(Perc=.95, H_Max, N)
}
\arguments{
  \item{Perc}{The desired interval. E.g., \code{Perc=.05} will generate the lower and upper bounds for \eqn{H(\bold{p})} that contain \eqn{95\%} of the cases (as determined by the ECT).}
  \item{H_Max}{The maximum entropy value. In the binary-binary setting, this can be computed using the function \code{\link{MaxEntICABinBin}}.}
  \item{N}{The sample size.}
  }
\value{An object of class \code{ECT} with components,
\item{Lower_H}{The lower bound of the requested interval.}
\item{Upper_H}{The upper bound of the requested interval, which equals \eqn{H_Max}.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2016). Surrogate markers validation: the continuous-binary setting from a causal inference perspective.
}
\author{
Wim Van der Elst, Paul Meyvisch, & Ariel Alonso
}
\seealso{
\code{\link{MaxEntICABinBin}}, \code{\link{ICA.BinBin}}
}
\examples{
ECT_fit <- ECT(Perc = .05, H_Max = 1.981811, N=454)
summary(ECT_fit)
}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Binary endpoint}
\keyword{Maximum entropy}
\keyword{Entropy Concentration Theorem}