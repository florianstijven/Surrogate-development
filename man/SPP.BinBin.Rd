\name{SPF.BinBin}
\alias{SPF.BinBin}
\title{
Evaluate the surrogate predictive function (SPF) in the binary-binary setting
}
\description{
Computes the surrogate predictive function (SPF), i.e., \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)}, in the setting where both \eqn{S} and \eqn{T} are binary endpoints. For example, \eqn{r(-1,1)} quantifies the probability that the treatment has a negative effect on the true endpoint (\eqn{\Delta T=-1}) given that it has a positive effect on the surrogate (\eqn{\Delta S=1}).
}
\usage{
SPF.BinBin(x)
}
\arguments{
  \item{x}{A fitted object of class \code{ICA.BinBin}, \code{ICA.BinBin.Grid.Full}, or \code{ICA.BinBin.Grid.Sample}.} 
}

\value{
\item{r_1_1}{The vector of values for \eqn{r(1, 1)}, i.e., \eqn{P(\Delta T=1|\Delta S=1}).}
\item{r_min1_1}{The vector of values for \eqn{r(-1, 1)}.}
\item{r_0_1}{The vector of values for \eqn{r(0, 1)}.}
\item{r_1_0}{The vector of values for \eqn{r(1, 0)}.}
\item{r_min1_0}{The vector of values for \eqn{r(-1, 0)}.}
\item{r_0_0}{The vector of values for \eqn{r(0, 0)}.}
\item{r_1_min1}{The vector of values for \eqn{r(1, -1)}.}
\item{r_min1_min1}{The vector of values for \eqn{r(-1, -1)}.}
\item{r_0_min1}{The vector of values for \eqn{r(0, -1)}.}
\item{Monotonicity}{The assumption regarding monotonicity under which the result was obtained.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2015). Assessing a surrogate effect predictive value in a causal inference framework.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinBin}}, \code{\link{ICA.BinBin.Grid.Sample}}, \code{\link{plot.SPF.BinBin}}
}
\examples{
\dontrun{
# Use ICA.BinBin.Grid.Sample to obtain plausible values for pi
ICA_BINBIN_Grid_Sample <- ICA.BinBin.Grid.Sample(pi1_1_=0.341, pi0_1_=0.119,
pi1_0_=0.254, pi_1_1=0.686, pi_1_0=0.088, pi_0_1=0.078, Seed=1,
Monotonicity=c("General"), M=5000)

# Obtain SPF
SPF.BinBin(ICA_BINBIN_Grid_Sample)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{BinBin}
\keyword{Sensitivity}
\keyword{SPF}