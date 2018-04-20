\name{Bootstrap.MEP.BinBin}
\alias{Bootstrap.MEP.BinBin}
\title{
Bootstrap 95\% CI around the maximum-entropy ICA and SPF (surrogate predictive function)}
\description{
 Computes a 95\% bootstrap-based CI around the maximum-entropy ICA and SPF (surrogate predictive function) in the binary-binary setting}
\usage{
Bootstrap.MEP.BinBin(Data, Surr, True, Treat, M=100, Seed=123)
}
\arguments{
  \item{Data}{The dataset to be used.}
  \item{Surr}{The name of the surrogate variable.}
  \item{True}{The name of the true endpoint.}
  \item{Treat}{The name of the treatment indicator.}
  \item{M}{The number of bootstrap samples taken. Default \code{M=1000}.}
  \item{Seed}{The seed to be used. Default \code{Seed=123}.}
}
\value{
\item{R2H}{The vector the bootstrapped MEP ICA values.}
\item{r_1_1}{The vector of the bootstrapped bootstrapped MEP \eqn{r(1, 1)} values.}
\item{r_min1_1}{The vector of the bootstrapped bootstrapped MEP \eqn{r(-1, 1)}.}
\item{r_0_1}{The vector of the bootstrapped bootstrapped MEP \eqn{r(0, 1)}.}
\item{r_1_0}{The vector of the bootstrapped bootstrapped MEP \eqn{r(1, 0)}.}
\item{r_min1_0}{The vector of the bootstrapped bootstrapped MEP \eqn{r(-1, 0)}.}
\item{r_0_0}{The vector of the bootstrapped bootstrapped MEP \eqn{r(0, 0)}.}
\item{r_1_min1}{The vector of the bootstrapped bootstrapped MEP \eqn{r(1, -1)}.}
\item{r_min1_min1}{The vector of the bootstrapped bootstrapped MEP \eqn{r(-1, -1)}.}
\item{r_0_min1}{The vector of the bootstrapped bootstrapped MEP \eqn{r(0, -1)}.}
\item{vector_p}{The matrix that contains all bootstrapped maximum entropy distributions of the vector of the potential outcomes.}
}
\references{
Alonso, A., & Van der Elst, W. (2015). A maximum-entropy approach for the evluation of surrogate endpoints based on causal inference.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinBin}}, \code{\link{ICA.BinBin.Grid.Sample}}, \code{\link{ICA.BinBin.Grid.Full}}, \code{\link{plot MaxEntSPF BinBin}}
}
\examples{
\dontrun{ # time consuming code part
MEP_CI <- Bootstrap.MEP.BinBin(Data = Schizo_Bin, Surr = "BPRS_Bin", True = "PANSS_Bin",
                     Treat = "Treat", M = 500, Seed=123)
summary(MEP_CI)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{BinBin}
\keyword{Maximum Entropy}
\keyword{Bootstrap}