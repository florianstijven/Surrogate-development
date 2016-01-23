\name{SPF.BinBin}
\alias{SPF.BinBin}
\title{
Evaluate the surrogate predictive function (SPF) in the binary-binary setting (sensitivity-analysis based approach)
}
\description{
Computes the surrogate predictive function (SPF) based on sensitivity-analyis, i.e., \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)}, in the setting where both \eqn{S} and \eqn{T} are binary endpoints. For example, \eqn{r(-1,1)} quantifies the probability that the treatment has a negative effect on the true endpoint (\eqn{\Delta T=-1}) given that it has a positive effect on the surrogate (\eqn{\Delta S=1}). All quantities of interest are derived from the vectors of 'plausible values' for \eqn{\pi} (i.e., vectors \eqn{\pi} that are compatible with the observable data at hand). See \bold{Details} below.
}
\usage{
SPF.BinBin(x)
}
\arguments{
  \item{x}{A fitted object of class \code{ICA.BinBin}, \code{ICA.BinBin.Grid.Full}, or \code{ICA.BinBin.Grid.Sample}.} 
}
\details{

All \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} are derived from \eqn{\pi} (vector of potential outcomes). Denote by \eqn{\bold{Y}'=(T_0,T_1,S_0,S_1)} the vector of potential outcomes. The vector \eqn{\bold{Y}} can take 16 values and the set of parameters \eqn{\pi_{ijpq}=P(T_0=i,T_1=j,S_0=p,S_1=q)} (with \eqn{i,j,p,q=0/1}) fully characterizes its distribution.

Based on the data and assuming SUTVA, the marginal probabilites \eqn{\pi_{1 \cdot 1 \cdot}}, \eqn{\pi_{1 \cdot 0 \cdot}}, \eqn{\pi_{\cdot 1 \cdot 1}}, \eqn{\pi_{\cdot 1 \cdot 0}}, \eqn{\pi_{0 \cdot 1 \cdot}}, and \eqn{\pi_{\cdot 0 \cdot 1}} can be computed (by hand or using the function \code{\link{MarginalProbs}}). Define the vector 
\deqn{\bold{b}'=(1, \pi_{1 \cdot 1 \cdot}, \pi_{1 \cdot 0 \cdot}, \pi_{\cdot 1 \cdot 1}, \pi_{\cdot 1 \cdot 0}, \pi_{0 \cdot 1 \cdot}, \pi_{\cdot 0 \cdot 1})} and \eqn{\bold{A}} is a contrast matrix such that the identified restrictions can be written as a system of linear equation
\deqn{\bold{A \pi} = \bold{b}.}

The matrix \eqn{\bold{A}} has rank \eqn{7} and can be partitioned as \eqn{\bold{A=(A_r | A_f)}}, and similarly the vector \eqn{\bold{\pi}} can be partitioned as \eqn{\bold{\pi^{'}=(\pi_r^{'} | \pi_f^{'})}} (where \eqn{f} refers to the submatrix/vector given by the \eqn{9} last columns/components of \eqn{\bold{A/\pi}}). Using these partitions the previous system of linear equations can be rewritten as \deqn{\bold{A_r \pi_r + A_f \pi_f = b}.}

The functions \code{\link{ICA.BinBin}}, \code{\link{ICA.BinBin.Grid.Sample}}, and  \code{\link{ICA.BinBin.Grid.Full}} contain algorithms that generate plausible distributions for \eqn{\bold{Y}} (for details, see the documentation of these functions). Based on the output of these functions, \code{SPF.BinBin} computes the surrogate predictive function.
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
Wim Van der Elst, Paul Meyvisch, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinBin}}, \code{\link{ICA.BinBin.Grid.Sample}}, \code{\link{ICA.BinBin.Grid.Full}}, \code{\link{plot.SPF.BinBin}}
}
\examples{
# Use ICA.BinBin.Grid.Sample to obtain plausible values for pi
ICA_BINBIN_Grid_Sample <- ICA.BinBin.Grid.Sample(pi1_1_=0.341, pi0_1_=0.119,
pi1_0_=0.254, pi_1_1=0.686, pi_1_0=0.088, pi_0_1=0.078, Seed=1,
Monotonicity=c("General"), M=2500)

# Obtain SPF
SPF <- SPF.BinBin(ICA_BINBIN_Grid_Sample)

# examine results
summary(SPF)
plot(SPF)
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{BinBin}
\keyword{Sensitivity}
\keyword{SPF}