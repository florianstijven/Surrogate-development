\name{MaxEntSPFBinBin}
\alias{MaxEntSPFBinBin}
\title{
Use the maximum-entropy approach to compute SPF (surrogate predictive function) in the binary-binary setting
}
\description{
In a surrogate evaluation setting where both \eqn{S} and \eqn{T} are binary
endpoints, a sensitivity-based approach where multiple 'plausible values' for vector \eqn{\pi} (i.e., vectors \eqn{\pi} that are compatible with the observable data at hand) can be used (for details, see \code{\link{SPF.BinBin}}). Alternatively, the maximum entropy distribution for vector \eqn{\pi}
can be considered (Alonso et al., 2015). The use of the distribution that maximizes the entropy can be justified
based on the fact that any other distribution would necessarily 
(i) assume information that we do not have, or (ii) contradict information
that we do have. The function \code{MaxEntSPFBinBin} implements the latter approach.

Based on vector \eqn{\pi}, the surrogate predictive function (SPF) is computed, i.e., \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)}. For example, \eqn{r(-1,1)} quantifies the probability that the treatment has a negative effect on the true endpoint (\eqn{\Delta T=-1}) given that it has a positive effect on the surrogate (\eqn{\Delta S=1}). }
\usage{
MaxEntSPFBinBin(pi1_1_, pi1_0_, pi_1_1,
pi_1_0, pi0_1_, pi_0_1, Method="BFGS", 
Fitted.ICA=NULL)
}
\arguments{
  \item{pi1_1_}{A scalar that contains the estimated value for \eqn{P(T=1,S=1|Z=0)}, i.e., the probability that \eqn{S=T=1} when under treatment \eqn{Z=0}.}
  \item{pi1_0_}{A scalar that contains the estimated value for \eqn{P(T=1,S=0|Z=0)}.}
  \item{pi_1_1}{A scalar that contains the estimated value for \eqn{P(T=1,S=1|Z=1)}.}
  \item{pi_1_0}{A scalar that contains the estimated value for \eqn{P(T=1,S=0|Z=1)}.}
  \item{pi0_1_}{A scalar that contains the estimated value for \eqn{P(T=0,S=1|Z=0)}.}
  \item{pi_0_1}{A scalar that contains the estimated value for \eqn{P(T=0,S=1|Z=1)}.}
  \item{Method}{The maximum entropy frequency vector \eqn{p^{*}} is calculated based on the optimal solution to an unconstrained dual convex programming problem (for details, see Alonso et al., 2015). Two different optimization methods can be specified, i.e., \code{Method="BFGS"} and \code{Method="CG"}, which implement the quasi-Newton BFGS (Broyden, Fletcher, Goldfarb, and Shanno) and the conjugent gradient (CG) methods (for details on these methods, see the help files of the \code{optim()} function and the references theirin). 
  Alternatively, the \eqn{\pi} vector (obtained when the functions \code{ICA.BinBin}, \code{ICA.BinBin.Grid.Full}, or \code{ICA.BinBin.Grid.Sample} are executed) that is 'closest' to the vector \eqn{\pi} can be retained. Here, the 'closest' vector is defined as the vector where the sum of the squared differences between
the components in the vectors \eqn{\pi} and \eqn{\pi} is smallest. The latter 'Minimum Difference' method can re requested by specifying the argument \code{Method="MD"} in the function call. Default \code{Method="BFGS"}.}
  \item{Fitted.ICA}{A fitted object of class \code{ICA.BinBin}, \code{ICA.BinBin.Grid.Full}, or \code{ICA.BinBin.Grid.Sample}. Only required when \code{Method="MD"} is used. } 
}

\value{
\item{Vector_p}{The maximum entropy frequency vector \eqn{p^{*}}}
\item{r_1_1}{The vector of values for \eqn{r(1, 1)}, i.e., \eqn{P(\Delta T=1|\Delta S=1}).}
\item{r_min1_1}{The vector of values for \eqn{r(-1, 1)}.}
\item{r_0_1}{The vector of values for \eqn{r(0, 1)}.}
\item{r_1_0}{The vector of values for \eqn{r(1, 0)}.}
\item{r_min1_0}{The vector of values for \eqn{r(-1, 0)}.}
\item{r_0_0}{The vector of values for \eqn{r(0, 0)}.}
\item{r_1_min1}{The vector of values for \eqn{r(1, -1)}.}
\item{r_min1_min1}{The vector of values for \eqn{r(-1, -1)}.}
\item{r_0_min1}{The vector of values for \eqn{r(0, -1)}.}
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
# Sensitivity-based ICA results using ICA.BinBin.Grid.Sample
ICA <- ICA.BinBin.Grid.Sample(pi1_1_=0.341, pi0_1_=0.119, pi1_0_=0.254,
pi_1_1=0.686, pi_1_0=0.088, pi_0_1=0.078, Seed=1, 
Monotonicity=c("No"), M=5000)

# Sensitivity-based SPF
SPFSens <- SPF.BinBin(ICA)

# Maximum-entropy based SPF
SPFMaxEnt <- MaxEntSPFBinBin(pi1_1_=0.341, pi0_1_=0.119, pi1_0_=0.254,
pi_1_1=0.686, pi_1_0=0.088, pi_0_1=0.078)

# Explore maximum-entropy results
summary(SPFMaxEnt)

# Plot results
plot(x=SPFMaxEnt, SPF.Fit=SPFSens)
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{BinBin}
\keyword{Maximum Entropy}