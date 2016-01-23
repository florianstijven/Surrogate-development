\name{plot MaxEntSPF BinBin}
\alias{plot MaxEntSPF BinBin}
\alias{plot.MaxEntSPF.BinBin}
\title{
Plots the sensitivity-based and maximum entropy based surrogate predictive function (SPF) when S and T are binary outcomes.
}
\description{
Plots the sensitivity-based (Alonso et al., 2015a) and maximum entropy based (Alonso et al., 2015b) surrogate predictive function (SPF), i.e., \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)}, in the setting where both \eqn{S} and \eqn{T} are binary endpoints. For example, \eqn{r(-1,1)} quantifies the probability that the treatment has a negative effect on the true endpoint (\eqn{\Delta T=-1}) given that it has a positive effect on the surrogate (\eqn{\Delta S=1}).}

\usage{
\method{plot}{MaxEntSPF.BinBin}(x, SPF.Fit, Type="All.Histograms", Col="grey", ...)
}
\arguments{
\item{x}{A fitted object of class \code{MaxEntSPF.BinBin}. See \code{\link{MaxEntSPFBinBin}}.}
\item{SPF.Fit}{A fitted object of class \code{SPF.BinBin}. See \code{\link{SPF.BinBin}}.}
\item{Type}{The type of plot that is requested. Possible choices are: \code{Type="All.Histograms"}, the histograms of all \eqn{9} \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors arranged in a \eqn{3} by \eqn{3} grid; \code{Type="All.Densities"}, plots of densities of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors. Default \code{Type="All.Densities"}.}
\item{Col}{The color of the bins or lines when histograms or density plots are requested. Default \code{"grey"}.}
\item{\dots}{Other arguments to be passed to the \code{plot()} function.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2015a). Assessing a surrogate effect predictive value in a causal inference framework.

Alonso, A., & Van der Elst, W. (2015b). A maximum-entropy approach for the evluation of surrogate endpoints based on causal inference.
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{SPF.BinBin}}
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

# Plot results
plot(x=SPFMaxEnt, SPF.Fit=SPFSens)
}
\keyword{Plot SPF}
\keyword{Causal-Inference framework}
\keyword{Sensitivity}
\keyword{Maximum Entropy}