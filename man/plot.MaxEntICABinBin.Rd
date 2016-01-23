\name{plot MaxEntICA BinBin}
\alias{plot MaxEntICA BinBin}
\alias{plot.MaxEntICA.BinBin}
\title{
Plots the sensitivity-based and maximum entropy based Individual Causal Association when S and T are binary outcomes
}
\description{
This function provides a plot that displays the frequencies or densities of the individual causal association (ICA; \eqn{R^2_{H}}) as identified based on the sensitivity- (using the functions \code{\link{ICA.BinBin}}, \code{\link{ICA.BinBin.Grid.Sample}}, or \code{\link{ICA.BinBin.Grid.Full}}) and maximum entropy-based (using the function \code{\link{MaxEntICABinBin}}) approaches. }
\usage{
\method{plot}{MaxEntICA.BinBin}(x, ICA.Fit, 
Type="Density", Xlab, col, Main, ...)}
\arguments{
\item{x}{An object of class \code{MaxEntICABinBin}. See \code{\link{MaxEntICABinBin}}.}
\item{ICA.Fit}{An object of class \code{ICA.BinBin}. See \code{\link{ICA.BinBin}}.}
\item{Type}{The type of plot that is produced. When \code{Type="Freq"}, the Y-axis shows frequencies of \eqn{R^2_{H}}.  When \code{Type="Density"}, the density is shown.}
\item{Xlab}{The legend of the X-axis of the plot.}
  \item{col}{The color of the bins (frequeny plot) or line (density plot). Default \code{col <- c(8)}.}
 \item{Main}{The title of the plot.}
 \item{...}{Other arguments to be passed to \code{plot()}}
}
\references{
Alonso, A., & Van der Elst, W. (2015). A maximum-entropy approach for the evluation of surrogate endpoints based on causal inference.}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{ICA.BinBin}, \link{MaxEntICABinBin}}
\examples{
# Sensitivity-based ICA results using ICA.BinBin.Grid.Sample
ICA <- ICA.BinBin.Grid.Sample(pi1_1_=0.341, pi0_1_=0.119, pi1_0_=0.254,
pi_1_1=0.686, pi_1_0=0.088, pi_0_1=0.078, Seed=1, 
Monotonicity=c("No"), M=5000)

# Maximum-entropy based ICA
MaxEnt <- MaxEntICABinBin(pi1_1_=0.341, pi0_1_=0.119, pi1_0_=0.254,
pi_1_1=0.686, pi_1_0=0.088, pi_0_1=0.078)

# Plot results
plot(x=MaxEnt, ICA.Fit=ICA)
}
\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{BinBin}
\keyword{Maximum Entropy}