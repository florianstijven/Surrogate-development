\name{plot Causal-Inference BinBin}
\alias{plot Causal-Inference BinBin}
\alias{plot.ICA.BinBin}
\title{
Plots the (Meta-Analytic) Individual Causal Association and related metrics when S and T are binary outcomes
}
\description{
This function provides a plot that displays the frequencies, percentages, cumulative percentages or densities of the individual causal association (ICA; \eqn{R^2_{H}} or \eqn{R_{H}}), and/or the odds ratios for \eqn{S} and \eqn{T} (\eqn{\theta_{S}} and \eqn{\theta_{T}}). }
\usage{
\method{plot}{ICA.BinBin}(x, R2_H=TRUE, R_H=FALSE, Theta_T=FALSE, 
Theta_S=FALSE, Type="Percent", Labels=FALSE, Xlab.R2_H, 
Main.R2_H, Xlab.R_H, Main.R_H, Xlab.Theta_S, Main.Theta_S, Xlab.Theta_T, 
Main.Theta_T, Cex.Legend=1, Cex.Position="topright", All.Densities=FALSE, 
col, Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), \dots)
}
\arguments{
\item{x}{An object of class \code{ICA.BinBin}. See \code{\link{ICA.BinBin}}.}
\item{R2_H}{Logical. When \code{R2_H=TRUE}, a plot of the \eqn{R^2_{H}} is provided. Default \code{TRUE}.}
\item{R_H}{Logical. When \code{R_H=TRUE}, a plot of the \eqn{R_{H}} is provided. Default \code{FALSE}.}
\item{Theta_T}{Logical. When \code{Theta_T=TRUE}, a plot of the \eqn{\theta_{T}} is provided. Default \code{FALSE}.}
\item{Theta_S}{Logical. When \code{Theta_S=TRUE}, a plot of the \eqn{\theta_{S}} is provided. Default \code{FALSE}.}
\item{Type}{The type of plot that is produced. When \code{Type="Freq"} or \code{Type="Percent"}, the Y-axis shows frequencies or percentages of \eqn{R^2_{H}}, \eqn{R_{H}}, \eqn{\theta_{T}}, or \eqn{\theta_{S}}. When \code{Type="CumPerc"}, the Y-axis shows cumulative percentages.  When \code{Type="Density"}, the density is shown. When \code{Type="All.Densities"} and the fitted object of class \code{ICA.BinBin} was obtained using a general analysis (i.e., combining the results that are obtained when no monotonicity is assumed, monotonicity is assumed for \eqn{S} alone, for \eqn{T} alone, and for both \eqn{S} and \eqn{T}, so using \code{Monotonicity=c("General")} in the function call of \code{ICA.BinBin}), the density plots are shown for the four scenarios where different assumptions regarding monotonicity are made. Default \code{"Percent"}.}  
\item{Labels}{Logical. When \code{Labels=TRUE}, the percentage of \eqn{R^2_{H}}, \eqn{R_{H}}, \eqn{\theta_{T}}, or \eqn{\theta_{S}} values that are equal to or larger than the midpoint value of each of the bins are displayed (on top of each bin). Default \code{FALSE}.}
\item{Xlab.R2_H}{The legend of the X-axis of the \eqn{R^2_{H}} plot.}
 \item{Main.R2_H}{The title of the \eqn{R^2_{H}} plot.}
\item{Xlab.R_H}{The legend of the X-axis of the \eqn{R_{H}} plot.}
 \item{Main.R_H}{The title of the \eqn{R_{H}} plot.} 
\item{Xlab.Theta_S}{The legend of the X-axis of the \eqn{\theta_{S}} plot.}
 \item{Main.Theta_S}{The title of the \eqn{\theta_{S}} plot.}
\item{Xlab.Theta_T}{The legend of the X-axis of the \eqn{\theta_{T}} plot.}
 \item{Main.Theta_T}{The title of the \eqn{\theta_{T}} plot.}
\item{Cex.Legend}{The size of the legend when \code{Type="All.Densities"} is used. Default \code{Cex.Legend=1}.}
\item{Cex.Position}{The position of the legend, \code{Cex.Position="topright"} or \code{Cex.Position="topleft"}. Default \code{Cex.Position="topright"}.}
\item{All.Densities}{
When a general analysis is requested (i.e., when \code{Monotonicity=c("General")} is used in the call of \code{\link{ICA.BinBin}}), the results of the different settings where other assumptions regarding monotonicity are made are combined (no monotonicity, monotonicity for \eqn{S} alone, for \eqn{T} alone, and for both for \eqn{S} and \eqn{T}. When \code{All.densities=TRUE} is used, a plot is provided which provides the densities of \eqn{R^2_{H}} assuming no monotonicity, assuming monotonicity for \eqn{S} alone, assuming monotonicity for \eqn{T} alone, and assuming monotonicity for both for \eqn{S} and \eqn{T}. }
  \item{col}{The color of the bins. Default \code{col <- c(8)}.}
  \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
 \item{\dots}{Extra graphical parameters to be passed to \code{hist()}.}


}
\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). A causal-inference approach for the validation of surrogate endpoints based on information theory and sensitivity analysis.}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{ICA.BinBin}}
\examples{
\dontrun{
# Generate the vector of R2_H values using the specified marginal probabilities 
# (general analysis, making no specific assumptions regarding monotonicity)

ICA <- ICA.BinBin(pi1_1_=0.2619048, pi1_0_=0.2857143, pi_1_1=0.6372549, 
pi_1_0=0.07843137, pi0_1_=0.1349206, pi_0_1=0.127451, Monotonicity=c("General"), 
Sum_Pi_f = seq(from=0.01, to=.99, by=.01), M=500, Seed=1)
           
# Plot the results (density of R2_H):
plot(ICA, Type="Density", R2_H=TRUE, R_H=FALSE, 
Theta_T=FALSE, Theta_S=FALSE)

# Plot the results as a function of assumptions regarding montonicity
plot(ICA, Type="All.Densities", R2_H=TRUE, R_H=FALSE, 
Theta_T=FALSE, Theta_S=FALSE)

}
}
\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Sensitivity}