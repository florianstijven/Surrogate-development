\name{plot.Fano.BinBin}
\alias{plot.Fano.BinBin}
\title{Plots the distribution of \eqn{R^2_{HL}} either as a density or as function of \eqn{\pi_{10}} in the setting where both \eqn{S} and \eqn{T} are binary endpoints}
\description{
The function \code{plot.Fano.BinBin} plots the distribution of \eqn{R^2_{HL}} which is fully identifiable for given values of \eqn{\pi_{10}}.  See \bold{Details} below.}
\usage{
\method{plot}{Fano.BinBin}(x,Type="Density",Xlab.R2_HL,main.R2_HL,
ylab="density",Par=par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5.1,4.1,4.1,2.1)),
Cex.Legend=1,Cex.Position="top", lwd=3,linety=c(5,6,7),color=c(8,9,3),...)
}
\arguments{
  \item{x}{An object of class \code{Fano.BinBin}. See \code{\link{Fano.BinBin}}.}
  \item{Type}{The type of plot that is produced. When \code{Type="Freq"}, a histogram of \eqn{R^2_{HL}} is produced. When \code{Type="Density"}, the density of \eqn{R^2_{HL}} is produced. When \code{Type="Scatter"}, a scatter plot of \eqn{R^2_{HL}} is produced as a function of \eqn{\pi_{10}}. Default \code{Type="Scatter"}.} 
  \item{Xlab.R2_HL}{The label of the X-axis when density plots or histograms are produced.}
  \item{main.R2_HL}{Title of the density plot or histogram.}
  \item{ylab}{The label of the Y-axis when density plots or histograms are produced. Default \code{ylab="density"}. }
  \item{Par}{Graphical parameters for the plot. Default \code{par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5.1,4.1,4.1,2.1))}.}
  \item{Cex.Legend}{The size of the legend. Default \code{Cex.Legend=1}.}
  \item{Cex.Position}{The position of the legend.  Default \code{Cex.Position="top"}.}
  \item{lwd}{The line width for the density plot . Default \code{lwd=3}.}
  \item{linety}{The line types corresponding to each level of \code{fano_delta} . Default \code{linety=c(5,6,7)}.}
  \item{color}{The color corresponding to each level of \code{fano_delta} . Default \code{color=c(8,9,3)}.}
  \item{...}{Other arguments to be passed.}
}

\details{
Values for \eqn{\pi_{10}} have to be uniformly sampled from the interval \eqn{[0,\min(\pi_{1\cdot},\pi_{\cdot0})]}. Any sampled value for \eqn{\pi_{10}} will fully determine the bivariate distribution of potential outcomes for the true endpoint.

The vector \eqn{\bold{\pi_{km}}} fully determines \eqn{R^2_{HL}}.
}

\value{An object of class \code{Fano.BinBin} with components,
\item{R2_HL}{The sampled values for \eqn{R^2_{HL}}.}
\item{H_Delta_T}{The sampled values for \eqn{H{\Delta T}}.}
\item{minpi10}{The minimum value for \eqn{\pi_{10}}.}
\item{maxpi10}{The maximum value for \eqn{\pi_{10}}.}
\item{samplepi10}{The sampled value for \eqn{\pi_{10}}.}
\item{delta}{The specified vector of upper bounds for the prediction errors.}
\item{uncertainty}{Indexes the sampling of \eqn{pi1\_}.}
\item{pi_00}{The sampled values for \eqn{\pi_{00}}.}
\item{pi_11}{The sampled values for \eqn{\pi_{11}}.}
\item{pi_01}{The sampled values for \eqn{\pi_{01}}.}
\item{pi_10}{The sampled values for \eqn{\pi_{10}}.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2014). Validation of surrogate endpoints: the binary-binary setting from a causal inference perspective.
}
\author{
Paul Meyvisch, Wim Van der Elst, Ariel Alonso
}
\seealso{
\code{\link{Fano.BinBin}}
}
\examples{
# Conduct the analysis assuming no montonicity
# for the true endpoint, using a range of
# upper bounds for prediction errors 
FANO<-Fano.BinBin(pi1_ = 0.5951 ,  pi_1 = 0.7745, 
fano_delta=c(0.05, 0.1, 0.2), M=1000)

plot(FANO, Type="Scatter",color=c(3,4,5),Cex.Position="bottom")
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Fano}
\keyword{ICA}
\keyword{MarginalProbs}