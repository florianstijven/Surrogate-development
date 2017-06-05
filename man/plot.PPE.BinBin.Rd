\name{plot.PPE.BinBin}
\alias{plot.PPE.BinBin}
\title{Plots the distribution of either \eqn{PPE}, \eqn{RPE} or \eqn{R^2_{H}} either as a density or as a histogram in the setting where both \eqn{S} and \eqn{T} are binary endpoints}
\description{
The function \code{plot.PPE.BinBin} plots the distribution of \eqn{PPE}, \eqn{RPE} or \eqn{R^2_{H}} in the setting where both surrogate and true endpoints are binary in the single-trial causal-inference framework.  See \bold{Details} below.}
\usage{
\method{plot}{PPE.BinBin}(x,Type="Density",Param="PPE",
Xlab.PE,main.PE,ylab="density",
Par=par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5.1,4.1,4.1,2.1)),
Cex.Legend=1,Cex.Position="bottomright", lwd=3,linety=1,color=1,
Breaks=0.05,xlimits=c(0,1), ...)
}
\arguments{
  \item{x}{An object of class \code{PPE.BinBin}. See \code{\link{PPE.BinBin}}.}
  \item{Type}{The type of plot that is produced. When \code{Type="Freq"}, a histogram is produced. When \code{Type="Density"}, a density is produced.  Default \code{Type="Density"}.} 
  \item{Param}{The parameter that is to be plotted:  \code{Param="PPE"}, \code{Param="RPE"} or \code{Param="ICA"}. Default \code{Param="PPE"}.}
  \item{Xlab.PE}{The label of the X-axis when density plots or histograms are produced.}
  \item{main.PE}{Title of the density plot or histogram.}
  \item{ylab}{The label of the Y-axis for the density plots. Default \code{ylab="density"}. }
  \item{Par}{Graphical parameters for the plot. Default \code{par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5.1,4.1,4.1,2.1))}.}
  \item{Cex.Legend}{The size of the legend. Default \code{Cex.Legend=1}.}
  \item{Cex.Position}{The position of the legend.  Default \code{Cex.Position="bottomright"}.}
  \item{lwd}{The line width for the density plot. Default \code{lwd=3}.}
  \item{linety}{The line types for the density. Default \code{linety=1}.}
  \item{color}{The color of the density or histogram. Default \code{color=1}.}
  \item{Breaks}{The breaks for the histogram. Default \code{Breaks=0.05}.}
  \item{xlimits}{The limits for the X-axis. Default \code{xlimits=c(0,1)}.}
  \item{...}{Other arguments to be passed.}
}

\details{
In the continuous normal setting, surroagacy can be assessed by studying the association between the individual causal effects on \eqn{S} and \eqn{T} (see \code{\link{ICA.ContCont}}). In that setting, the Pearson correlation is the obvious measure of association.

When \eqn{S} and \eqn{T} are binary endpoints, multiple alternatives exist. Alonso et al. (2016) proposed the individual causal association (ICA; \eqn{R_{H}^{2}}), which captures the association between the individual causal effects of the treatment on \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T}) using information-theoretic principles. 

The function \code{PPE.BinBin} computes \eqn{R_{H}^{2}} using a grid-based approach where all possible combinations of the specified grids for the parameters that are allowed that are allowed to vary freely are considered. It additionally computes the minimal probability of a prediction error (PPE) and the reduction on the PPE using information that \eqn{S} conveys on \eqn{T}. Both measures provide complementary information over the \eqn{R_{H}^{2}} and facilitate more straightforward clinical interpretation. 
}

\value{An object of class \code{PPE.BinBin} with components,
\item{index}{count variable}
\item{Monotonicity}{The vector of Monotonicity assumptions}
\item{PPE}{The vector of the PPE values.}
\item{RPE}{The vector of the RPE values.}
\item{PPE_T}{The vector of the \eqn{PPE_T} values indicating the probability on a prediction error without using information on \eqn{S}.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{H_Delta_T}{The vector of the entropies of \eqn{\Delta_T}.}
\item{H_Delta_S}{The vector of the entropies of \eqn{\Delta_S}.}
\item{I_Delta_T_Delta_S}{The vector of the mutual information of \eqn{\Delta_S} and \eqn{\Delta_T}.}
\item{Pi.Vectors}{An object of class \code{data.frame} that contains the valid \eqn{\pi} vectors.}

}
\references{
Alonso A,  Van der Elst W, Molenberghs G, Buyse M and Burzykowski T. (2016). An information-theoretic approach for the evaluation of surrogate endpoints based on causal inference. 

Alonso A, Van der Elst W and Meyvisch P (2016). Assessing a surrogate predictive value: A causal inference approach. 
}
\author{
Paul Meyvisch, Wim Van der Elst, Ariel Alonso
}
\seealso{
\code{\link{PPE.BinBin}}
}
\examples{

\dontrun{ # time consuming code part
PANSS <- PPE.BinBin(pi1_1_=0.4215, pi0_1_=0.0538, pi1_0_=0.0538,
                   pi_1_1=0.5088, pi_1_0=0.0307,pi_0_1=0.0482, 
                   Seed=1,Monotonicity=c( "No"), M=1000000) 

plot(PANSS,Type="Freq",Param="ICA",color="grey",Breaks=0.02,xlimits=c(0,0.8),main="PANSS")
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{ICA}
\keyword{PPE}