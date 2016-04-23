\name{plot Causal-Inference BinCont}
\alias{plot Causal-Inference BinCont}
\alias{plot.ICA.BinCont}
\title{Plots the (Meta-Analytic) Individual Causal Association and related metrics when S is continuous and T is binary}
\description{
This function provides a plot that displays the frequencies, percentages, cumulative percentages or densities of the individual causal association (ICA; \eqn{R^2_{H}}) in the setting where S is continuous and T is binary. }
\usage{
\method{plot}{ICA.BinCont}(x, Xlab, Main=" ", Type="Density", Labels=TRUE, \dots)
}
\arguments{
\item{x}{An object of class \code{ICA.BinCont}. See \code{\link{ICA.BinCont}}.}
\item{Xlab}{The legend of the X-axis of the plot.}
\item{Main}{The title of the plot.}
\item{Type}{The type of plot that is produced. When \code{Type="Freq"} or \code{Type="Percent"}, the Y-axis shows frequencies or percentages of \eqn{R^2_{H}}. When \code{Type="CumPerc"}, the Y-axis shows cumulative percentages.  When \code{Type="Density"}, the density is shown}. 
\item{Labels}{Logical. When \code{Labels=TRUE}, the percentage of \eqn{R^2_{H}} values that are equal to or larger than the midpoint value of each of the bins are displayed (on top of each bin). Default \code{FALSE}.}
 \item{\dots}{Extra graphical parameters to be passed to \code{hist()}.}
}
\references{
Alonso, A., Van der Elst, W., & Meyvisch, P. (2016). Surrogate markers validation: the continuous-binary setting from a causal inference perspective.}
\author{
Wim Van der Elst, Paul Meyvisch, & Ariel Alonso
}
\seealso{\link{ICA.BinCont}}
\examples{
\dontrun{# Time consuming code part
Fit <- ICA.BinCont(Dataset = Schizo, Surr = BPRS, True = PANSS_Bin, 
Treat=Treat, M=50, Seed=1)

summary(Fit)
plot(Fit)
}
}
\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Sensitivity}