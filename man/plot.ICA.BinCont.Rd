\name{plot.ICA.BinCont}
\alias{plot.ICA.BinCont}
\title{Plot the individual causal association (ICA) in the causal-inference single-trial setting in the binary-continuous case.}
\description{
This function is used to a plot that displays the frequencies, percentages, cumulative percentages or densities of the individual causal association (ICA) in the single-trial setting within the causal-inference framework when the surrogate endpoint is continuous (normally distributed) and the true endpoint is a binary outcome. In addition, several plots to evaluate the goodness-of-fit of the mixture model used to fit the conditional distribution of potential outcomes on the surrogate endpoint can also be provided. For details, see Alonso Abad \emph{et al.} (2023).
}
\usage{
\method{plot}{ICA.BinCont}(x, Histogram.ICA=TRUE, Mixmean=TRUE, Mixvar=TRUE, Deviance=TRUE,
                             Type="Percent", Labels=FALSE, ...)
}
\arguments{
\item{x}{A fitted object of class \code{ICA.BinCont}. See \code{\link{ICA.BinCont}} or \code{\link{ICA.BinCont.BS}}.}
\item{Histogram.ICA}{Logical. Should a histogram of ICA be provided? Default \code{Histogram.ICA=TRUE}.}
\item{Mixmean}{Logical. Should a plot of the calculated means of the fitted mixtures for \eqn{S_{0}} and \eqn{S_{1}} across different iterations be provided? Default \code{Mixmean=TRUE}.}
\item{Mixvar}{Logical. Should a plot of the calculated variances of the fitted mixtures for \eqn{S_{0}} and \eqn{S_{1}} across different iterations be provided? Default \code{Mixvar=TRUE}.}
\item{Deviance}{Logical. Should a boxplot of the deviances for the fitted mixtures of \eqn{S_{0}} and \eqn{S_{1}} be provided? Default \code{Deviance=TRUE}.}
\item{Type}{The type of plot that is produced for the histogram of ICA. When \code{Type="Freq"} or \code{Type="Percent"}, the Y-axis shows frequencies or percentages of \eqn{R^2_{H}}. When \code{Type="CumPerc"}, the Y-axis shows cumulative percentages. When \code{Type="Density"}, the density is shown.}
\item{Labels}{Logical. When \code{Labels=TRUE}, the percentage of \eqn{R^2_{H}} values that are equal to or larger than the midpoint value of each of the bins are added in the histogram of ICA (on top of each bin). Default \code{Labels=FALSE}.}
 \item{\dots}{Extra graphical parameters to be passed to \code{plot()} or \code{hist()}.}
}
\references{
Alonso Abad, A., Ong, F., Stijven, F., Van der Elst, W., Molenberghs, G., Van Keilegom, I., Verbeke, G., & Callegaro, A. (2023). An information-theoretic approach for the assessment of a continuous outcome as a surrogate for a binary true endpoint based on causal inference: Application to vaccine evaluation.
}
\author{
Wim Van der Elst, Fenny Ong, Ariel Alonso, and Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinCont}}, \code{\link{ICA.BinCont.BS}}
}
\examples{
\dontrun{# Time consuming code part
data(Schizo)
Fit <- ICA.BinCont.BS(Dataset = Schizo, Surr = BPRS, True = PANSS_Bin, nb = 10,
Theta.S_0=c(-10,-5,5,10,10,10,10,10), Theta.S_1=c(-10,-5,5,10,10,10,10,10),
Treat=Treat, M=50, Seed=1)

summary(Fit)
plot(Fit)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ICA}
\keyword{Plot ICA}
\keyword{Continuous endpoint}
\keyword{Binary endpoint}
