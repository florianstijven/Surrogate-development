\name{plot SPF BinCont}
\alias{plot SPF BinCont}
\alias{plot.SPF.BinCont}
\title{
Plots the surrogate predictive function (SPF) in the binary-continuous setting.
}
\description{
Plots the surrogate predictive function (SPF) based on sensitivity-analyis, i.e., \eqn{P(\Delta T | \Delta S \in I[ab])}, in the setting where \eqn{S} is continuous and \eqn{T} is a binary endpoint. }
\usage{
\method{plot}{SPF.BinCont}(x, Type="Frequency", Col="grey", Main, Xlab=TRUE, ...)
}
\arguments{
\item{x}{A fitted object of class \code{SPF.BinCont}. See \code{\link{ICA.BinCont}}.}
\item{Type}{The type of plot that is requested. The argument \code{Type="Frequency"} requests histograms for \eqn{P(\Delta T | \Delta S \in I[ab])}. The argument \code{Type="Percentage"} requests relative frequenties  for \eqn{P(\Delta T | \Delta S \in I[ab])}. The argument \code{Type="Most.Likely.DeltaT"} requests a histogram of the most likely \eqn{\Delta T} values. For example, when in one run of the sensitivity analysis, \eqn{P(\Delta T =-1| \Delta S \in I[ab])=.6}, \eqn{P(\Delta T =0| \Delta S \in I[ab])=.3}, and \eqn{P(\Delta T =-1| \Delta S \in I[ab])=.1}, the most likely outcome in this run would be \eqn{P(\Delta T =-1}. The argument \code{Type="Most.Likely.DeltaT"} generates a plot with percentages for the most likely \eqn{P(\Delta T)} value across all obtained values in the sensitivity analysis.} 
\item{Col}{The color of the bins or lines when histograms or density plots are requested. Default \code{"grey"}.}
\item{Main}{The title of the plot.}
\item{Xlab}{Logical. Should labels on the X-axis be shown? Default \code{Xlab=TRUE}.}
\item{\dots}{Arguments to be passed to the plot, histogram, ... functions.}
}
\references{
Alonso, A., Van der Elst, W., Molenberghs, G., & Verbeke, G. (2017). Assessing the predictive value of a continuous surogate for a binary true endpoint based on causal inference.
}
\author{
Wim Van der Elst & Ariel Alonso
}
\seealso{\link{SPF.BinCont}}
\examples{
\dontrun{ # time consuming code part
data(Schizo_BinCont)
# Use ICA.BinCont to examine surrogacy
Result_BinCont <- ICA.BinCont(M = 1000, Dataset = Schizo_BinCont,
Surr = PANSS, True = CGI_Bin, Treat=Treat, Diff.Sigma=TRUE)

# Obtain SPF
Fit <- SPF.BinCont(x=Result_BinCont, a = -30, b = -3)

# examine results
summary(Fit1)
plot(Fit1)

plot(Fit1, Type="Most.Likely.DeltaT")
}
}
\keyword{Plot SPF}
\keyword{Causal-Inference framework}
\keyword{Sensitivity}
\keyword{BinCont}