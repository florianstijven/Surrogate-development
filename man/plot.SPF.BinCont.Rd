\name{plot.SPF.BinCont}
\alias{plot.SPF.BinCont}
\title{
Plot the surrogate predictive function (SPF) in the causal-inference single-trial setting in the binary-continuous case.
}
\description{
This function is used to create several plots related to the surrogate predictive function (SPF) in the single-trial setting within the causal-inference framework when the surrogate endpoint is continuous (normally distributed) and the true endpoint is a binary outcome. For details, see Alonso \emph{et al.} (2024).
}
\usage{
\method{plot}{SPF.BinCont}(x, Histogram.SPF=TRUE, Causal.necessity=TRUE, Best.pred=TRUE, Max.psi=TRUE, ...)
}
\arguments{
\item{x}{A fitted object of class \code{SPF.BinCont}. See \code{\link{SPF.BinCont}}.}
\item{Histogram.SPF}{Logical. Should histograms of SPF be provided? When it is requested, a matrix of histograms illustrating various combination of the SPF, i.e., the \eqn{P[\Delta T | \Delta S \in I_{ab}]}, will be produced. Default \code{Histogram.SPF=TRUE}.}
\item{Causal.necessity}{Logical. Should a histogram showing the \eqn{P[\Delta T = 0 | \Delta S = 0]} be provided? Default \code{Causal.necessity=TRUE}.}
\item{Best.pred}{Logical. Should a bar plot showing the frequency of \eqn{\tilde{\psi}_{ab}=i} for each interval \eqn{(x,y)} be provided? Default \code{Best.pred=TRUE}.}
\item{Max.psi}{Logical. Should a histogram showing the \eqn{P[\Delta T = \tilde{\psi}_{ab}(\Delta S)]} be provided? Default \code{Max.psi=TRUE}.}
\item{\dots}{Extra graphical parameters to be passed to \code{hist()} or \code{barplot()}.}
}
\references{
Alonso, A., Ong, F., Van der Elst, W., Molenberghs, G., & Callegaro, A. (2024). Assessing a continuous surrogate predictive value for a binary true endpoint based on causal inference and information theory in vaccine trial.
}
\author{
Fenny Ong, Wim Van der Elst, Ariel Alonso, and Geert Molenberghs
}
\seealso{
\code{\link{SPF.BinCont}}
}
\examples{
\dontrun{# Time consuming code part
data(Schizo)
fit.ica <- ICA.BinCont.BS(Dataset = Schizo, Surr = BPRS, True = PANSS_Bin, nb = 10,
Theta.S_0=c(-10,-5,5,10,10,10,10,10), Theta.S_1=c(-10,-5,5,10,10,10,10,10),
Treat=Treat, M=50, Seed=1)

fit.spf <- SPF.BinCont(fit.ica, a=-5, b=5)

summary(fit.spf)
plot(fit.spf)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{SPF}
\keyword{Plot SPF}
\keyword{Continuous endpoint}
\keyword{Binary endpoint}
