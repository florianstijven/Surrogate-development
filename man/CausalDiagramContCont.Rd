\name{CausalDiagramContCont}
\alias{CausalDiagramContCont}
\title{
Draws a causal diagram depicting the median correlations between the counterfactuals for a specified range of values of ICA or MICA in the continuous-continuous setting
}
\description{
This function provides a diagram that depicts the medians of the correlations between the counterfactuals for a specified range of values of the individual causal association (ICA; \eqn{\rho_{\Delta}}) or the meta-analytic individual causal association (MICA; \eqn{\rho_{M}}). 
}
\usage{
CausalDiagramContCont(x, Min=-1, Max=1, Cex.Letters=3, Cex.Corrs=2, 
Lines.Rel.Width=TRUE, Col.Pos.Neg=TRUE)
}
\arguments{
  \item{x}{An object of class \code{ICA.ContCont} or \code{MICA.ContCont}. See \code{\link{ICA.ContCont}} or \code{\link{MICA.ContCont}}.}
  \item{Min}{The minimum values of (M)ICA that should be considered. Default=\eqn{-1}.}
  \item{Max}{The maximum values of (M)ICA that should be considered. Default=\eqn{1}.}
  \item{Cex.Letters}{The size of the symbols for the counterfactuals (\eqn{S_0}, \eqn{S_1}), \eqn{T_0}, \eqn{T_1}). Default=\eqn{3}.}
  \item{Cex.Corrs}{The size of the text depicting the median correlations between the counterfactuals. Default=\eqn{2}.}
  \item{Lines.Rel.Width}{Logical. When \code{Lines.Rel.Width=TRUE}, the widths of the lines that represent the correlations between the counterfactuals are relative to the size of the correlations (i.e., a smaller line is used for correlations closer to zero whereas a thicker line is used for (absolute) correlations closer to \eqn{1}). When \code{Lines.Rel.Width=FALSE}, the width of all lines representing the correlations between the counterfactuals is identical. Default=\code{TRUE}.}
  \item{Col.Pos.Neg}{Logical. When \code{Col.Pos.Neg=TRUE}, the color of the lines that represent the correlations between the counterfactuals is red for negative correlations and black for positive ones. When \code{Col.Pos.Neg=FALSE}, all lines are in black. Default=\code{TRUE}.}
}
\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). On the relationship between the causal inference and meta-analytic paradigms for the validation of surrogate markers. 

Van der Elst, W., Alonso, A., & Molenberghs, G. (submitted). An exploration of the relationship between causal inference and meta-analytic measures of surrogacy.}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{\link{ICA.ContCont}, \link{MICA.ContCont}}
\examples{
\dontrun{ #Time consuming (>5 sec) code parts
# Generate the vector of ICA values when rho_T0S0=.91, rho_T1S1=.91, and when the
# grid of values {0, .1, ..., 1} is considered for the correlations
# between the counterfactuals:
SurICA <- ICA.ContCont(T0S0=.95, T1S1=.91, T0T1=seq(0, 1, by=.1), T0S1=seq(0, 1, by=.1), 
T1S0=seq(0, 1, by=.1), S0S1=seq(0, 1, by=.1))

#obtain a plot of ICA

# Obtain a causal diagram that provides the medians of the 
# correlations between the counterfactuals for the range
# of ICA values between .9 and 1 (i.e., which assumed 
# correlations between the counterfactuals lead to a 
# high ICA?)
CausalDiagramContCont(SurICA, Min=.9, Max=1)

# Same, for low values of ICA
CausalDiagramContCont(SurICA, Min=0, Max=.5)}
}
\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Multiple-trial setting}
