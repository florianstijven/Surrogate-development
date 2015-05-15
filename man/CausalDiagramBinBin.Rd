\name{CausalDiagramBinBin}
\alias{CausalDiagramBinBin}
\title{
Draws a causal diagram depicting the median informational coefficients of correlation (or odds ratios) between the counterfactuals for a specified range of values of the ICA in the binary-binary setting.
}
\description{
This function provides a diagram that depicts the medians of the informational coefficients of correlation (or odds ratios) between the counterfactuals for a specified range of values of the individual causal association in the binary-binary setting (\eqn{R_{H}^{2}}). 
}
\usage{
CausalDiagramBinBin(x, Values="Corrs", Theta_T0S0, Theta_T1S1, 
Min=0, Max=1, Cex.Letters=3, Cex.Corrs=2, Lines.Rel.Width=TRUE, 
Col.Pos.Neg=TRUE, Monotonicity)
}
\arguments{
  \item{x}{An object of class \code{ICA.BinBin}. See \code{\link{ICA.BinBin}}.}
  \item{Values}{Specifies whether the median informational coefficients of correlation or median odds ratios between the counterfactuals should be depicted, i.e., \code{Values="Corrs"} or \code{Values="ORs"}.}
\item{Theta_T0S0}{The odds ratio between \eqn{T} and \eqn{S} in the control group. This quantity is estimable based on the observed data. Only has to be provided when \code{Values="ORs"}.}
  \item{Theta_T1S1}{The odds ratio between \eqn{T} and \eqn{S} in the experimental treatment group. This quantity is estimable based on the observed data. Only has to be provided when \code{Values="ORs"}.}
  \item{Min}{The minimum value of \eqn{R_{H}^{2}} that should be considered. Default=\eqn{-1}.}
  \item{Max}{The maximum value of \eqn{R_{H}^{2}} that should be considered. Default=\eqn{1}.}
  \item{Cex.Letters}{The size of the symbols for the counterfactuals (\eqn{S_0}, \eqn{S_1}), \eqn{T_0}, \eqn{T_1}). Default=\eqn{3}.}
  \item{Cex.Corrs}{The size of the text depicting the median odds ratios between the counterfactuals. Default=\eqn{2}.}
  \item{Lines.Rel.Width}{Logical. When \code{Lines.Rel.Width=TRUE}, the widths of the lines that represent the odds ratios between the counterfactuals are relative to the size of the odds ratios (i.e., a smaller/thicker line is used for smaller/higher odds ratios. When \code{Lines.Rel.Width=FALSE}, the width of all lines representing the odds ratios between the counterfactuals is identical. Default=\code{TRUE}.  Only considered when \code{Values="ORs"}.}
  \item{Col.Pos.Neg}{Logical. When \code{Col.Pos.Neg=TRUE}, the color of the lines that represent the odds ratios between the counterfactuals is red for odds ratios below \eqn{1} and black for the ones above \eqn{1}. When \code{Col.Pos.Neg=FALSE}, all lines are in black. Default=\code{TRUE}. Only considered when \code{Values="ORs"}.}
\item{Monotonicity}{Specifies the monotonicity scenario that should be considered (i.e., \code{Monotonicity=c("No")}, \code{Monotonicity=c("True.Endp")}, \code{Monotonicity=c("Surr.Endp")}, or \code{Monotonicity=c("Surr.True.Endp")}).}
}

\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). On the relationship between the causal inference and meta-analytic paradigms for the validation of surrogate markers. 

Van der Elst, W., Alonso, A., & Molenberghs, G. (submitted). An exploration of the relationship between causal inference and meta-analytic measures of surrogacy.}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{\link{ICA.BinBin}}
\examples{
# Compute R2_H given the marginals specified as the pi's
ICA <- ICA.BinBin.Grid.Sample(pi1_1_=0.2619048, pi1_0_=0.2857143, 
pi_1_1=0.6372549, pi_1_0=0.07843137, pi0_1_=0.1349206, pi_0_1=0.127451,
Seed=1, Monotonicity=c("General"), M=1000)

# Obtain a causal diagram that provides the medians of the 
# correlations between the counterfactuals for the range
# of R2_H values between 0.1 and 1
   # Assume no monotonicty 
CausalDiagramBinBin(x=ICA, Min=0.1, Max=1, Monotonicity="No") 

   # Assume monotonicty for S 
CausalDiagramBinBin(x=ICA, Min=0.1, Max=1, Monotonicity="Surr.Endp") 

# Now only consider the results that were obtained when 
# monotonicity was assumed for the true endpoint
CausalDiagramBinBin(x=ICA, Values="ORs", Theta_T0S0=2.156, Theta_T1S1=10, 
Min=0, Max=1,  Monotonicity="True.Endp") 
}
\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}