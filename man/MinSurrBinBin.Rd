\name{MinSurrBinBin}
\alias{MinSurrBinBin}
\title{Examine the plausibility of finding a good surrogate endpoint in the binary-binary case}

\description{
The function \code{MinSurrBinBin} examines the plausibility of finding a good surrogate endpoint in the binary-binary setting. For details, see Alonso et al. (submitted).}

\usage{MinSurrBinBin(Object, Monotonicity=c("No"), delta=0.5)}
\arguments{
  \item{Object}{A fitted object of class \code{\link{ICA.BinBin}}.}
  \item{Monotonicity}{The assumption that is made regarding monotonicity (which determines the cardinality of the support for |\eqn{\Delta T}|). In the most general scenario (no monotonicity assumed), |\eqn{\Delta T}|=3. Default \code{Monotonicity=c("No")}.}  
  \item{delta}{A scalar that specifies an upper bound for the prediction error
when predicting the individual causal effect of the treatment on the true endpoint
based on the individual causal effect of the treatment on the surrogate. Default \code{delta=0.5}.}
  }
\value{An object of class \code{MinSurrBinBin} with components, 
\item{R2_HL_all}{A scalar or vector that contains the \eqn{R^2_{HL}} estimates.}
}

\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). A causal-inference approach for the validation ofsurrogate endpoints based on information theory and sensitivity analysis.  
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{
\code{\link{ICA.BinBin}}, \code{\link{plot.ICA.BinBin}}}

\examples{
\dontrun{ICA <- ICA.BinBin(pi1_1_=0.2619048, pi1_0_=0.2857143, pi_1_1=0.6372549, 
pi_1_0=0.07843137, pi0_1_=0.1349206, pi_0_1=0.127451, Seed=1, 
Monotonicity=c("No"), Sum_Pi_f = seq(from=0.01, to=.99, by=.01), M=1000)

# examine the plausibility of finding a good surrogate, using
# delta=.5
MinSurrBinBin(ICA, Monotonicity=c("No"), delta=0.1)
}
}
\keyword{Plausibility of a good surrogate}

