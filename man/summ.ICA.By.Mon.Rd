\name{summ.ICA.By.Mon}
\alias{summ.ICA.By.Mon}

\title{
Provide summary measures for ICA in the bin-bin setting for each of the monotonicity scenarios 
}
\description{
This functions provide summary measures (e.g., mean, median, SD, ...) for ICA in the bin-bin setting for each of the monotonicity scenarios.
}
\usage{
summ.ICA.By.Mon(x)
}
\arguments{
  \item{x}{An object of class \code{ICA.BinBin}}
}
\author{
Wim Van der Elst, Ariel Alonso, Marc Buyse, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinBin}}, \code{\link{ICA.BinBin.Grid.Full}}, \code{\link{ICA.BinBin.Grid.Sample}}
}
\examples{
\dontrun{
# Generate the vector of R2_H values using the specified marginal probabilities 
# (general analysis, making no specific assumptions regarding monotonicity)

ICA <- ICA.BinBin.Grid.Sample(pi1_1_=0.2619048, pi1_0_=0.2857143, pi_1_1=0.6372549, 
pi_1_0=0.07843137, pi0_1_=0.1349206, pi_0_1=0.127451, Monotonicity=c("General"), 
M=5000, Seed=1)

#obtain summary of results, per monotonicity scenario seperatly
summ.ICA.By.Mon(ICA)
}}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Sensitivity}