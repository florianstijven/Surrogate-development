\name{MaxEntContCont}
\alias{MaxEntContCont}
\title{
Use the maximum-entropy approach to compute ICA in the continuous-continuous sinlge-trial setting
}
\description{
In a surrogate evaluation setting where both \eqn{S} and \eqn{T} are continuous
endpoints, a sensitivity-based approach where multiple 'plausible values' for ICA are retained can be used (see functions \code{ICA.ContCont}). The function \code{MaxEntContCont} identifies the estimate which has the maximuum entropy.}
\usage{
MaxEntContCont(x, T0T0, T1T1, S0S0, S1S1)
}
\arguments{
  \item{x}{A fitted object of class \code{ICA.ContCont}.}
  \item{T0T0}{A scalar that specifies the variance of the true endpoint in the control treatment condition.}
  \item{T1T1}{A scalar that specifies the variance of the true endpoint in the experimental treatment condition.}
  \item{S0S0}{A scalar that specifies the variance of the surrogate endpoint in the control treatment condition.}
  \item{S1S1}{A scalar that specifies the variance of the surrogate endpoint in the experimental treatment condition.}
}

\value{
\item{ICA.Max.Ent}{The ICA value with maximum entropy.}
\item{Max.Ent}{The maximum entropy.}
\item{Entropy}{The vector of entropies corresponding to the vector of 'plausible values' for ICA.}
\item{Table.ICA.Entropy}{A \code{data.frame} that contains the vector of ICA, their entropies, and the correlations between the counterfactuals.}
\item{ICA.Fit}{The fitted \code{ICA.ContCont} object.}
}
\references{
Add
}
\author{
Wim Van der Elst, Ariel Alonso, Paul Meyvisch, & Geert Molenberghs
}
\seealso{
\code{\link{ICA.ContCont}}, \code{\link{MaxEntICABinBin}}
}
\examples{
\dontrun{ #time-consuming code parts
# Compute ICA for ARMD dataset, using the grid  
# G={-1, -.80, ..., 1} for the undidentifiable correlations

ICA <- ICA.ContCont(T0S0 = 0.769, T1S1 = 0.712, S0S0 = 188.926, 
S1S1 = 132.638, T0T0 = 264.797, T1T1 = 231.771, 
T0T1 = seq(-1, 1, by = 0.2), T0S1 = seq(-1, 1, by = 0.2), 
T1S0 = seq(-1, 1, by = 0.2), S0S1 = seq(-1, 1, by = 0.2))

# Identify the maximum entropy ICA
MaxEnt_ARMD <- MaxEntContCont(x = ICA, S0S0 = 188.926, 
S1S1 = 132.638, T0T0 = 264.797, T1T1 = 231.771)

  # Explore results using summary() and plot() functions
summary(MaxEnt_ARMD)
plot(MaxEnt_ARMD)
plot(MaxEnt_ARMD, Entropy.By.ICA = TRUE)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{ContCont}
\keyword{Maximum Entropy}