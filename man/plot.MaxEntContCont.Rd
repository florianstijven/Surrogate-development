\name{plot MaxEnt ContCont}
\alias{plot MaxEnt ContCont}
\alias{plot.MaxEntContCont}
\title{
Plots the sensitivity-based and maximum entropy based Individual Causal Association when S and T are continuous outcomes in the single-trial setting
}
\description{
This function provides a plot that displays the frequencies or densities of the individual causal association (ICA; \eqn{rho[Delta]}) as identified based on the sensitivity- (using the functions \code{\link{ICA.ContCont}}) and maximum entropy-based (using the function \code{\link{MaxEntContCont}}) approaches. }
\usage{
\method{plot}{MaxEntContCont}(x, Type="Freq", Xlab, col, 
Main, Entropy.By.ICA=FALSE, ...)}
\arguments{
\item{x}{An object of class \code{MaxEntContCont}. See \code{\link{MaxEntContCont}}.}
\item{Type}{The type of plot that is produced. When \code{Type="Freq"}, the Y-axis shows frequencies of ICA.  When \code{Type="Density"}, the density is shown. Default \code{Type="Freq"}.}
\item{Xlab}{The legend of the X-axis of the plot.}
\item{col}{The color of the bins (frequeny plot) or line (density plot). Default \code{col <- c(8)}.}
\item{Main}{The title of the plot.}
\item{Entropy.By.ICA}{Plot with ICA on Y-axis and entropy on X-axis.}
\item{...}{Other arguments to be passed to \code{plot()}}
}
\references{
Add
}
\author{
Wim Van der Elst, Ariel Alonso, Paul Meyvisch, & Geert Molenberghs
}
\seealso{\link{ICA.ContCont}, \link{MaxEntContCont}}
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
\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}
\keyword{Sensitivity}
\keyword{ContCont}
\keyword{Maximum Entropy}