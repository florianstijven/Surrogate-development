\name{MinSurrContCont}
\alias{MinSurrContCont}
\title{Examine the plausibility of finding a good surrogate endpoint in the Continuous-continuous case}

\description{
The function \code{MinSurrContCont} examines the plausibility of finding a good surrogate endpoint in the continuous-continuous setting. For details, see Alonso et al. (submitted).}

\usage{MinSurrContCont(T0T0, T1T1, Delta, T0T1=seq(from=0, to=1, by=.01))}
\arguments{
  \item{T0T0}{A scalar that specifies the variance of the true endpoint in the control treatment condition.}
  \item{T1T1}{A scalar that specifies the variance of the true endpoint in the experimental treatment condition.}
  \item{Delta}{A scalar that specifies an upper bound for the prediction mean squared error when predicting the individual causal effect of the treatment on the true endpoint based on the individual causal effect of the treatment on the surrogate.}
  \item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1 that should be considered in the computation of     \eqn{\rho_{min}^{2}}. Default \code{seq(0, 1, by=.1)}, i.e., the values \eqn{0}, \eqn{0.10}, \eqn{0.20}, \dots, \eqn{1}.}
  }


\value{An object of class \code{MinSurrContCont} with components, 
\item{T0T1}{A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1 that were considered (i.e., \eqn{\rho_{T_{0}T_{1}}}).}
\item{Sigma.Delta.T}{A scalar or vector that contains the standard deviations of the individual causal treatment effects on the true endpoint as a function of \eqn{\rho_{T_{0}T_{1}}}.}
\item{Rho2.Min}{A scalar or vector that contains the \eqn{\rho_{min}^{2}} values as a function of \eqn{\rho_{T_{0}T_{1}}}.}
}

\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). On the relationship between the causal-inference and meta-analytic paradigms for the validation of surrogate markers. 
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{
\code{\link{ICA.ContCont}}, \code{\link{plot Causal-Inference}, \link{plot MinSurrContCont}}
}

\examples{
# Assess the plausibility of finding a good surrogate when
# sigma_T0T0 = sigma_T1T1 = 8 and Delta = 1
MinSurr <- MinSurrContCont(T0T0 = 8, T1T1 = 8, Delta = 1)
summary(MinSurr)
\dontrun{plot(MinSurr)}
}
\keyword{Plausibility of a good surrogate}

