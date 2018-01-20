\name{SPF.BinCont}
\alias{SPF.BinCont}
\title{
Evaluate the surrogate predictive function (SPF) in the binary-continuous setting (sensitivity-analysis based approach)
}
\description{
Computes the surrogate predictive function (SPF) based on sensitivity-analyis, i.e., \eqn{P(\Delta T | \Delta S \in I[ab])}, in the setting where \eqn{S} is continuous and \eqn{T} is a binary endpoint. 
}
\usage{
SPF.BinCont(x, a, b)
}
\arguments{
  \item{x}{A fitted object of class \code{ICA.BinCont}.} 
  \item{a}{The lower interval \eqn{a} in \eqn{P(\Delta T | \Delta S \in I[ab])}.}
  \item{b}{The upper interval \eqn{b} in \eqn{P(\Delta T | \Delta S \in I[ab])}.}
}
%\details{}
\value{
\item{a}{The lower interval \eqn{a} in \eqn{P(\Delta T | \Delta S \in I[ab])}.}
\item{b}{The upper interval \eqn{b} in \eqn{P(\Delta T | \Delta S \in I[ab])}.}
\item{P_Delta_T_min1}{The vector of values for \eqn{P(\Delta T = -1| \Delta S \in I[ab])}.}
\item{P_Delta_T_0}{The vector of values for \eqn{P(\Delta T = 0| \Delta S \in I[ab])}.}
\item{P_Delta_T_1}{The vector of values for \eqn{P(\Delta T = 1| \Delta S \in I[ab])}.}
}
\references{
Alonso, A., Van der Elst, W., Molenberghs, G., & Verbeke, G. (2017). Assessing the predictive value of a continuous surogate for a binary true endpoint based on causal inference.
}
\author{
Wim Van der Elst & Ariel Alonso
}
\seealso{
\code{\link{ICA.BinBin}}, \code{\link{plot.SPF.BinCont}}
}
\examples{
\dontrun{ # time consuming code part
# Use ICA.BinCont to examine surrogacy
data(Schizo_BinCont)
Result_BinCont <- ICA.BinCont(M = 1000, Dataset = Schizo_BinCont,
Surr = PANSS, True = CGI_Bin, Treat=Treat, Diff.Sigma=TRUE)

# Obtain SPF
Fit <- SPF.BinCont(x=Result_BinCont, a = -30, b = -3)

# examine results
summary(Fit1)
plot(Fit1)
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{BinCont}
\keyword{Sensitivity}
\keyword{SPF}