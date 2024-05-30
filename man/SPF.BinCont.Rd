\name{SPF.BinCont}
\alias{SPF.BinCont}
\title{Evaluate the surrogate predictive function (SPF) in the causal-inference single-trial setting in the binary-continuous case}
\description{
The function \code{SPF.BinCont} computes the surrogate predictive function (SPF), i.e., the \eqn{P[\Delta T | \Delta S \in I_{ab}]} in the single-trial setting within the causal-inference framework when the surrogate endpoint is continuous (normally distributed) and the true endpoint is a binary outcome. For details, see Alonso \emph{et al.} (2024).
}
\usage{
SPF.BinCont(x, a, b)
}
\arguments{
  \item{x}{A fitted object of class \code{ICA.BinCont}.}
  \item{a}{The lower interval \eqn{a} in \eqn{P[\Delta T | \Delta S \in I_{ab}]}.}
  \item{b}{The upper interval \eqn{b} in \eqn{P[\Delta T | \Delta S \in I_{ab}]}.}
}
\value{An object of class \code{SPF.BinCont} with important or relevant components:
\item{a}{The lower interval \eqn{a} in \eqn{P[\Delta T | \Delta S \in I_{ab}]}.}
\item{b}{The upper interval \eqn{b} in \eqn{P[\Delta T | \Delta S \in I_{ab}]}.}
\item{r_min1_min1}{The vector of \eqn{P[\Delta T = -1 | \Delta S \in I_{(-\infty,a)}]}.}
\item{r_0_min1}{The vector of \eqn{P[\Delta T = 0 | \Delta S \in I_{(-\infty,a)}]}.}
\item{r_1_min1}{The vector of \eqn{P[\Delta T = 1 | \Delta S \in I_{(-\infty,a)}]}.}
\item{r_min1_0}{The vector of \eqn{P[\Delta T = -1 | \Delta S \in I_{(a,b)}]}.}
\item{r_0_0}{The vector of \eqn{P[\Delta T = 0 | \Delta S \in I_{(a,b)}]}.}
\item{r_1_0}{The vector of \eqn{P[\Delta T = 1 | \Delta S \in I_{(a,b)}]}.}
\item{r_min1_1}{The vector of \eqn{P[\Delta T = -1 | \Delta S \in I_{(b,\infty)}]}.}
\item{r_0_1}{The vector of \eqn{P[\Delta T = 0 | \Delta S \in I_{(b,\infty)}]}.}
\item{r_1_1}{The vector of \eqn{P[\Delta T = 1 | \Delta S \in I_{(b,\infty)}]}.}
\item{P_DT_0_DS_0}{The vector of \eqn{P[\Delta T = 0 | \Delta S = 0]}.}
\item{P_DT_psi_DS_max}{The vector of \eqn{P[\Delta T = \tilde{\psi}_{ab}(\Delta S)]}, where \eqn{\tilde{\psi}_{ab}(\Delta S)=arg max_{i}P[\Delta T=i|\Delta S \in (x,y)]}.}
\item{best.pred.min1}{The vector of \eqn{\tilde{\psi}_{ab}(\Delta S)=arg max_{i}P[\Delta T=i|\Delta S \in (x,y)]}, where \eqn{(x,y)=(-\infty,a)}.}
\item{best.pred.0}{The vector of \eqn{\tilde{\psi}_{ab}(\Delta S)=arg max_{i}P[\Delta T=i|\Delta S \in (x,y)]}, where \eqn{(x,y)=(a,b)}.}
\item{best.pred.1}{The vector of \eqn{\tilde{\psi}_{ab}(\Delta S)=arg max_{i}P[\Delta T=i|\Delta S \in (x,y)]}, where \eqn{(x,y)=(b,\infty)}.}
}
\references{
Alonso, A., Ong, F., Van der Elst, W., Molenberghs, G., & Callegaro, A. (2024). Assessing a continuous surrogate predictive value for a binary true endpoint based on causal inference and information theory in vaccine trial.
}
\author{
Fenny Ong, Wim Van der Elst, Ariel Alonso, and Geert Molenberghs
}
\seealso{
\code{\link{ICA.BinCont}}, \code{\link{ICA.BinCont.BS}}, \code{\link{plot.SPF.BinCont}}
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
\keyword{ICA}
\keyword{SPF}
\keyword{Continuous endpoint}
\keyword{Binary endpoint}
