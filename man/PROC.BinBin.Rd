\name{PROC.BinBin}
\alias{PROC.BinBin}
\title{Evaluate the individual causal association (ICA) and reduction in probability of a prediction error (RPE) in the setting where both \eqn{S} and \eqn{T} are binary endpoints}
\description{
The function \code{PROC.BinBin} assesses the ICA and RPE in the single-trial causal-inference framework when both the surrogate and the true endpoints are binary outcomes. It additionally allows to account for sampling variability by means of bootstrap. See \bold{Details} below.}
\usage{
PROC.BinBin(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat, 
BS=FALSE, seqs=250, MC_samples=1000, Seed=1)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a binary surrogate value, a binary true endpoint value, and a treatment indicator. }
  \item{Surr}{The name of the variable in \code{Dataset} that contains the binary surrogate endpoint values. Should be coded as \eqn{0} and \eqn{1}.}
  \item{True}{The name of the variable in \code{Dataset} that contains the binary true endpoint values. Should be coded as \eqn{0} and \eqn{1}.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group. }
   \item{BS}{Logical. If \code{TRUE}, then Dataset will be bootstrapped to account for sampling variability. If \code{FALSE}, then no bootstrap is performed. See the \bold{Details} section below. Default \code{FALSE}. }
  \item{seqs}{The number of copies of the dataset that are produced or alternatively the number of bootstrap datasets that are produced. Default \code{seqs=250}. }
   \item{MC_samples}{The number of Monte Carlo samples that need to be obtained per copy of the data set. Default \code{MC_samples=1000}. }
  \item{Seed}{The seed to be used. Default \code{Seed=1}.}

}

\details{

In the continuous normal setting, surroagacy can be assessed by studying the association between the individual causal effects on \eqn{S} and \eqn{T} (see \code{\link{ICA.ContCont}}). In that setting, the Pearson correlation is the obvious measure of association.

When \eqn{S} and \eqn{T} are binary endpoints, multiple alternatives exist. Alonso et al. (2016) proposed the individual causal association (ICA; \eqn{R_{H}^{2}}), which captures the association between the individual causal effects of the treatment on \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T}) using information-theoretic principles. 

The function \code{PPE.BinBin} computes \eqn{R_{H}^{2}} using a grid-based approach where all possible combinations of the specified grids for the parameters that are allowed to vary freely are considered. It additionally computes the minimal probability of a prediction error (PPE) and the reduction on the PPE using information that \eqn{S} conveys on \eqn{T} (RPE). Both measures provide complementary information over the \eqn{R_{H}^{2}} and facilitate more straightforward clinical interpretation. No assumption about monotonicity can be made. The function \code{PROC.BinBin} makes direct use of the function \code{PPE.BinBin}. However, it is computationally much faster thanks to equally dividing the number of Monte Carlo samples over copies of the input data. In addition, it allows to account for sampling variability using a bootstrap procedure. Finally, the function \code{PROC.BinBin} computes the marginal probabilities directly from the input data set. 
}


\value{An object of class \code{PPE.BinBin} with components,
\item{PPE}{The vector of the PPE values.}
\item{RPE}{The vector of the RPE values.}
\item{PPE_T}{The vector of the \eqn{PPE_T} values indicating the probability on a prediction error without using information on \eqn{S}.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}


}
\references{
Alonso A,  Van der Elst W, Molenberghs G, Buyse M and Burzykowski T. (2016). An information-theoretic approach for the evaluation of surrogate endpoints based on causal inference. 

Meyvisch P., Alonso A.,Van der Elst W, Molenberghs G.. Assessing the predictive value of a binary surrogate for a binary true endpoint, based on the minimum probability of a prediction error.
}
\author{
Paul Meyvisch, Wim Van der Elst, Ariel Alonso, Geert Molenberghs
}
\seealso{
\code{\link{PPE.BinBin}}
}
\examples{
# Conduct the analysis 
 
\dontrun{ # time consuming code part
library(Surrogate)
# load the CIGTS data 
data(CIGTS)
CIGTS_25000<-PROC.BinBin(Dataset=CIGTS, Surr=IOP_12, True=IOP_96, 
Treat=Treat, BS=FALSE,seqs=250, MC_samples=100, Seed=1)
}

}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{ICA}
\keyword{RPE}
\keyword{PPE}