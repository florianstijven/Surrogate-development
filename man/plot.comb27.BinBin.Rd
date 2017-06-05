\name{plot.comb27.BinBin}
\alias{plot.comb27.BinBin}
\title{Plots the distribution of prediction error functions in decreasing order of appearance.}
\description{
The function \code{plot.comb27.BinBin} plots each of the selected prediction functions in decreasing order in the single-trial causal-inference framework when both the surrogate and the true endpoints are binary outcomes. The distribution of frequencies at which each of the 27 possible predicton functions are selected provides additional insights regarding the association between \eqn{S} (\eqn{\Delta_S}) and \eqn{T} (\eqn{\Delta_T})..  See \bold{Details} below.}
\usage{
\method{plot}{comb27.BinBin}(x,lab,...)
}
\arguments{
  \item{x}{An object of class \code{comb27.BinBin}. See \code{\link{comb27.BinBin}}.}
  \item{lab}{a supplementary label to the graph.}
  \item{...}{Other arguments to be passed}
}

\details{
Each of the 27 prediction functions is coded as x/y/z with x, y and z taking values in \eqn{{-1,0,1}}. As an example, the combination 0/0/0 represents the prediction function that projects every value of \eqn{\Delta_S} to 0. Similarly, the combination -1/0/1 is the identity function projecting every value of \eqn{\Delta_S} to the same value for \eqn{\Delta_T}.  

}

\value{An object of class \code{comb27.BinBin} with components,
\item{index}{count variable}
\item{Monotonicity}{The vector of Monotonicity assumptions}
\item{Pe}{The vector of the prediction error values.}
\item{combo}{The vector containing the codes for the each of the 27 prediction functions.}
\item{R2_H}{The vector of the \eqn{R_H^2} values.}
\item{H_Delta_T}{The vector of the entropies of \eqn{\Delta_T}.}
\item{H_Delta_S}{The vector of the entropies of \eqn{\Delta_S}.}
\item{I_Delta_T_Delta_S}{The vector of the mutual information of \eqn{\Delta_S} and \eqn{\Delta_T}.}


}
\references{
Alonso A,  Van der Elst W, Molenberghs G, Buyse M and Burzykowski T. (2016). An information-theoretic approach for the evaluation of surrogate endpoints based on causal inference. 

Alonso A, Van der Elst W and Meyvisch P (2016). Assessing a surrogate predictive value: A causal inference approach. 
}
\author{
Paul Meyvisch, Wim Van der Elst, Ariel Alonso
}
\seealso{
\code{\link{comb27.BinBin}}
}
\examples{
\dontrun{ # time consuming code part
CIGTS_27 <- comb27.BinBin(pi1_1_ = 0.3412, pi1_0_ = 0.2539, pi0_1_ = 0.119, 
                       pi_1_1 = 0.6863, pi_1_0 = 0.0882, pi_0_1 = 0.0784,  
                       Seed=1,Monotonicity=c("No"), M=500000) 
plot.comb27.BinBin(CIGTS_27,lab="CIGTS")
}
}
\keyword{Causal-Inference framework}
\keyword{Counterfactuals}
\keyword{Single-trial setting}
\keyword{ICA}
\keyword{PPE}