\name{Sim.Data.Counterfactuals}
\alias{Sim.Data.Counterfactuals}
\title{Simulate a dataset that contains counterfactuals}
\description{
The function \code{Sim.Data.Counterfactuals} simulates a dataset that contains four (continuous) counterfactuals (i.e., potential outcomes) and a (binary) treatment indicator. The counterfactuals \eqn{T_0} and \eqn{T_1} denote the true endpoints of a patient under the control and the experimental treatments, respectively, and the counterfactuals \eqn{S_0} and \eqn{S_1} denote the surrogate endpoints of the patient under the control and the experimental treatments, respectively. The user can specify the number of patients, the desired mean values for the counterfactuals (i.e., \eqn{\bold{\mu}_c}), and the desired correlations between the counterfactuals (i.e., the off-diagonal values in the standardized \eqn{\bold{\Sigma}_c} matrix). For details, see the papers of Alonso et al. (submitted) and Van der Elst et al. (submitted).
}  
\usage{Sim.Data.Counterfactuals(N.Total=2000, mu_c=c(0, 0, 0, 0), T0S0=0, T1S1=0, 
T0T1=0, T0S1=0, T1S0=0, S0S1=0, Seed=sample(1:1000, size=1))
}
\arguments{
  \item{N.Total}{The total number of patients in the simulated dataset. Default \eqn{2000}.}
  \item{mu_c}{A vector that specifies the desired means for the counterfactuals \eqn{S_0}, \eqn{S_1}, \eqn{T_0}, and \eqn{T_1}, respectively. Default \code{c(0, 0, 0, 0)}.}
  \item{T0S0}{A scalar that specifies the desired correlation between the counterfactuals T0 and S0 that should be used in the generation of the data. Default \eqn{0}.}
  \item{T1S1}{A scalar that specifies the desired correlation between the counterfactuals T1 and S1 that should be used in the generation of the data. Default \eqn{0}.}
  \item{T0T1}{A scalar that specifies the desired correlation between the counterfactuals T0 and T1 that should be used in the generation of the data. Default \eqn{0}.}
  \item{T0S1}{A scalar that specifies the desired correlation between the counterfactuals T0 and S1 that should be used in the generation of the data. Default \eqn{0}.}
  \item{T1S0}{A scalar that specifies the desired correlation between the counterfactuals T1 and S0 that should be used in the generation of the data. Default \eqn{0}.}
  \item{S0S1}{A scalar that specifies the desired correlation between the counterfactuals T0 and T1 that should be used in the generation of the data. Default \eqn{0}.}
  \item{Seed}{A seed that is used to generate the dataset. Default \code{sample(x=1:1000, size=1)}, i.e., a random number between 1 and 1000.}
  }
\details{
The generated object \code{Data.Counterfactuals} (of class \code{data.frame}) is placed in the workspace.

The specified values for T0S0, T1S1, T0T1, T0S1, T1S0, and S0S1 in the function call should form a matrix that is positive definite (i.e., they should form a valid correlation matrix). When the user specifies values that form a matrix that is not positive definite, an error message is given and the object  \code{Data.Counterfactuals} is not generated. The function \code{\link{Pos.Def.Matrices}} can be used to examine beforehand whether a \eqn{4} by \eqn{4} matrix is positive definite. 
}

\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). On the relationship between the causal inference and meta-analytic paradigms for the validation of surrogate markers. 

Van der Elst, W., Alonso, A., & Molenberghs, G. (submitted). An exploration of the relationship between causal inference and meta-analytic measures of surrogacy. 
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs}
\seealso{
\code{\link{Sim.Data.MTS}}, \code{\link{Sim.Data.STS}}
}
\examples{
## Generate a dataset with 2000 patients, cor(S0,T0)=cor(S1,T1)=.5, 
## cor(T0,T1)=cor(T0,S1)=cor(T1,S0)=cor(S0,S1)=0, with means 
## 5, 9, 12, and 15 for S0, S1, T0, and T1, respectively:
Sim.Data.Counterfactuals(N=2000, T0S0=.5, T1S1=.5, T0T1=0, T0S1=0, T1S0=0, S0S1=0, 
mu_c=c(5, 9, 12, 15), Seed=1)
}
\keyword{Causal-Inference framework}
\keyword{Simulate data}
\keyword{Counterfactuals}