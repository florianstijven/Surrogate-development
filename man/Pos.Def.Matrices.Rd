\name{Pos.Def.Matrices}
\alias{Pos.Def.Matrices}
\title{Generate 4 by 4 correlation matrices and flag the positive definite ones}
\description{
Based on vectors (or scalars) for the six off-diagonal correlations of a \eqn{4} by \eqn{4} matrix, the function \code{Pos.Def.Matrices} constructs all possible matrices that can be formed by combining the specified values, computes the minimum eigenvalues for each of these matrices, and flags the positive definite ones (i.e., valid correlation matrices). 
}
\usage{
Pos.Def.Matrices(T0T1=seq(0, 1, by=.2), T0S0=seq(0, 1, by=.2), T0S1=seq(0, 1, 
by=.2), T1S0=seq(0, 1, by=.2), T1S1=seq(0, 1, by=.2), S0S1=seq(0, 1, by=.2))
}
\arguments{
  \item{T0T1}{A vector or scalar that specifies the correlation(s) between T0 and T1 that should be considered to construct all possible \eqn{4} by \eqn{4} matrices. Default \code{seq(0, 1, by=.2)}, i.e., the values \eqn{0}, \eqn{0.20}, \dots, \eqn{1}.}
  \item{T0S0}{A vector or scalar that specifies the correlation(s) between T0 and S0 that should be considered to construct all possible \eqn{4} by \eqn{4} matrices. Default \code{seq(0, 1, by=.2)}.}
  \item{T0S1}{A vector or scalar that specifies the correlation(s) between T0 and S1 that should be considered to construct all possible \eqn{4} by \eqn{4} matrices. Default \code{seq(0, 1, by=.2)}.}
  \item{T1S0}{A vector or scalar that specifies the correlation(s) between T1 and S0 that should be considered to construct all possible \eqn{4} by \eqn{4} matrices. Default \code{seq(0, 1, by=.2)}.}
  \item{T1S1}{A vector or scalar that specifies the correlation(s) between T1 and S1 that should be considered to construct all possible \eqn{4} by \eqn{4} matrices. Default \code{seq(0, 1, by=.2)}.}
  \item{S0S1}{A vector or scalar that specifies the correlation(s) between S0 and S1 that should be considered to construct all possible \eqn{4} by \eqn{4} matrices. Default \code{seq(0, 1, by=.2)}.}
}
\details{
The generated object \code{Generated.Matrices} (of class \code{data.frame}) is placed in the workspace (for easy access).}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs}
\seealso{
\code{\link{Sim.Data.Counterfactuals}}
}
\examples{
## Generate all 4x4 matrices that can be formed using rho(T0,S0)=rho(T1,S1)=.5
## and the grid of values 0, .2, ..., 1 for the other off-diagonal correlations: 
Pos.Def.Matrices(T0T1=seq(0, 1, by=.2), T0S0=.5, T0S1=seq(0, 1, by=.2), 
T1S0=seq(0, 1, by=.2), T1S1=.5, S0S1=seq(0, 1, by=.2))

## Examine the first 10 rows of the the object Generated.Matrices:
Generated.Matrices[1:10,]

## Check how many of the generated matrices are positive definite
## (counts and percentages):
table(Generated.Matrices$Pos.Def.Status)
table(Generated.Matrices$Pos.Def.Status)/nrow(Generated.Matrices)

## Make an object PosDef which contains the positive definite matrices:
PosDef <- Generated.Matrices[Generated.Matrices$Pos.Def.Status==1,]

## Shows the 10 first matrices that are positive definite:
PosDef[1:10,]
}