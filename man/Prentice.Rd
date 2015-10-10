\name{Prentice}
\alias{Prentice}
\title{Evaluates surrogacy based on the Prentice criteria for continuous endpoints (single-trial setting)}
\description{
The function \code{Prentice} evaluates the validity of a potential surrogate based on the Prentice criteria (Prentice, 1989) in the setting where the candidate surrogate and the true endpoint are normally distributed endpoints.

\bold{Warning} The Prentice approach is included in the \emph{Surrogate} package for illustrative purposes (as it was the first formal approach to assess surrogacy), but this method has some severe problems that renders its use problematic (see \bold{Details} below). It is recommended to replace the Prentice approach by a more statistically-sound approach to evaluate a surrogate (e.g., the meta-analytic methods; see the functions \code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}}, \code{\link{BimixedContCont}}).  
}

\usage{
Prentice(Dataset, Surr, True, Treat, Pat.ID, Alpha=.05)
}
\arguments{
  \item{Dataset}{A \code{data.frame} that should consist of one line per patient. Each line contains (at least) a surrogate value, a true endpoint value, a treatment indicator, a patient ID, and a trial ID.}
  \item{Surr}{The name of the variable in \code{Dataset} that contains the surrogate values.}
  \item{True}{The name of the variable in \code{Dataset} that contains the true endpoint values.}
  \item{Treat}{The name of the variable in \code{Dataset} that contains the treatment indicators. The treatment indicator should either be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group, or as \eqn{1} for the experimental group and \eqn{0} for the control group.}
  \item{Pat.ID}{The name of the variable in \code{Dataset} that contains the patient's ID.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to examine whether the Prentice criteria are fulfilled. Default \eqn{0.05}.}  
}
\details{
The Prentice criteria are examined by fitting the following regression models (when the surrogate and true endpoints are continuous variables):
\deqn{S_{j}=\mu_{S}+\alpha Z_{j}+\varepsilon_{Sj}, (1)} 
\deqn{T_{j}=\mu_{T}+\beta Z_{j}+\varepsilon_{Tj}, (2)} 
\deqn{T_{j}=\mu+\gamma Z_{j}+\varepsilon_{j}, (3)} 
\deqn{T_{j}=\tilde{\mu}_{T}+\beta_{S} Z_{j}+\gamma_{Z} S_{j}+\tilde{\varepsilon}_{Tj}, (4)} 

where the error terms of (1) and (2) have a joint zero-mean normal distribution with variance-covariance matrix

\deqn{\boldsymbol{\Sigma}=\left(\begin{array}{cc}
\sigma_{SS}\\
\sigma_{ST} & \sigma_{TT}
\end{array}\right)},

and where \eqn{j} is the subject indicator, \eqn{S_{j}} and \eqn{T_{j}} are the surrogate and true endpoint values of subject \eqn{j}, and \eqn{Z_{j}} is the treatment indicator for subject \eqn{j}. 

To be in line with the Prentice criteria, Z should have a significant effect on S in model 1 (Prentice criterion 1), Z should have a significant effect on T in model 2 (Prentice criterion 2), S should have a significant effect on T in model 3 (Prentice criterion criterion 3), and the effect of Z on T should be fully captured by S in model 4 (Prentice criterion 4).  
 
The Prentice approach to assess surrogavy has some fundamental limitations. For example, the fourth Prentice criterion requires that the statistical test for the  \eqn{\beta_S} in model 4 is non-significant. This criterion is useful to reject a poor surrogate, but it is not suitable to validate a good surrogate (i.e., a non-significant result may always be attributable to a lack of statistical power). Even when lack of power would not be an issue, the result of the statistical test to evaluate the fourth Prentice criterion cannot prove that the effect of the treatment on the true endpoint is fully captured by the surrogate. 

The use of the Prentice approach to evaluate a surrogate is not recommended. Instead, consider using the single-trial meta-anlytic method (if no multiple clinical trials are available or if there is no other clustering unit in the data; see function \code{\link{Single.Trial.RE.AA}}) or the multiple-trial meta-analytic methods (see \code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}}, and \code{\link{BimixedContCont}}). 
 
}

\value{
\item{Prentice.Model.1}{An object of class \code{lm} that contains the fitted model 1 (using the Prentice approach).}
\item{Prentice.Model.2}{An object of class \code{lm} that contains the fitted model 2 (using the Prentice approach).}
\item{Prentice.Model.3}{An object of class \code{lm} that contains the fitted model 3 (using the Prentice approach).}
\item{Prentice.Model.4}{An object of class \code{lm} that contains the fitted model 4 (using the Prentice approach).}
\item{Prentice.Passed}{\code{Logical}. If all four Prentice criteria are fulfilled, \code{Prentice.Passed=TRUE}. If at least one criterion is not fulfilled, \code{Prentice.Passed=FALSE}.}
}

\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Prentice, R. L. (1989). Surrogate endpoints in clinical trials: definitions and
operational criteria. \emph{Statistics in Medicine, 8,} 431-440.

}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\examples{

## Load the ARMD dataset
data(ARMD)

## Evaluate the Prentice criteria in the ARMD dataset 
Prent <- Prentice(Dataset=ARMD, Surr=Diff24, True=Diff52, Treat=Treat, Pat.ID=Id)

# Summary of results
summary(Prent)
}
\keyword{Single-trial setting}
\keyword{Prentice criteria}