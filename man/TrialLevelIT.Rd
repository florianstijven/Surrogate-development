\name{TrialLevelIT}
\alias{TrialLevelIT}
\title{Estimates trial-level surrogacy in the information-theoretic framework}
\description{
The function \code{TrialLevelIT} estimates trial-level surrogacy based on the vectors of treatment effects on \eqn{S} (i.e., \eqn{\alpha_{i}}), intercepts on \eqn{S} (i.e., \eqn{\mu_{i}}) and \eqn{T} (i.e., \eqn{\beta_{i}}) in the different trials. See the \bold{Details} section below.  
}
\usage{
TrialLevelIT(Alpha.Vector, Mu_S.Vector=NULL, 
Beta.Vector, N.Trial, Model="Reduced", Alpha=.05)
}
\arguments{
\item{Alpha.Vector}{The vector of treatment effects on \eqn{S} in the different trials, i.e., \eqn{\alpha_{i}}.}
\item{Mu_S.Vector}{The vector of intercepts for \eqn{S} in the different trials, i.e., \eqn{\mu_{Si}}. Only required when a full model is requested.}
\item{Beta.Vector}{The vector of treatment effects on \eqn{T} in the different trials, i.e., \eqn{\beta_{i}}.}
\item{N.Trial}{The total number of available trials.}
\item{Model}{The type of model that should be fitted, i.e., \code{Model=c("Full")} or \code{Model=c("Reduced")}. See the \bold{Details} section below. Default \code{Model=c("Reduced")}.}
  \item{Alpha}{The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}} and \eqn{R_{trial}}. Default \eqn{0.05}.}
}
\details{
When a full model is requested (by using the argument \code{Model=c("Full")} in the function call), trial-level surrogacy is assessed by fitting the following univariate model: 

\deqn{{\beta}_{i}=\lambda_{0}+\lambda_{1}{\mu_{Si}}+\lambda_{2}{\alpha}_{i}+ \varepsilon_{i},   (1)}

where \eqn{\beta_i} = the trial-specific treatment effects on \eqn{T}, \eqn{\mu_{Si}} =  the trial-specific intercepts for \eqn{S}, and \eqn{\alpha_i} =  the trial-specific treatment effects on \eqn{S}. The \eqn{-2} log likelihood value of model (1) (\eqn{L_1}) is subsequently compared to the \eqn{-2} log likelihood value of an intercept-only model (\eqn{{\beta}_{i}=\lambda_{3}}; \eqn{L_0}), and \eqn{R^2_{ht}} is computed based based on the Variance Reduction Factor (for details, see Alonso & Molenberghs, 2007):

\deqn{R^2_{ht}= 1 - exp \left(-\frac{L_1-L_0}{N} \right),}

where \eqn{N} is the number of trials. 

When a reduced model is requested (by using the argument \code{Model=c("Reduced")} in the function call), the following model is fitted:

\deqn{{\beta}_{i}=\lambda_{0}+\lambda_{1}{\alpha}_{i}+\varepsilon_{i}.}

The \eqn{-2} log likelihood value of this model (\eqn{L_1} for the reduced model) is subsequently compared to the \eqn{-2} log likelihood value of an intercept-only model (\eqn{{\beta}_{i}=\lambda_{3}}; \eqn{L_0}), and \eqn{R^2_{ht}} is computed based on the reduction in the likelihood (as described above).
 
}
\value{
An object of class \code{TrialLevelIT} with components,
\item{Alpha.Vector}{The vector of treatment effects on \eqn{S} in the different trials.}
\item{Beta.Vector}{The vector of treatment effects on \eqn{T} in the different trials.}
\item{N.Trial}{The total number of trials.}
\item{R2.ht}{A \code{data.frame} that contains the trial-level coefficient of determination (\eqn{R^2_{ht}}), its standard error and confidence interval.}
}
\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.

Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{UnimixedContCont}}, \code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{BimixedContCont}}, \code{\link{plot.TrialLevelIT}}
}
\examples{
# Generate vector treatment effects on S
set.seed(seed = 1)
Alpha.Vector <- seq(from = 5, to = 10, by=.1) + runif(min = -.5, max = .5, n = 51)

# Generate vector treatment effects on T
set.seed(seed=2)
Beta.Vector <- (Alpha.Vector * 3) + runif(min = -5, max = 5, n = 51)

# Apply the function to estimate R^2_{h.t}
Fit <- TrialLevelIT(Alpha.Vector=Alpha.Vector,
Beta.Vector=Beta.Vector, N.Trial=50, Model="Reduced")

summary(Fit)
plot(Fit)
}
\keyword{Multiple-trial setting}
\keyword{Information-theoretic framework}
\keyword{Trial-level surrogacy}