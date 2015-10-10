\name{Pred.TrialT.ContCont}
\alias{Pred.TrialT.ContCont}
\title{
Compute the expected treatment effect on the true endpoint in a new trial (when both S and T are normally distributed continuous endpoints)
}
\description{
The key motivation to evaluate a surrogate endpoint is to be able to predict the  treatment effect on the true endpoint \eqn{T} based on the treatment effect on \eqn{S} in a new trial \eqn{i=0}. The function \code{Pred.TrialT.ContCont} allows for making such predictions based on fitted models of class \code{\link{BimixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}} and \code{\link{UnifixedContCont}}.  
}
\usage{
Pred.TrialT.ContCont(Object, mu_S0, alpha_0, alpha.CI=0.05)
}
\arguments{
  \item{Object}{A fitted object of class \code{\link{BimixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}} and \code{\link{UnifixedContCont}}. Some of the components in these fitted objects are needed to estimate \eqn{E(\beta + b_0)} and its variance.}
  \item{mu_S0}{The intercept of a regression model in the new trial \eqn{i=0} where the surrogate endpoint is regressed on the true endpoint, i.e., \eqn{S_[0j]=\mu_{S0} + \alpha_0 Z_{0j} + \varepsilon_{S0j}}, where \eqn{S} is the surrogate endpoint, \eqn{j} is the patient indicator, and \eqn{Z} is the treatment. This argument only needs to be specified when a full model was used to examine surroacy.}
  \item{alpha_0}{The regression weight of the treatment in the regression model specified under argument \code{mu_S0}.}
  \item{alpha.CI}{The \eqn{\alpha}-level to be used to determine the confidence interval around \eqn{E(\beta + b_0)}. Default \code{alpha.CI=0.05}.}
}
\details{
The key motivation to evaluate a surrogate endpoint is to be able to predict the  treatment effect on the true endpoint \eqn{T} based on the treatment effect on \eqn{S} in a new trial \eqn{i=0}.

When a so-called full (fixed or mixed) bi- or univariate model was fitted in the surrogate evaluation phase (for details, see \code{\link{BimixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}} and \code{\link{UnifixedContCont}}), this prediction is made as:

\deqn{E(\beta + b_0 | m_{S0}, a_0) = \beta + \left(\begin{array}{c}
d_{Sb}\\
d_{ab}
\end{array}\right)^T \left(\begin{array}{cc}
d_{SS} & D_{Sa}\\
d_{Sa} & d_{aa}
\end{array}\right)^{-1} \left(\begin{array}{c}
\mu_{S0} - \mu_S\\
\alpha_0 - \alpha
\end{array}\right) }

\deqn{Var(\beta + b_0 | m_{S0}, a_0) = d_{bb} + \left(\begin{array}{c}
d_{Sb}\\
d_{ab}
\end{array}\right)^T \left(\begin{array}{cc}
d_{SS} & D_{Sa}\\
d_{Sa} & d_{aa}
\end{array}\right)^{-1} \left(\begin{array}{c}
d_{Sb}\\
d_{ab}
\end{array}\right),}
 
where all components are defined as in \code{\link{BimixedContCont}}. When the univariate mixed-effects models are used or the (univariate or bivariate) fixed effects models, the fitted components contained in \code{D.Equiv} are used instead of those in \code{D}.


When a reduced-model approach was used in the surrogate evaluation phase, the prediction is made as:

\deqn{E(\beta + b_0 | a_0) = \beta + \frac{}d_{ab}{d_{aa}}+ (\alpha_0 - \alpha),}

\deqn{Var(\beta + b_0 | a_0) = d_{bb} - \frac{d_{ab}^2}{d_{aa}},}
 
where all components are defined as in \code{\link{BimixedContCont}}. When the univariate mixed-effects models are used or the (univariate or bivariate) fixed effects models, the fitted components contained in \code{D.Equiv} are used instead of those in \code{D}.

}
\value{
\item{Beta_0}{The predicted \eqn{\beta_0}.}
\item{Variance}{The variance of the prediction.}
\item{Lower}{The lower bound of the confidence interval around the expected \eqn{\beta_0}.}
\item{Upper}{The upper bound of the confidence interval around the expected \eqn{\beta_0}.}
\item{alpha.CI}{The \eqn{\alpha}-level used to establish the confidence interval.}
\item{Surr.Model}{The model that was used to compute \eqn{\beta_0}}.
\item{alpha_0}{The slope of the regression model specified in the \code{Arguments} section.}
}
\references{
Burzykowski, T., Molenberghs, G., & Buyse, M. (2005). \emph{The evaluation of surrogate endpoints}. New York: Springer-Verlag.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{
\code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}}
}
\examples{
# Generate dataset
Sim.Data.MTS(N.Total=2000, N.Trial=15, R.Trial.Target=.8, 
R.Indiv.Target=.8, D.aa=10, D.bb=50, Fixed.Effects=c(1, 2, 30, 90), 
Seed=1)

# Evaluate surrogacy using a reduced bivariate mixed-effects model
BimixedFit <- BimixedContCont(Dataset = Data.Observed.MTS, Surr = Surr, 
True = True, Treat = Treat, Trial.ID = Trial.ID, Pat.ID = Pat.ID, 
Model="Reduced")

# Suppose that in a new trial, it was estimated alpha_0 = 30
# predict beta_0 in this trial
Pred_Beta <- Pred.TrialT.ContCont(Object = BimixedFit, 
alpha_0 = 30)

# Examine the results
summary(Pred_Beta)

# Plot the results
plot(Pred_Beta)
}
\keyword{New trial}
\keyword{Predict treatment effect T}
