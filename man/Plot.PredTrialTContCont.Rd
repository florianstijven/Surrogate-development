\name{plot PredTrialTContCont}
\alias{plot PredTrialTContCont}
\alias{plot.PredTrialTContCont}
\title{
Plots the expected treatment effect on the true endpoint in a new trial (when both S and T are normally distributed continuous endpoints)
}
\description{
The key motivation to evaluate a surrogate endpoint is to be able to predict the  treatment effect on the true endpoint \eqn{T} based on the treatment effect on \eqn{S} in a new trial \eqn{i=0}. The function \code{Pred.TrialT.ContCont} allows for making such predictions. The present plot function shows the results graphically.  
}
\usage{
\method{plot}{PredTrialTContCont}(x, Size.New.Trial=5, CI.Segment=1, ...)
}
\arguments{
  \item{x}{A fitted object of class \code{Pred.TrialT.ContCont}, for details see \code{\link{Pred.TrialT.ContCont}}.}
  \item{Size.New.Trial}{The expected treatment effect on \eqn{T} is drawn as a black circle with size specified by \code{Size.New.Trial}. Default \code{Size.New.Trial=5}.}
  \item{CI.Segment}{The confidence interval around the expected treatment effect on \eqn{T} is depicted by a dashed horizontal line. By default, the width of the horizontal line of the horizontal section of the confidence interval indicator is \eqn{2} times the values specified by \code{CI.Segment}. Default \eqn{CI.Segment=1}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{plot()}.}
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{Pred.TrialT.ContCont}
}
\examples{
# Generate dataset
Sim.Data.MTS(N.Total=2000, N.Trial=15, R.Trial.Target=.95, 
R.Indiv.Target=.8, D.aa=10, D.bb=50, 
Fixed.Effects=c(1, 2, 30, 90), Seed=1)

# Evaluate surrogacy using a reduced bivariate mixed-effects model
BimixedFit <- BimixedContCont(Dataset = Data.Observed.MTS, 
Surr = Surr, True = True, Treat = Treat, Trial.ID = Trial.ID, 
Pat.ID = Pat.ID, Model="Reduced")

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
