\name{plot TwoStageSurvSurv}
\alias{plot.TwoStageSurvSurv}
\title{Plots trial-level surrogacy in the meta-analytic framework when two survival endpoints are considered.
}
\description{
Produces a plot that graphically depicts trial-level surrogacy when the surrogate and true endpoints are survival endpoints.    
}
\usage{
\method{plot}{TwoStageSurvSurv}(x, Weighted=TRUE, xlab, ylab, main,
Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), \dots)
}

\arguments{
  \item{x}{An object of class \code{TwoStageContCont}.}
  \item{Weighted}{Logical. If \code{Weighted=TRUE}, the circles that depict the trial-specific treatment effects on the true endpoint against the surrogate endpoint are proportional to the number of patients in the trial. If \code{Weighted=FALSE}, all circles have the same size. Default \code{TRUE}.}
  \item{xlab}{The legend of the X-axis, default "Treatment effect on the surrogate endpoint (\eqn{\alpha_{i}})".}
  \item{ylab}{The legend of the Y-axis, default "Treatment effect on the true endpoint (\eqn{\beta_{i}})".}
  \item{main}{The title of the plot, default "Trial-level surrogacy".}
  \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{plot()}.}
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{TwoStageSurvSurv}
}
\examples{
# Open Ovarian dataset
data(Ovarian)
# Conduct analysis
Results <- TwoStageSurvSurv(Dataset = Ovarian, Surr = Pfs, SurrCens = PfsInd,
True = Surv, TrueCens = SurvInd, Treat = Treat, Trial.ID = Center)
# Examine results of analysis
summary(Results)
plot(Results)
}
\keyword{Plot surrogacy}
\keyword{Information-theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Multiple-trial setting}
\keyword{Single-trial setting}
\keyword{Survival endpoints}