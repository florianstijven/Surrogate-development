\name{plot FixedDiscrDiscrIT}
\alias{plot FixedDiscrDiscrIT}
\alias{plot.FixedDiscrDiscrIT}
\title{Provides plots of trial-level surrogacy in the Information-Theoretic framework}
\description{
Produces plots that provide a graphical representation of trial level surrogacy \eqn{R^2_{ht}} based on the Information-Theoretic approach of Alonso & Molenberghs (2007).   
}
\usage{
 \method{plot}{FixedDiscrDiscrIT}(x, Weighted=TRUE, Xlab.Trial, Ylab.Trial, Main.Trial,
	 Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), \dots)
}
\arguments{
    \item{x}{An object of class \code{FixedDiscrDiscrIT}.}
    \item{Weighted}{Logical. This argument only has effect when the user requests a trial-level surrogacy plot (i.e., when \code{Trial.Level=TRUE} in the function call). If \code{Weighted=TRUE}, the circles that depict the trial-specific treatment effects on the true endpoint against the surrogate endpoint are proportional to the number of patients in the trial. If \code{Weighted=FALSE}, all circles have the same size. Default \code{TRUE}.}
    \item{Xlab.Trial}{The legend of the X-axis of the plot that depicts trial-level surrogacy. Default "Treatment effect on the surrogate endpoint (\eqn{\alpha_{i}})".}
    \item{Ylab.Trial}{The legend of the Y-axis of the plot that depicts trial-level surrogacy. Default "Treatment effect on the true endpoint (\eqn{\beta_{i}})".}
    \item{Main.Trial}{The title of the plot that depicts trial-level surrogacy. Default "Trial-level surrogacy".}
      \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), 
mar=c(5.1, 4.1, 4.1, 2.1))}.}
    \item{\dots}{Extra graphical parameters to be passed to \code{plot()}.}
  }
 \references{
  Alonso, A, & Molenberghs, G. (2007). Surrogate marker evaluation from an information theory perspective. \emph{Biometrics, 63,} 180-186.
  }
  \author{
Hannah M. Ensor & Christopher J. Weir
  }
  \seealso{\link{FixedDiscrDiscrIT}
  }
\examples{
\dontrun{ # Time consuming (>5sec) code part
# Simulate the data:
Sim.Data.MTS(N.Total=2000, N.Trial=100, R.Trial.Target=.8, R.Indiv.Target=.8,
             Seed=123, Model="Full")
             
# create a binary true and ordinal surrogate outcome
Data.Observed.MTS$True<-findInterval(Data.Observed.MTS$True, 
        c(quantile(Data.Observed.MTS$True,0.5)))
Data.Observed.MTS$Surr<-findInterval(Data.Observed.MTS$Surr, 
        c(quantile(Data.Observed.MTS$Surr,0.333),quantile(Data.Observed.MTS$Surr,0.666)))

# Assess surrogacy based on a full fixed-effect model
# in the information-theoretic framework for a binary surrogate and ordinal true outcome:
SurEval <- FixedDiscrDiscrIT(Dataset=Data.Observed.MTS, Surr=Surr, True=True, Treat=Treat,
Trial.ID=Trial.ID, Setting="ordbin")

## Request trial-level surrogacy plot. In the trial-level plot,
## make the size of the circles proportional to the number of patients in a trial:
plot(SurEval, Weighted=FALSE)

}}
\keyword{Plot surrogacy}
\keyword{Information-Theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Multiple-trial setting}