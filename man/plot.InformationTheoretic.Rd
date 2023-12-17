\name{plot Information-Theoretic}
\alias{plot Information-Theoretic}
\alias{plot.FixedContContIT}
\alias{plot.MixedContContIT}
\title{Provides plots of trial- and individual-level surrogacy in the Information-Theoretic framework}
\description{
Produces plots that provide a graphical representation of trial- and/or individual-level surrogacy (R2_ht and R2_h) based on the Information-Theoretic approach of Alonso & Molenberghs (2007).   
}
\usage{
\method{plot}{FixedContContIT}(x, Trial.Level=TRUE, Weighted=TRUE, Indiv.Level=TRUE, 
Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, Main.Trial, Main.Indiv, 
Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), \dots)

\method{plot}{MixedContContIT}(x, Trial.Level=TRUE, Weighted=TRUE, Indiv.Level=TRUE, 
Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, Main.Trial, Main.Indiv, 
Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), \dots)

}
\arguments{
  \item{x}{An object of class \code{MixedContContIT} or \code{FixedContContIT}.}
  \item{Trial.Level}{Logical. If \code{Trial.Level=TRUE}, a plot of the trial-specific treatment effects on the true endpoint against the trial-specific treatment effect on the surrogate endpoints is provided (as a graphical representation of \eqn{R_{ht}}). Default \code{TRUE}.}
  \item{Weighted}{Logical. This argument only has effect when the user requests a trial-level surrogacy plot (i.e., when \code{Trial.Level=TRUE} in the function call). If \code{Weighted=TRUE}, the circles that depict the trial-specific treatment effects on the true endpoint against the surrogate endpoint are proportional to the number of patients in the trial. If \code{Weighted=FALSE}, all circles have the same size. Default \code{TRUE}.}
  \item{Indiv.Level}{Logical. If \code{Indiv.Level=TRUE}, a plot of the trial- and treatment-corrected residuals of the true and surrogate endpoints is provided. This plot provides a graphical representation of \eqn{R_{h}}. Default \code{TRUE}.}
  \item{Xlab.Indiv}{The legend of the X-axis of the plot that depicts individual-level surrogacy. Default "Residuals for the surrogate endpoint (\eqn{\varepsilon_{Sij}})".}
  \item{Ylab.Indiv}{The legend of the Y-axis of the plot that depicts individual-level surrogacy. Default "Residuals for the true endpoint (\eqn{\varepsilon_{Tij}})".}
  \item{Xlab.Trial}{The legend of the X-axis of the plot that depicts trial-level surrogacy. Default "Treatment effect on the surrogate endpoint (\eqn{\alpha_{i}})".}
  \item{Ylab.Trial}{The legend of the Y-axis of the plot that depicts trial-level surrogacy. Default "Treatment effect on the true endpoint (\eqn{\beta_{i}})".}
  \item{Main.Indiv}{The title of the plot that depicts individual-level surrogacy. Default "Individual-level surrogacy".}
  \item{Main.Trial}{The title of the plot that depicts trial-level surrogacy. Default "Trial-level surrogacy".}
    \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{plot()}.}
}
\references{
Alonso, A, & Molenberghs, G. (2007). Surrogate marker evaluation from an information theory perspective. \emph{Biometrics, 63,} 180-186.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{MixedContContIT}, \link{FixedContContIT}
}


\keyword{Plot surrogacy}
\keyword{Information-Theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Multiple-trial setting}