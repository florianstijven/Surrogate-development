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
\examples{
## Load ARMD dataset
data(ARMD)

## Conduct a surrogacy analysis, using a weighted reduced univariate fixed effect model:
Sur <- MixedContContIT(Dataset=ARMD, Surr=Diff24, True=Diff52, Treat=Treat, Trial.ID=Center, 
Pat.ID=Id, Model=c("Full"))

## Request both trial- and individual-level surrogacy plots. In the trial-level plot,
## make the size of the circles proportional to the number of patients in a trial:
\dontrun{plot(Sur, Trial.Level=TRUE, Weighted=TRUE, Indiv.Level=TRUE)}

## Make a trial-level surrogacy plot using filled blue circles that 
## are transparent (to make sure that the results of overlapping trials remain
## visible), and modify the title and the axes labels of the plot: 
\dontrun{plot(Sur, pch=16, col=rgb(.3, .2, 1, 0.3), Indiv.Level=FALSE, Trial.Level=TRUE, 
Weighted=TRUE, Main.Trial=c("Trial-level surrogacy (ARMD dataset)"), 
Xlab.Trial=c("Difference in vision after 6 months (Surrogate)"),
Ylab.Trial=c("Difference in vision after 12 months (True enpoint)"))}

## Add the estimated R2_ht value in the previous plot at position (X=-2.2, Y=0)  
## (the previous plot should not have been closed):
\dontrun{R2ht <- format(round(as.numeric(Sur$R2ht[1]), 3))
text(x=-2.2, y=0, cex=1.4, labels=(bquote(paste("R"[ht]^{2}, "="~.(R2ht)))))}

## Make an Individual-level surrogacy plot with red squares to depict individuals
## (rather than black circles):
\dontrun{plot(Sur, pch=15, col="red", Indiv.Level=TRUE, Trial.Level=FALSE)}

## Same plot as before, but now with smaller squares, a y-axis with range [-40; 40], 
## and the estimated R2_h value in the title of the plot:
R2h <- format(round(as.numeric(Sur$R2h[1]), 3))
\dontrun{plot(Sur, pch=15, col="red", Indiv.Level=TRUE, Trial.Level=FALSE, cex=.5, 
ylim=c(-40, 40), Main.Indiv=bquote(paste("R"[h]^{2}, "="~.(R2h))))}


}

\keyword{Plot surrogacy}
\keyword{Information-Theoretic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Multiple-trial setting}