\name{plot Meta-Analytic}
\alias{plot Meta-Analytic}
\alias{plot.BifixedContCont}
\alias{plot.BimixedContCont}
\alias{plot.UnifixedContCont}
\alias{plot.UnimixedContCont}
\title{Provides plots of trial- and individual-level surrogacy in the meta-analytic framework
}
\description{
Produces plots that provide a graphical representation of trial- and/or individual-level surrogacy based on the meta-analytic approach of Buyse & Molenberghs (2000) in the single- and multiple-trial settings.   
}
\usage{
\method{plot}{BifixedContCont}(x, Trial.Level=TRUE, Weighted=TRUE, 
Indiv.Level=TRUE, Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, 
Main.Trial, Main.Indiv, Par=par(oma=c(0, 0, 0, 0), 
mar=c(5.1, 4.1, 4.1, 2.1)), \dots)

\method{plot}{BimixedContCont}(x, Trial.Level=TRUE, Weighted=TRUE, 
Indiv.Level=TRUE, Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, 
Main.Trial, Main.Indiv, Par=par(oma=c(0, 0, 0, 0), 
mar=c(5.1, 4.1, 4.1, 2.1)), \dots)

\method{plot}{UnifixedContCont}(x, Trial.Level=TRUE, Weighted=TRUE, 
Indiv.Level=TRUE, Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, 
Main.Trial, Main.Indiv, Par=par(oma=c(0, 0, 0, 0), 
mar=c(5.1, 4.1, 4.1, 2.1)), \dots)

\method{plot}{UnimixedContCont}(x, Trial.Level=TRUE, Weighted=TRUE, 
Indiv.Level=TRUE, Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, 
Main.Trial, Main.Indiv, Par=par(oma=c(0, 0, 0, 0), 
mar=c(5.1, 4.1, 4.1, 2.1)), \dots)

}
\arguments{
  \item{x}{An object of class \code{UnifixedContCont}, \code{BifixedContCont}, \code{UnimixedContCont}, \code{BimixedContCont}, or \code{Single.Trial.RE.AA}.}
  \item{Trial.Level}{Logical. If \code{Trial.Level=TRUE} and an object of class \code{UnifixedContCont}, \code{BifixedContCont}, \code{UnimixedContCont},  or \code{BimixedContCont} is considered, a plot of the trial-specific treatment effects on the true endpoint against the trial-specific treatment effect on the surrogate endpoints is provided (as a graphical representation of \eqn{R_{trial}}). If \code{Trial.Level=TRUE} and an object of class \code{Single.Trial.RE.AA} is considered, a plot of the treatment effect on the true endpoint against the treatment effect on the surrogate endpoint is provided, and a regression line that goes through the origin with slope RE is added to the plot (to depict the constant RE assumption, see \code{\link{Single.Trial.RE.AA}} for details).  If \code{Trial.Level=FALSE}, this plot is not provided. Default \code{TRUE}.}
  \item{Weighted}{Logical. This argument only has effect when the user requests a trial-level surrogacy plot (i.e., when \code{Trial.Level=TRUE} in the function call) and when an object of class \code{UnifixedContCont}, \code{BifixedContCont}, \code{UnimixedContCont}, or \code{BimixedContCont} is considered (not when an object of class \code{Single.Trial.RE.AA} is considered). If \code{Weighted=TRUE}, the circles that depict the trial-specific treatment effects on the true endpoint against the surrogate endpoint are proportional to the number of patients in the trial. If \code{Weighted=FALSE}, all circles have the same size. Default \code{TRUE}.}
  \item{Indiv.Level}{Logical. If \code{Indiv.Level=TRUE}, a plot of the trial- and treatment-corrected residuals of the true and surrogate endpoints is provided (when an object of class \code{UnifixedContCont}, \code{BifixedContCont}, \code{UnimixedContCont}, or \code{BimixedContCont} is considered), or a plot of the treatment-corrected residuals (when an object of class \code{Single.Trial.RE.AA} is considered). This plot provides a graphical representation of \eqn{R_{indiv}}. If \code{Indiv.Level=FALSE}, this plot is not provided. Default \code{TRUE}.}
  \item{Xlab.Indiv}{The legend of the X-axis of the plot that depicts individual-level surrogacy. Default "Residuals for the surrogate endpoint (\eqn{\varepsilon_{Sij}})" (without the \eqn{i} subscript when an object of class \code{Single.Trial.RE.AA} is considered).}
  \item{Ylab.Indiv}{The legend of the Y-axis of the plot that depicts individual-level surrogacy. Default "Residuals for the true endpoint (\eqn{\varepsilon_{Tij}})" (without the \eqn{i} subscript when an object of class \code{Single.Trial.RE.AA} is considered).}
  \item{Xlab.Trial}{The legend of the X-axis of the plot that depicts trial-level surrogacy. Default "Treatment effect on the surrogate endpoint (\eqn{\alpha_{i}})" (without the \eqn{i} subscript when an object of class \code{Single.Trial.RE.AA} is considered).}
  \item{Ylab.Trial}{The legend of the Y-axis of the plot that depicts trial-level surrogacy. Default "Treatment effect on the true endpoint (\eqn{\beta_{i}})" (without the \eqn{i} subscript when an object of class \code{Single.Trial.RE.AA} is considered).}
  \item{Main.Indiv}{The title of the plot that depicts individual-level surrogacy. Default "Individual-level surrogacy" when an object of class \code{UnifixedContCont}, \code{BifixedContCont}, \code{UnimixedContCont}, or \code{BimixedContCont} is considered, and "Adjusted Association (\eqn{rho_{Z}}) when an object of class \code{Single.Trial.RE.AA} is considered.}
  \item{Main.Trial}{The title of the plot that depicts trial-level surrogacy. Default "Trial-level surrogacy" (when an object of class \code{UnifixedContCont}, \code{BifixedContCont}, \code{UnimixedContCont},  or \code{BimixedContCont} is considered) or "Relative Effect (RE)" (when an object of class \code{Single.Trial.RE.AA} is considered).}
  \item{Par}{Graphical parameters for the plot. Default \code{par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1))}.}
  \item{\dots}{Extra graphical parameters to be passed to \code{plot()}.}
}
\references{
Buyse, M., Molenberghs, G., Burzykowski, T., Renard, D., & Geys, H. (2000). The validation of surrogate endpoints in meta-analysis of randomized experiments. \emph{Biostatistics, 1,} 49-67.
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{UnifixedContCont}, \link{BifixedContCont}, \link{UnifixedContCont}, \link{BimixedContCont}, \link{Single.Trial.RE.AA}
}
\examples{
##### Multiple-trial setting

## Load ARMD dataset
data(ARMD)

## Conduct a surrogacy analysis, using a weighted reduced univariate fixed effect model:
Sur <- UnifixedContCont(Dataset=ARMD, Surr=Diff24, True=Diff52, Treat=Treat, Trial.ID=Center, 
Pat.ID=Id, Number.Bootstraps=100, Model=c("Reduced"), Weighted=TRUE)

## Request both trial- and individual-level surrogacy plots. In the trial-level plot,
## make the size of the circles proportional to the number of patients in a trial:
plot(Sur, Trial.Level=TRUE, Weighted=TRUE, Indiv.Level=TRUE)

## Make a trial-level surrogacy plot using filled blue circles that 
## are transparent (to make sure that the results of overlapping trials remain
## visible), and modify the title and the axes labels of the plot: 
plot(Sur, pch=16, col=rgb(.3, .2, 1, 0.3), Indiv.Level=FALSE, Trial.Level=TRUE, 
Weighted=TRUE, Main.Trial=c("Trial-level surrogacy (ARMD dataset)"), 
Xlab.Trial=c("Difference in vision after 6 months (Surrogate)"),
Ylab.Trial=c("Difference in vision after 12 months (True enpoint)"))

## Add the estimated R2_trial value in the previous plot at position (X=-7, Y=11)  
## (the previous plot should not have been closed):
R2trial <- format(round(as.numeric(Sur$Trial.R2[1]), 3))
text(x=-7, y=11, cex=1.4, labels=(bquote(paste("R"[trial]^{2}, "="~.(R2trial)))))

## Make an Individual-level surrogacy plot with red squares to depict individuals
## (rather than black circles):
plot(Sur, pch=15, col="red", Indiv.Level=TRUE, Trial.Level=FALSE)

## Same plot as before, but now with smaller squares, a y-axis with range [-40; 40], 
## and the estimated R2_indiv value in the title of the plot:
R2ind <- format(round(as.numeric(Sur$Indiv.R2[1]), 3))
plot(Sur, pch=15, col="red", Indiv.Level=TRUE, Trial.Level=FALSE, cex=.5, 
ylim=c(-40, 40), Main.Indiv=bquote(paste("R"[indiv]^{2}, "="~.(R2ind))))


##### Single-trial setting

## Conduct a surrogacy analysis in the single-trial meta-analytic setting:
SurSTS <- Single.Trial.RE.AA(Dataset=ARMD, Surr=Diff24, True=Diff52, Treat=Treat, Pat.ID=Id)

# Request a plot of individual-level surrogacy and a plot that depicts the Relative effect 
# and the constant RE assumption:
plot(SurSTS, Trial.Level=TRUE, Indiv.Level=TRUE)
}

\keyword{Plot surrogacy}
\keyword{Meta-analytic framework}
\keyword{Trial-level surrogacy}
\keyword{Individual-level surrogacy}
\keyword{Multiple-trial setting}
\keyword{Single-trial setting}
\keyword{Continuous endpoint}