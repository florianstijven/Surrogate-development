\name{plot SPF BinBin}
\alias{plot SPF BinBin}
\alias{plot.SPF.BinBin}
\title{
Plots the surrogate predictive function (SPF).
}
\description{
Plots the surrogate predictive function (SPF), i.e., \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)}, in the setting where both \eqn{S} and \eqn{T} are binary endpoints. For example, \eqn{r(-1,1)} quantifies the probability that the treatment has a negative effect on the true endpoint (\eqn{\Delta T=-1}) given that it has a positive effect on the surrogate (\eqn{\Delta S=1}).}

\usage{
\method{plot}{SPF.BinBin}(x, Type="All.Histograms", Specific.Pi="r_0_0", Col="grey", 
Box.Plot.Outliers=FALSE, Legend.Pos="topleft", Legend.Cex=1, ...)
}
\arguments{
\item{x}{A fitted object of class \code{SPF.BinBin}. See \code{\link{ICA.BinBin}}.}
\item{Type}{The type of plot that is requested. Possible choices are: \code{Type="All.Histograms"}, the histograms of all \eqn{9}  \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors arranged in a \eqn{3} by \eqn{3} grid; \code{Type="All.Densities"}, plots of densities of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors; \code{Type="Histogram"}, the histogram of a particular \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vector (the \code{Specific.Pi=} argument has to be used to specify the desired \eqn{r(i,j)}); \code{Type="Density"}, the density of a particular \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vector (the \code{Specific.Pi=} argument has to be used to specify the desired \eqn{r(i,j)}); \code{Type="Box.Plot"}, a box plot of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors; \code{Type="Lines.Mean"}, a line plot the depicts the means of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors; \code{Type="Lines.Median"}, a line plot the depicts the medians of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors; \code{Type="Lines.Mode"}, a line plot the depicts the modes of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors; \code{Type="3D.Mean"}, a 3D bar plot the depicts the means of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors; \code{Type="3D.Median"}, a 3D bar plot the depicts the medians of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors; \code{Type="3D.Mode"}, a 3D bar plot the depicts the modes of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors; \code{Type="3D.Spinning.Mean"}, a spinning 3D plot that depicts the means of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors that can be rotated; \code{Type="3D.Spinning.Median"}, a spinning 3D plot that depicts the medians of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors that can be rotated; \code{Type="3D.Spinning.Mode"}, a spinning 3D plot that depicts the modes of all \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vectors that can be rotated.}
\item{Specific.Pi}{When \code{Type="Histogram"} or \code{Type="Density"} , the histogram/density of a particular \eqn{r(i,j)=P(\Delta T=i|\Delta S=j)} vector is shown. The \code{Specific.Pi=} argument is used to specify the desired \eqn{r(i,j)}). Default \code{r_0_0}.}
\item{Col}{The color of the bins or lines when histograms or density plots are requested. Default \code{"grey"}.}
\item{Box.Plot.Outliers}{Logical. Should outliers be depicted in the box plots?. Default \code{FALSE}.}
\item{Legend.Pos}{Position of the legend when a \code{type="Box.Plot"}, \code{type="Lines.Mean"}, \code{type="Lines.Median"}, or \code{type="Lines.Mode"} is requested. Default \code{"topleft"}.}
\item{Legend.Cex}{Size of the legend when a \code{type="Box.Plot"}, \code{type="Lines.Mean"}, \code{type="Lines.Median"}, or \code{type="Lines.Mode"} is requested. Default \code{1}.}
\item{\dots}{Arguments to be passed to the plot, histogram, ... functions.}
}
\references{
Alonso, A., Van der Elst, W., & Molenberghs, G. (2015). Assessing a surrogate effect predictive value in a causal inference framework.}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\seealso{\link{SPF.BinBin}}
\examples{
\dontrun{
# Generate plausible values for Pi  
ICA <- ICA.BinBin.Grid.Sample(pi1_1_=0.341, pi0_1_=0.119,
pi1_0_=0.254, pi_1_1=0.686, pi_1_0=0.088, pi_0_1=0.078, Seed=1,
Monotonicity=c("General"), M=5000)
           
# Compute the surrogate predictive function (SPF)
SPF <- SPF.BinBin(ICA)

# Explore the results
summary(SPF)

# Examples of plots 
plot(SPF, Type="All.Histograms")
plot(SPF, Type="All.Densities")
plot(SPF, Type="Histogram", Specific.Pi="r_0_0")
plot(SPF, Type="Box.Plot", Legend.Pos="topleft", Legend.Cex=.7)
plot(SPF, Type="Lines.Mean")
plot(SPF, Type="Lines.Median")
plot(SPF, Type="Lines.Mode")
plot(SPF, Type="3D.Mean")
plot(SPF, Type="3D.Median")
plot(SPF, Type="3D.Mode")
plot(SPF, Type="3D.Spinning.Mean")
plot(SPF, Type="3D.Spinning.Median")
plot(SPF, Type="3D.Spinning.Mode")
} 
}
\keyword{Plot SPF}
\keyword{Causal-Inference framework}
\keyword{Sensitivity}