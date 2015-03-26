\name{ICA.BinBin.CounterAssum}
\alias{ICA.BinBin.CounterAssum}
\title{ICA (binary-binary setting) that is obtaied when the counterfactual correlations are assumed to fall within some prespecified ranges.}
\description{
Shows the results of ICA (binary-binary setting) in the subgroup of results where the counterfactual correlations are assumed to fall within some prespecified ranges.
}
\usage{
ICA.BinBin.CounterAssum(x, r_h_S0S1_min, r_h_S0S1_max, r_h_S0T1_min, 
r_h_S0T1_max, r_h_T0T1_min, r_h_T0T1_max, r_h_T0S1_min, r_h_T0S1_max, 
Monotonicity="General", Type="Freq", MainPlot=" ", Cex.Legend=1, 
Cex.Position="topright", ...)
}
\arguments{
  \item{x}{An object of class \code{ICA.BinBin}. See \code{\link{ICA.BinBin}}.}
  \item{r_h_S0S1_min}{The minimum value to be considered for the counterfactual correlation \eqn{r_h(S0,S1)}.}
  \item{r_h_S0S1_max}{The maximum value to be considered for the counterfactual correlation \eqn{r_h(S0,S1)}.}
  \item{r_h_S0T1_min}{The minimum value to be considered for the counterfactual correlation \eqn{r_h(S0,T1)}.}
  \item{r_h_S0T1_max}{The maximum value to be considered for the counterfactual correlation \eqn{r_h(S0,T1)}.}
  \item{r_h_T0T1_min}{The minimum value to be considered for the counterfactual correlation \eqn{r_h(T0,T1)}.}
  \item{r_h_T0T1_max}{The maximum value to be considered for the counterfactual correlation \eqn{r_h(T0,T1)}.}
  \item{r_h_T0S1_min}{The minimum value to be considered for the counterfactual correlation \eqn{r_h(T0,S1)}.}
  \item{r_h_T0S1_max}{The maximum value to be considered for the counterfactual correlation \eqn{r_h(T0,S1)}.}
\item{Monotonicity}{Specifies whether the all results in the fitted object \code{ICA.BinBin} should be shown (i.e., \code{Monotonicity=c("General")}), or a subset of the results arising under specific assumptions (i.e., \code{Monotonicity=c("No")}, \code{Monotonicity=c("True.Endp")}, \code{Monotonicity=c("Surr.Endp")}, or \code{Monotonicity=c("Surr.True.Endp")}). Default \code{Monotonicity=c("General")}.}
\item{Type}{The type of plot that is produced. When \code{Type="Freq"} or \code{Type="Density"}, the Y-axis shows frequencies or densities of \eqn{R^2_{H}}.  When \code{Type="All.Densities"} and the fitted object of class \code{ICA.BinBin} was obtained using a general analysis (i.e., conducting the analyses assuming no monotonicity, monotonicity for \eqn{S} alone, monotonicity for \eqn{T} alone, and for both \eqn{S} and \eqn{T}, so using \code{Monotonicity=c("General")} in the function call of \code{ICA.BinBin}), the density plots are shown for the four scenarios where different assumptions regarding monotonicity are made. Default \code{"Freq"}.}  
\item{MainPlot}{The title of the plot. Default \code{" "}.}
\item{Cex.Legend}{The size of the legend when \code{Type="All.Densities"} is used. Default \code{Cex.Legend=1}.}
\item{Cex.Position}{The position of the legend, \code{Cex.Position="topright"} or \code{Cex.Position="topleft"}. Default \code{Cex.Position="topright"}.}
\item{...}{Other arguments to be passed to the \code{plot()} function.}
}

\references{
Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., & Burzykowski, T. (submitted). On the relationship between the causal inference and meta-analytic paradigms for the validation of surrogate markers. 

Van der Elst, W., Alonso, A., & Molenberghs, G. (submitted). An exploration of the relationship between causal inference and meta-analytic measures of surrogacy.}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{\link{ICA.BinBin}}
\examples{
\dontrun{
# Compute R2_H given the marginals specified as the pi's, making no 
# assumptions regarding monotonicity (general case)
ICA <- ICA.BinBin(pi1_1_=0.2619048, pi1_0_=0.2857143, pi_1_1=0.6372549, 
pi_1_0=0.07843137, pi0_1_=0.1349206, pi_0_1=0.127451, Seed=1, 
Monotonicity=c("General"), Sum_Pi_f = seq(from=0.01, to=.99, by=.01), M=1000)

# Obtain a density plot of R2_H, assuming that 
# r_h_S0S1>=.2, r_h_S0T1>=0, r_h_T0T1>=.2, and r_h_T0S1>=0
ICA.BinBin.CounterAssum(ICA, r_h_S0S1_min=.2, r_h_S0S1_max=1, 
r_h_S0T1_min=0, r_h_S0T1_max=1, r_h_T0T1_min=0.2, r_h_T0T1_max=1, 
r_h_T0S1_min=0, r_h_T0S1_max=1, Monotonicity="General",
Type="Density") 

# Now show the densities of R2_H under the different 
# monotonicity assumptions 
ICA.BinBin.CounterAssum(ICA, r_h_S0S1_min=.2, r_h_S0S1_max=1, 
r_h_S0T1_min=0, r_h_S0T1_max=1, r_h_T0T1_min=0.2, r_h_T0T1_max=1, 
r_h_T0S1_min=0, r_h_T0S1_max=1, Monotonicity="General",
Type="All.Densities", MainPlot=" ", Cex.Legend=1, 
Cex.Position="topright", ylim=c(0, 20)) 
}
}
\keyword{Plot surrogacy}
\keyword{Causal-Inference framework}
\keyword{Single-trial setting}