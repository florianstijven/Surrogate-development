\name{Sim.Data.MTS}
\alias{Sim.Data.MTS}
\title{Simulates a dataset that can be used to assess surrogacy in the multiple-trial setting}
\description{
The function \code{Sim.Data.MTS} simulates a dataset that contains the variables Treat, Trial.ID, Surr, True, and Pat.ID. The user can specify the number of patients and the number of trials that should be included in the simulated dataset, the desired \eqn{R_{trial}} and \eqn{R_{indiv}} values, the desired variability of the trial-specific treatment effects for the surrogate and the true endpoints (i.e., \eqn{d_{aa}} and \eqn{d_{bb}}, respectively), and the desired fixed-effect parameters of the intercepts and treatment effects for the surrogate and the true endpoints. 
}  
\usage{Sim.Data.MTS(N.Total=2000, N.Trial=50, R.Trial.Target=.8, R.Indiv.Target=.8, 
Fixed.Effects=c(0, 0, 0, 0), D.aa=10, D.bb=10, Seed=sample(1:1000, size=1), 
Model=c("Full"))
}                                                             
\arguments{
  \item{N.Total}{The total number of patients in the simulated dataset. Default \eqn{2000}.}
  \item{N.Trial}{The number of trials. Default \eqn{50}.}
  \item{R.Trial.Target}{The desired \eqn{R_{trial}} value in the sumilated dataset. Default \eqn{0.80}}
  \item{R.Indiv.Target}{The desired \eqn{R_{indiv}} value in the simulated dataset. Default \eqn{0.80}.}
  \item{Fixed.Effects}{A vector that specifies the desired fixed-effect intercept for the surrogate, fixed-effect intercept for the true endpoint, fixed treatment effect for the surrogate, and fixed treatment effect for the true endpoint, respectively. Default \code{c(0, 0, 0, 0)}.}
  \item{D.aa}{The desired variability of the trial-specific treatment effects on the surrogate endpoint. Default \eqn{10}.}
  \item{D.bb}{The desired variability of the trial-specific treatment effects on the true endpoint. Default \eqn{10}.}
  \item{Model}{The type of model that will be fitted on the data when surrogacy is assessed, i.e., a full, semireduced, or reduced model (for details, see \code{\link{UnifixedContCont}}, \code{\link{UnimixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{BimixedContCont}}).}
\item{Seed}{The seed that is used to generate the dataset. Default \code{sample(x=1:1000, size=1)}, i.e., a random number between 1 and 1000.}
}
\details{
The generated object \code{Data.Observed.MTS} (of class \code{data.frame}) is placed in the workspace (for easy access).

The number of patients per trial in the simulated dataset is identical in each trial, and equals the requested total number of patients divided by the requested number of trials (=\code{N.Total/N.Trial}). If this is not a whole number, a warning is given and the number of patients per trial is automatically rounded up to the nearest whole number. See \bold{Examples} below.

Treatment allocation is balanced when the number of patients per trial is an odd number. If this is not the case, treatment allocation is balanced up to one patient (the remaining patient is randomly allocated to the exprimental or the control treatment groups in each of the trials).     
}

\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}

\seealso{\code{\link{UnifixedContCont}}, \code{\link{BifixedContCont}}, \code{\link{UnimixedContCont}}, \code{\link{BimixedContCont}}, \code{\link{Sim.Data.STS}}
}

\examples{
# Simulate a dataset with 2000 patients, 50 trials, Rindiv=Rtrial=.8, D.aa=10,
# D.bb=50, and fixed effect values 1, 2, 30, and 90:
Sim.Data.MTS(N.Total=2000, N.Trial=50, R.Trial.Target=.8, R.Indiv.Target=.8, D.aa=10, 
D.bb=50, Fixed.Effects=c(1, 2, 30, 90), Seed=1)  

# Sample output, the first 10 rows of Data.Observed.MTS:
Data.Observed.MTS[1:10,]

# Note: When the following code is used to generate a dataset:
Sim.Data.MTS(N.Total=2000, N.Trial=99, R.Trial.Target=.5, R.Indiv.Target=.8, 
D.aa=10, D.bb=50, Fixed.Effects=c(1, 2, 30, 90), Seed=1)  

# R gives the following warning: 

# > NOTE: The number of patients per trial requested in the function call 
# > equals 20.20202 (=N.Total/N.Trial), which is not a whole number.  
# > To obtain a dataset where the number of patients per trial is balanced for 
# > all trials, the number of patients per trial was rounded to 21 to generate 
# > the dataset. Data.Observed.MTS thus contains a total of 2079 patients rather 
# > than the requested 2000 in the function call.
}
\keyword{Simulate data}
\keyword{Multiple-trial setting}
