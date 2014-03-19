\name{LongToWide}
\alias{LongToWide}
\title{
Reshapes a dataset from the 'long' format (i.e., multiple lines per patient) into the 'wide' format (i.e., one line per patient)
}
\description{
Reshapes a dataset that is in the 'long' format into the 'wide' format. The dataset should contain a single surrogate endpoint and a single true endpoint value per subject.  
}
\usage{
LongToWide(Dataset, OutcomeIndicator, IdIndicator, TreatIndicator, OutcomeValue)
}
\arguments{
  \item{Dataset}{A \code{data.frame} in the 'long' format that contains (at least) five columns, i.e., one that contains the subject ID, one that contains the trial ID,  one that contains the endpoint indicator, one that contains the treatment indicator, and one that contains the endpoint values.}
  \item{OutcomeIndicator}{The name of the variable in \code{Dataset} that contains the indicator that distinguishes between the surrogate and true endpoints.}
  \item{IdIndicator}{The name of the variable in \code{Dataset} that contains the subject ID.}
  \item{TreatIndicator}{The name of the variable in \code{Dataset} that contains the treatment indicator. For the subsequent surrogacy analyses, the treatment indicator should either be coded as \eqn{1} for the experimental group and \eqn{-1} for the control group, or as \eqn{1} for the experimental group and \eqn{0} for the control group. The \eqn{-1/1} coding is recommended. }
  \item{OutcomeValue}{The name of the variable in \code{Dataset} that contains the endpoint values.}
}
\value{
A \code{data.frame} in the 'wide' format, i.e., a \code{data.frame} that contains one line per subject. Each line contains a surrogate value, a true endpoint value, a treatment indicator, a patient ID, and a trial ID.}
\author{
Wim Van der Elst, Ariel Alonso, and Geert Molenberghs}
\examples{
# Generate a dataset in the 'long' format that contains 
# S and T values for 100 patients
Outcome <- rep(x=c(0, 1), times=100)
ID <- rep(seq(1:100), each=2)
Treat <- rep(seq(c(0,1)), each=100)
Outcomes <- as.numeric(matrix(rnorm(1*200, mean=100, sd=10), 
                                      ncol=200))
Data <- data.frame(cbind(Outcome, ID, Treat, Outcomes))

# Reshapes the Data object 
LongToWide(Dataset=Data, OutcomeIndicator=Outcome, IdIndicator=ID, 
           TreatIndicator=Treat, OutcomeValue=Outcomes)
}
\keyword{Transpose dataset}
