\name{RandVec}
\alias{RandVec}
\title{
Generate random vectors with a fixed sum
}
\description{
This function generates an n by m array x, each of whose m columns contains n random values lying in the interval [a,b], subject to the condition that their sum be equal to s. The distribution of values is uniform in the sense that it has the conditional probability distribution of a uniform distribution over the whole n-cube, given that the sum of the x's is s.
The function uses the \code{randfixedsum} algorithm, written by Roger Stafford and implemented in MatLab. For details, see http://www.mathworks.com/matlabcentral/fileexchange/9700-random-vectors-with-fixed-sum/content/randfixedsum.m 
}
\usage{
RandVec(a=0, b=1, s=1, n=9, m=1, Seed=sample(1:1000, size = 1))
}
\arguments{
  \item{a}{The function \code{RandVec} generates an \code{n} by \code{m} matrix x. Each of the \code{m} columns contain \code{n} random values lying in the interval [a,b]. The argument \code{a} specifies the lower limit of the interval. Default \code{0}.}
  \item{b}{The argument \code{b} specifies the upper limit of the interval. Default \code{1}.}
  \item{s}{The argument \code{s} specifies the value to which each of the \code{m} generated columns should sum to. Default \code{1}.}
  \item{n}{The number of requested elements per column. Default \code{9}.}
  \item{m}{The number of requested columns. Default \code{1}.}
  \item{Seed}{The seed that is used. Default \code{sample(1:1000, size = 1)}.}
}
\value{An object of class \code{RandVec} with components,
\item{RandVecOutput}{The randomly generated vectors.}
}
\references{
The function is an R adaptation of a matlab program written by Roger Stafford. For details on the original Matlab algorithm, see: 
http://www.mathworks.com/matlabcentral/fileexchange/9700-random-vectors-with-fixed-sum/content/randfixedsum.m
}
\author{
Wim Van der Elst, Ariel Alonso, & Geert Molenberghs
}
\examples{
# generate two vectors with 10 values ranging between 0 and 1
# where each vector sums to 1
# (uniform distribution over the whole n-cube)
Vectors <- RandVec(a=0, b=1, s=1, n=10, m=2)
sum(Vectors$RandVecOutput[,1])
sum(Vectors$RandVecOutput[,2])
}
\keyword{RandVec}
