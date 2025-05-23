% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_model_ContCont_copula.R
\name{continuous_continuous_loglik}
\alias{continuous_continuous_loglik}
\title{Loglikelihood function for continuous-continuous copula model}
\usage{
continuous_continuous_loglik(
  para,
  X,
  Y,
  copula_family,
  marginal_X,
  marginal_Y,
  return_sum = TRUE
)
}
\arguments{
\item{para}{Parameter vector. The parameters are ordered as follows:
\itemize{
\item \code{para[1:p1]}: Parameters for the distribution of \code{X} as specified in
\code{marginal_X}.
\item \code{para[(p1 + 1):(p1 + p2)]}: Parameters for the distribution of \code{Y} as specified in
\code{marginal_Y}.
\item \code{para[p1 + p2 + 1]}: copula parameter
}}

\item{X}{First variable (Continuous)}

\item{Y}{Second variable (Continuous)}

\item{copula_family}{Copula family, one of the following:
\itemize{
\item \code{"clayton"}
\item \code{"frank"}
\item \code{"gumbel"}
\item \code{"gaussian"}
}}

\item{marginal_X, marginal_Y}{List with the following three elements (in order):
\itemize{
\item Density function with first argument \code{x} and second argument \code{para} the parameter
vector for this distribution.
\item Distribution function with first argument \code{x} and second argument \code{para} the parameter
vector for this distribution.
\item Inverse distribution function with first argument \code{p} and second argument \code{para} the parameter
vector for this distribution.
\item The number of elements in \code{para}.
\item A vector of starting values for \code{para}.
}}

\item{return_sum}{Return the sum of the individual loglikelihoods? If \code{FALSE},
a vector with the individual loglikelihood contributions is returned.}
}
\value{
(numeric) loglikelihood value evaluated in \code{para}.
}
\description{
\code{\link[=continuous_continuous_loglik]{continuous_continuous_loglik()}} computes the observed-data loglikelihood for a
bivariate copula model with two continuous endpoints.
}
