% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensitivity_analysis_copula.R
\name{compute_ICA}
\alias{compute_ICA}
\title{Compute Individual Causal Association for a given D-vine copula model in the
setting of choice.}
\usage{
compute_ICA(endpoint_types, ...)
}
\arguments{
\item{endpoint_types}{(character) vector with two elements indicating the
endpoint types: \code{"continuous"} or \code{"ordinal"}.}

\item{...}{Arguments to pass onto \code{\link[=compute_ICA_ContCont]{compute_ICA_ContCont()}},
\code{\link[=compute_ICA_OrdCont]{compute_ICA_OrdCont()}}, or \code{\link[=compute_ICA_OrdOrd]{compute_ICA_OrdOrd()}}}
}
\value{
(numeric) A Named vector with the following elements:
\itemize{
\item ICA
\item Spearman's rho, \eqn{\rho_s (\Delta S, \Delta T)} (if asked)
\item Marginal association parameters in terms of Spearman's rho (if asked):
\deqn{\rho_{s}(T_0, S_0), \rho_{s}(T_0, S_1), \rho_{s}(T_0, T_1),
 \rho_{s}(S_0, S_1), \rho_{s}(S_0, T_1),
 \rho_{s}(S_1, T_1)}
}
}
\description{
The \code{\link[=compute_ICA]{compute_ICA()}} function computes the individual causal
association for a fully identified D-vine copula model. See details for the
default definition of the ICA in each setting.
}
