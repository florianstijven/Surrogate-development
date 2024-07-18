#' Loglikelihood function for binary-continuous copula model
#'
#' [ordinal_continuous_loglik()] computes the observed-data loglikelihood for a
#' bivariate copula model with a continuous and an ordinal endpoint. The model
#' is based on a latent variable representation of the ordinal endpoint.
#'
#' @details
#' # Vine Copula Model for Ordinal Endpoints
#'
#' Following the Neyman-Rubin potential outcomes framework, we assume that each
#' patient has four potential outcomes, two for each arm, represented by
#' \eqn{\boldsymbol{Y} = (T_0, S_0, S_1, T_1)'}. Here, \eqn{\boldsymbol{Y_z} =
#' (S_z, T_z)'} are the potential surrogate and true endpoints under treatment
#' \eqn{Z = z}. We will further assume that \eqn{T} is ordinal and \eqn{S} is
#' continuous; consequently, the function argument `X` corresponds to \eqn{T} and
#' `Y` to \eqn{S}. (The roles of \eqn{S} and \eqn{T} can be interchanged without
#' loss of generality.)
#'
#'
#' We introduce latent variables to model \eqn{\boldsymbol{Y}}. Latent variables
#' will be denoted by a tilde. For instance, if \eqn{T_z} is ordinal with \eqn{K_T}
#' categories, then \eqn{T_z} is a function of the latent
#' \eqn{\tilde{T}_z \sim N(0, 1)} as follows:
#' \deqn{
#' T_z = g_{T_z}(\tilde{T}_z; \boldsymbol{c}^{T_z}) = \begin{cases}
#' 1 & \text{ if } -\infty = c_0^{T_z} < \tilde{T_z} \le c_1^{T_z} \\
#' \vdots \\
#' k & \text{ if } c_{k - 1}^{T_z} < \tilde{T_z} \le c_k^{T_z} \\
#' \vdots \\
#' K & \text{ if } c_{K_{T} - 1}^{T_z} < \tilde{T_z} \le c_{K_{T}}^{T_z} = \infty, \\
#' \end{cases}
#' }
#' where \eqn{\boldsymbol{c}^{T_z} = (c_1^{T_z}, \cdots, c_{K_T - 1}^{T_z})}.
#' The latent counterpart of \eqn{\boldsymbol{Y}} is again denoted by a tilde;
#' for example, \eqn{\tilde{\boldsymbol{Y}} = (\tilde{T}_0, S_0, S_1, \tilde{T}_1)'}
#' if \eqn{T_z} is ordinal and \eqn{S_z} is continuous.
#'
#' The vector of latent potential outcome \eqn{\tilde{\boldsymbol{Y}}} is modeled
#' with a D-vine copula as follows:
#' \deqn{
#' f_{\tilde{\boldsymbol{Y}}} =  f_{\tilde{T}_0} \, f_{S_0} \, f_{S_1} \, f_{\tilde{T}_1}
#' \cdot c_{\tilde{T}_0, S_0 } \, c_{S_0, S_1} \, c_{S_1, \tilde{T}_1}
#' \cdot c_{\tilde{T}_0, S_1; S_0} \, c_{S_0, \tilde{T}_1; S_1}
#' \cdot c_{\tilde{T}_0, \tilde{T}_1; S_0, S_1},
#' }
#' where (i) \eqn{f_{T_0}}, \eqn{f_{S_0}}, \eqn{f_{S_1}}, and \eqn{f_{T_1}} are
#' univariate density functions, (ii) \eqn{c_{T_0, S_0}}, \eqn{c_{S_0, S_1}},
#' and \eqn{c_{S_1, T_1}} are unconditional bivariate copula densities, and (iii)
#' \eqn{c_{T_0, S_1; S_0}}, \eqn{c_{S_0, T_1; S_1}}, and \eqn{c_{T_0, T_1; S_0, S_1}}
#' are conditional bivariate copula densities (e.g., \eqn{c_{T_0, S_1; S_0}}
#' is the copula density of \eqn{(T_0, S_1)' \mid S_0}. We also make the
#' simplifying assumption for all copulas.
#'
#' # Observed-Data Likelihood
#'
#' In practice, we only observe \eqn{(S_0, T_0)'} or \eqn{(S_1, T_1)'}. Hence, to
#' estimate the (identifiable) parameters of the D-vine copula model, we need
#' to derive the observed-data likelihood. The observed-data loglikelihood for
#' \eqn{(S_z, T_z)'} is as follows:
#' \deqn{
#' f_{\boldsymbol{Y_z}}(s, t; \boldsymbol{\beta}) =
#' \int_{c^{T_z}_{t - 1}}^{+ \infty} f_{\boldsymbol{\tilde{Y}_z}}(s, x; \boldsymbol{\beta}) \, dx - \int_{c^{T_z}_{t}}^{+ \infty} f_{\boldsymbol{\tilde{Y}_z}}(s, x; \boldsymbol{\beta}) \, dx.
#' }
#' The above expression is used in [ordinal_continuous_loglik()] to compute the
#' loglikelihood for the observed values for \eqn{Z = 0} or \eqn{Z = 1}. In this
#' function, `X` and `Y` correspond to \eqn{T_z} and \eqn{S_z} if \eqn{T_z} is
#' ordinal and \eqn{S_z} continuous. Otherwise, `X` and `Y` correspond to
#' \eqn{S_z} and \eqn{T_z}.
#'
#'
#' @param para Parameter vector. The parameters are ordered as follows:
#' * `para[1:p1]`: Cutpoints for the latent distribution of X corresponding to
#' \eqn{c_1, \dots, c_{K - 1}} (see Details).
#' * `para[(p1 + 1):(p1 + p2)]`: Parameters for surrogate distribution, more details in
#'  `?Surrogate::cdf_fun` for the specific implementations.
#' * `para[p1 + p2 + 1]`: copula parameter
#' @param X First variable (Ordinal with \eqn{K} categories)
#' @param Y Second variable (Continuous)
#' @param marginal_Y List with the following three elements (in order):
#' * Density function with first argument `x` and second argument `para` the parameter
#' vector for this distribution.
#' * Distribution function with first argument `x` and second argument `para` the parameter
#' vector for this distribution.
#' * The number of elements in `para`.
#' @param K Number of categories in `X`.
#' @inheritParams loglik_copula_scale
#'
#' @return (numeric) loglikelihood value evaluated in `para`.
ordinal_continuous_loglik <- function(para, X, Y, copula_family, marginal_Y, K){
  # Number of independent cut-points for the distribution of X.
  p1 = K - 1
  # Parameters for the distribution of X.
  para_X = para[1:p1]
  # Number of parameters for the distribution of Y.
  p2 = marginal_Y[[3]]
  # Parameters for the distribution of Y.
  para_Y = para[(p1 + 1):(p1 + p2)]
  # Vector of copula parameters.
  theta = para[p1 + p2 + 1]

  # Construct marginal distribution and density functions.
  cdf_X = pnorm
  pdf_X = dnorm

  pdf_Y = function(x) {
    marginal_Y[[1]](x, para_Y)
  }
  cdf_Y = function(x) {
    marginal_Y[[2]](x, para_Y)
  }

  # We compute the individual likelihood contributions in terms of the difference
  # of two right-censored contributions.
  lik1 = exp(
    log_likelihood_copula_model(
      theta = theta,
      X = ordinal_to_cutpoints(X, para_X, TRUE),
      Y = Y,
      d1 = rep(0, length(X)),
      d2 = rep(1, length(X)),
      copula_family = copula_family,
      cdf_X = cdf_X,
      cdf_Y = cdf_Y,
      pdf_X = pdf_X,
      pdf_Y = pdf_Y,
      return_sum = FALSE
    )
  )
  lik2 = exp(
    log_likelihood_copula_model(
      theta = theta,
      X = ordinal_to_cutpoints(X, para_X, FALSE),
      Y = Y,
      d1 = rep(0, length(X)),
      d2 = rep(1, length(X)),
      copula_family = copula_family,
      cdf_X = cdf_X,
      cdf_Y = cdf_Y,
      pdf_X = pdf_X,
      pdf_Y = pdf_Y,
      return_sum = FALSE
    )
  )
  # The individual likelihood contributions (corresponding to interval-censored
  # data) is the difference of the two right-censored contributions.
  lik = lik1 - lik2
  # We take the log and sum to compute the observed-data log-likelihood.
  loglik = sum(log(lik))

  return(loglik)
}

#' Convert Ordinal Observations to Latent Cutpoints
#'
#' [ordinal_to_cutpoints()] converts the ordinal endpoints to the corresponding
#' cutpoints of the underlying latent continuous variable. Let
#' \eqn{P(x \le k) = G(c_k)} where \eqn{G} is the distribution function of the
#' latent variable. [ordinal_to_cutpoints()] converts \eqn{x} to \eqn{c_k} (or to
#' \eqn{c_{k - 1}} if `strict = TRUE`.
#'
#' @param x Integer vector with values in `1:(length(cutpoints) + 1)`.
#' @param cutpoints The cutpoints on the latent scale corresponding to
#' \eqn{\boldsymbol{c} = c(c_1, \cdots, c_{K - 1})}.
#'
#' @return Numeric vector with cutpoints corresponding to the values in `x`.
ordinal_to_cutpoints = function(x, cutpoints, strict) {
  if (strict) {
    # Add c_0 to the cutpoints
    cutpoints = c(-Inf, cutpoints)
  }
  else {
    # Add c_K to the cutpoints
    cutpoints = c(cutpoints, +Inf)
  }
  return(cutpoints[x])
}


