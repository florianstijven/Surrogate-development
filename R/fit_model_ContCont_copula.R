#' Fit ordinal-continuous copula submodel
#'
#' The `fit_copula_submodel_ContCont()` function fits the copula (sub)model for a
#' continuous surrogate and true endpoint with maximum likelihood.
#'
#' @param name description
#' @inheritParams continuous_continuous_loglik
#' @inherit fit_copula_submodel_OrdCont return
fit_copula_submodel_ContCont = function(X,
                                        Y,
                                        copula_family,
                                        marginal_X,
                                        marginal_Y,
                                        start_X,
                                        start_Y,
                                        start_copula,
                                        method = "BFGS",
                                        twostep = FALSE)
{
  # Number of parameters for X.
  p1 = marginal_X[[4]]
  # Number of parameters for Y.
  p2 = marginal_Y[[4]]

  # loglikelihood function to be maximized.
  log_lik_function = function(para) {
    continuous_continuous_loglik(
      para = para,
      X = X,
      Y = Y,
      copula_family = copula_family,
      marginal_X = marginal_X,
      marginal_Y = marginal_Y
    )
  }
  # Estimate marginal distribution of X.
  param_X = estimate_marginal(X, marginal_X, start_X)
  names(param_X) = paste0(
    "X[",
    1:p1,
    "]"
  )
  # Estimate marginal distribution of Y.
  param_Y = estimate_marginal(Y, marginal_Y, start_Y)
  names(param_Y) = paste0(
    "Y[",
    1:p2,
    "]"
  )

  # If twostep is TRUE, we compute the estimate in two stages. The marginal
  # distributions have already been estimated (param_X and param_Y). The
  # estimates are fixed in the twostage estimator.

  if (twostep) {
    fixed = 1:(p1 + p2)
  }
  else {
    fixed = NULL
    # Starting values for the full estimator are determined by the twostage
    # estimates.
    two_step_fit = fit_copula_submodel_ContCont(
      X,
      Y,
      copula_family,
      marginal_X,
      marginal_Y,
      start_X,
      start_Y,
      start_copula,
      twostep = TRUE
    )
    start_copula = coef(two_step_fit$ml_fit)[p1 + p2 + 1]
  }
  start = c(param_X, param_Y, "theta (copula)" = start_copula)

  suppressWarnings({
    ml_fit = maxLik::maxLik(
      logLik = log_lik_function,
      start = start,
      method = method,
      fixed = fixed
    )
  })

  est_X = coef(ml_fit)[1:p1]
  est_Y = coef(ml_fit)[(p1 + 1):(p1 + p2)]
  # Return a list with the marginal surrogate distribution, fit object for
  # copula parameter, and element to indicate the copula family.
  submodel_fit = list(
    ml_fit = ml_fit,
    marginal_X = marginal_cont_constructor(marginal_X, est_X),
    marginal_Y = marginal_cont_constructor(marginal_Y, est_Y),
    copula_family = copula_family
  )
  return(submodel_fit)
}

#' Loglikelihood function for continuous-continuous copula model
#'
#' [continuous_continuous_loglik()] computes the observed-data loglikelihood for a
#' bivariate copula model with two continuous endpoints.
#'
#'
#' @param para Parameter vector. The parameters are ordered as follows:
#' * `para[1:p1]`: Parameters for the distribution of `X` as specified in
#' `marginal_X`.
#' * `para[(p1 + 1):(p1 + p2)]`: Parameters for the distribution of `Y` as specified in
#' `marginal_Y`.
#' * `para[p1 + p2 + 1]`: copula parameter
#' @param X First variable (Continuous)
#' @param Y Second variable (Continuous)
#' @param marginal_X List with the following three elements (in order):
#' * Density function with first argument `x` and second argument `para` the parameter
#' vector for this distribution.
#' * Distribution function with first argument `x` and second argument `para` the parameter
#' vector for this distribution.
#' * The number of elements in `para`.
#' @param marginal_Y See `marginal_X`
#' @inheritParams loglik_copula_scale
#'
#' @return (numeric) loglikelihood value evaluated in `para`.
continuous_continuous_loglik <- function(para, X, Y, copula_family, marginal_X, marginal_Y){
  # Number of parameters for the distribution of X.
  p1 = marginal_X[[4]]
  # Parameters for the distribution of X.
  para_X = para[1:p1]
  # Number of parameters for the distribution of Y.
  p2 = marginal_Y[[4]]
  # Parameters for the distribution of Y.
  para_Y = para[(p1 + 1):(p1 + p2)]
  # Vector of copula parameters.
  theta = para[p1 + p2 + 1]

  # Construct marginal distribution and density functions.
  pdf_X = function(x) {
    marginal_X[[1]](x, para_X)
  }
  cdf_X = function(x) {
    marginal_X[[2]](x, para_X)
  }

  pdf_Y = function(x) {
    marginal_Y[[1]](x, para_Y)
  }
  cdf_Y = function(x) {
    marginal_Y[[2]](x, para_Y)
  }

  # We compute the observed-data loglikelihood.
  loglik =
    log_likelihood_copula_model(
      theta = theta,
      X = X,
      Y = Y,
      d1 = rep(1, length(X)),
      d2 = rep(1, length(X)),
      copula_family = copula_family,
      cdf_X = cdf_X,
      cdf_Y = cdf_Y,
      pdf_X = pdf_X,
      pdf_Y = pdf_Y
    )

  return(loglik)
}
