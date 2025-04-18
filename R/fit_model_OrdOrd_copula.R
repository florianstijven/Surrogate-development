#' Fit ordinal-ordinal vine copula model
#'
#' [fit_copula_OrdOrd()] fits the ordinal-ordinal vine copula model. See
#' Details for more information about this model.
#'
#' @param K_S,K_T Number of categories in the surrogate and true endpoints.
#' @inheritParams fit_copula_OrdCont
#' @inheritParams fit_copula_submodel_OrdOrd
#' @inherit ordinal_ordinal_loglik details
#'
#' @inherit fit_copula_OrdCont return seealso
#' @export
#'
#' @author Florian Stijven
fit_copula_OrdOrd = function(data,
                              copula_family,
                              K_S,
                              K_T,
                              start_copula,
                              method = "BFGS",
                              ...) {
  # If copula_family is length 1, we repeat the same copula family.
  if (length(copula_family) == 1) {
    copula_family = rep(copula_family, 2)
  }

  # Column names are added to make the intrepretation of the further code
  # easier. Pfs refers to the surrogate, Surv refers to the true endpoint.
  colnames(data) = c("surr", "true", "treat")

  # Split original dataset into two data sets, one for each treatment group.
  data0 = data[data$treat == 0, ]
  data1 = data[data$treat == 1, ]

  submodel_0 = fit_copula_submodel_OrdOrd(
    X = data0$surr,
    Y = data0$true,
    copula_family = copula_family[1],
    K_X = K_S,
    K_Y = K_T,
    start_copula = start_copula,
    method = method,
    ...
  )
  submodel_1 = fit_copula_submodel_OrdOrd(
    X = data1$surr,
    Y = data1$true,
    copula_family = copula_family[2],
    K_X = K_S,
    K_Y = K_T,
    start_copula = start_copula,
    method = method,
    ...
  )

  return(
    new_vine_copula_fit(submodel_0, submodel_1, c("ordinal", "ordinal"))
  )
}


#' Fit ordinal-continuous copula submodel
#'
#' The `fit_copula_submodel_OrdOrd()` function fits the copula (sub)model for an
#' ordinal surrogate and true endpoint with maximum likelihood.
#'
#' @inheritParams ordinal_ordinal_loglik
#' @inheritParams fit_copula_submodel_OrdCont
#' @inherit fit_copula_submodel_OrdCont return
fit_copula_submodel_OrdOrd = function(X,
                                      Y,
                                      copula_family,
                                      start_copula,
                                      method = "BFGS",
                                      K_X,
                                      K_Y,
                                      names_XY = c("Surr", "True"),
                                      twostep = FALSE,
                                      ...) {
  # Number of parameters for X.
  p1 = K_X - 1
  # Number of parameters for Y.
  p2 = K_Y - 1

  # loglikelihood function to be maximized.
  log_lik_function = function(para) {
    ordinal_ordinal_loglik(
      para = para,
      X = X,
      Y = Y,
      copula_family = copula_family,
      K_X = K_X,
      K_Y = K_Y,
      return_sum = FALSE
    )
  }
  # Compute empirical proportions for X.
  props_X = sapply(1:K_X, function(category) {
    mean(X == category)
  })
  # Compute cumulative probabilities. These are transformed to the normal cdf
  # scale.
  param_X = qnorm(cumsum(props_X[-length(props_X)]))
  names(param_X) = paste0(
    names_XY[1],
    "[",
    1:p1,
    "]"
  )

  # Do the same for Y.
  props_Y = sapply(1:K_Y, function(category) {
    mean(Y == category)
  })
  param_Y = qnorm(cumsum(props_Y[-length(props_Y)]))
  names(param_Y) = paste0(
    names_XY[2],
    "[",
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
    two_step_fit = fit_copula_submodel_OrdOrd(
      X,
      Y,
      copula_family,
      start_copula,
      K_X = K_X,
      K_Y = K_Y,
      twostep = TRUE
    )
    start_copula = coef(two_step_fit$ml_fit)[p1 + p2 + 1]
  }
  start = c(param_X, param_Y, start_copula)
  names(start)[length(start)] = "theta (copula)"
  suppressWarnings({
    ml_fit = maxLik::maxLik(
      logLik = log_lik_function,
      start = start,
      method = method,
      fixed = fixed,
      ...
    )
  })

  est_X = coef(ml_fit)[1:p1]
  est_Y = coef(ml_fit)[(p1 + 1):(p1 + p2)]
  # Return a list with the marginal surrogate distribution, fit object for
  # copula parameter, and element to indicate the copula family.
  submodel_fit = list(
    ml_fit = ml_fit,
    marginal_X = marginal_ord_constructor(est_X),
    marginal_Y = marginal_ord_constructor(est_Y),
    copula_family = copula_family,
    data = data.frame(X = X, Y = Y),
    names_XY = names_XY
  )
  return(submodel_fit)
}

#' Loglikelihood function for ordinal-ordinal copula model
#'
#' [ordinal_ordinal_loglik()] computes the observed-data loglikelihood for a
#' bivariate copula model with two ordinal endpoints. The model
#' is based on a latent variable representation of the ordinal endpoints.
#'
#' @details
#' ## Vine Copula Model for Ordinal Endpoints
#'
#' Following the Neyman-Rubin potential outcomes framework, we assume that each
#' patient has four potential outcomes, two for each arm, represented by
#' \eqn{\boldsymbol{Y} = (T_0, S_0, S_1, T_1)'}. Here, \eqn{\boldsymbol{Y_z} =
#' (S_z, T_z)'} are the potential surrogate and true endpoints under treatment
#' \eqn{Z = z}.
#'
#' The latent variable notation and D-vine copula model for \eqn{\boldsymbol{Y}}
#' is a straightforward extension of the notation in
#' [ordinal_continuous_loglik()].
#'
#' ## Observed-Data Likelihood
#'
#' In practice, we only observe \eqn{(S_0, T_0)'} or \eqn{(S_1, T_1)'}. Hence, to
#' estimate the (identifiable) parameters of the D-vine copula model, we need
#' to derive the observed-data likelihood. The observed-data loglikelihood for
#' \eqn{(S_z, T_z)'} is as follows:
#' \deqn{
#' f_{\boldsymbol{Y_z}}(s, t; \boldsymbol{\beta}) =
#' P \left( c^{S_z}_{s - 1} < \tilde{S}_z, c^{T_z}_{t - 1} < \tilde{T}_z  \right) - P \left( c^{S_z}_{s} < \tilde{S}_z, c^{T_z}_{t - 1} < \tilde{T}_z  \right)
#' - P \left( c^{S_z}_{s - 1} < \tilde{S}_z, c^{T_z}_{t} < \tilde{T}_z  \right) + P \left( c^{S_z}_{s} < \tilde{S}_z, c^{T_z}_{t} < \tilde{T}_z  \right).
#' }
#' The above expression is used in [ordinal_ordinal_loglik()] to compute the
#' loglikelihood for the observed values for \eqn{Z = 0} or \eqn{Z = 1}.
#'
#' @param para Parameter vector. The parameters are ordered as follows:
#' * `para[1:p1]`: Cutpoints for the latent distribution of X corresponding to
#' \eqn{c_1^X, \dots, c_{K_X - 1}^X} (see Details).
#' * `para[(p1 + 1):(p1 + p2)]`: Cutpoints for the latent distribution of Y corresponding to
#' \eqn{c_1^Y, \dots, c_{K_Y - 1}^Y} (see Details).
#' * `para[p1 + p2 + 1]`: copula parameter
#' @param X First variable (Ordinal with \eqn{K_X} categories)
#' @param Y Second variable (Ordinal with \eqn{K_Y} categories)
#' @param K_X Number of categories in `X`.
#' @param K_Y Number of categories in `Y`.
#' @inheritParams loglik_copula_scale
#'
#' @return (numeric) loglikelihood value evaluated in `para`.
ordinal_ordinal_loglik <- function(para, X, Y, copula_family, K_X, K_Y, return_sum = TRUE){
  # Number of independent cut-points for the distribution of X.
  p1 = K_X - 1
  # Parameters for the distribution of X.
  para_X = para[1:p1]
  # Number of parameters for the distribution of Y.
  p2 = K_Y - 1
  # Parameters for the distribution of Y.
  para_Y = para[(p1 + 1):(p1 + p2)]
  # Vector of copula parameters.
  theta = para[p1 + p2 + 1]

  # Construct marginal distribution and density functions.
  cdf_X = pnorm
  pdf_X = dnorm

  pdf_Y = dnorm
  cdf_Y = pnorm
  # We compute the individual likelihood contributions in terms of the linear
  # combination of 4 likelihoods.
  lik1 = exp(
    log_likelihood_copula_model(
      theta = theta,
      X = ordinal_to_cutpoints(X, para_X, TRUE),
      Y = ordinal_to_cutpoints(Y, para_Y, TRUE),
      d1 = rep(0, length(X)),
      d2 = rep(0, length(Y)),
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
      Y = ordinal_to_cutpoints(Y, para_Y, TRUE),
      d1 = rep(0, length(X)),
      d2 = rep(0, length(Y)),
      copula_family = copula_family,
      cdf_X = cdf_X,
      cdf_Y = cdf_Y,
      pdf_X = pdf_X,
      pdf_Y = pdf_Y,
      return_sum = FALSE
    )
  )
  lik3 = exp(
    log_likelihood_copula_model(
      theta = theta,
      X = ordinal_to_cutpoints(X, para_X, TRUE),
      Y = ordinal_to_cutpoints(Y, para_Y, FALSE),
      d1 = rep(0, length(X)),
      d2 = rep(0, length(Y)),
      copula_family = copula_family,
      cdf_X = cdf_X,
      cdf_Y = cdf_Y,
      pdf_X = pdf_X,
      pdf_Y = pdf_Y,
      return_sum = FALSE
    )
  )
  lik4 = exp(
    log_likelihood_copula_model(
      theta = theta,
      X = ordinal_to_cutpoints(X, para_X, FALSE),
      Y = ordinal_to_cutpoints(Y, para_Y, FALSE),
      d1 = rep(0, length(X)),
      d2 = rep(0, length(Y)),
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
  lik = lik1 - lik2 - lik3 + lik4
  loglik = log(lik)
  # We take the log and sum to compute the observed-data log-likelihood.
  if (return_sum) loglik = sum(loglik)

  return(loglik)
}
