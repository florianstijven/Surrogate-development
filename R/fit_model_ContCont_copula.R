#' Fit continuous-continuous vine copula model
#'
#' [fit_copula_ContCont()] fits the continuous-continuous vine copula model. See
#' Details for more information about this model.
#'
#' @param marginal_S0,marginal_S1,marginal_T0,marginal_T1 List with the
#'   following three elements (in order):
#' * Density function with first argument `x` and second argument `para` the parameter
#'   vector for this distribution.
#' * Distribution function with first argument `x` and second argument `para` the parameter
#'   vector for this distribution.
#' * Inverse distribution function with first argument `p` and second argument `para` the parameter
#'   vector for this distribution.
#' * The number of elements in `para`.
#' * A vector of starting values for `para`.
#' @inheritParams fit_copula_OrdCont
#' @inherit continuous_continuous_loglik details
#'
#' @return Returns an S3 object that can be used to perform the sensitivity
#'   analysis with [sensitivity_analysis_copula()].
#' @export
#'
#' @author Florian Stijven
#'
#' @seealso [sensitivity_analysis_copula()], [print.vine_copula_fitted()],
#'   [plot.vine_copula_fitted()]
fit_copula_ContCont = function(data,
                              copula_family,
                              marginal_S0,
                              marginal_S1,
                              marginal_T0,
                              marginal_T1,
                              start_copula,
                              method = "BFGS",
                              maxit = 500,
                              copula_transform = function(x) x,
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

  # Starting values may be specified in the marginal_S0 and marginal_S1 lists.
  if (length(marginal_S0) >= 5) {
    start_S0 = marginal_S0[[5]]
  }
  else start_S0 = rep(1, marginal_S0[[4]])
  if (length(marginal_S1) >= 5) {
    start_S1 = marginal_S1[[5]]
  }
  else start_S1 = rep(1, marginal_S1[[4]])

  if (length(marginal_T0) >= 5) {
    start_T0 = marginal_T0[[5]]
  }
  else start_T0 = rep(1, marginal_T0[[4]])
  if (length(marginal_T1) >= 5) {
    start_T1 = marginal_T1[[5]]
  }
  else start_T1 = rep(1, marginal_T1[[4]])

  submodel_0 = fit_copula_submodel_ContCont(
    X = data0$surr,
    Y = data0$true,
    copula_family = copula_family[1],
    marginal_X = marginal_S0,
    marginal_Y = marginal_T0,
    start_X = start_S0,
    start_Y = start_T0,
    start_copula = start_copula,
    method = method,
    copula_transform = copula_transform,,
    maxit = maxit,
    ...
  )
  submodel_1 = fit_copula_submodel_ContCont(
    X = data1$surr,
    Y = data1$true,
    copula_family = copula_family[2],
    marginal_X = marginal_S1,
    marginal_Y = marginal_T1,
    start_X = start_S1,
    start_Y = start_T1,
    start_copula = start_copula,
    method = method,
    copula_transform = copula_transform,
    maxit = maxit,
    ...
  )


  return(
    new_vine_copula_fit(submodel_0, submodel_1, c("continuous", "continuous"))
  )
}

#' Fit ordinal-continuous copula submodel
#'
#' The `fit_copula_submodel_ContCont()` function fits the copula (sub)model for
#' a continuous surrogate and true endpoint with maximum likelihood.
#'
#' @param start_X,start_Y Starting values corresponding to `marginal_X` and
#'   `marginal_Y`.
#' @param copula_transform Used for reparameterizing the copula parameter.
#'   [copula_transform()] backtransforms the transformed copula parameter to the
#'   original scale. Note that `start_copula` should be specified on the
#'   transformed scale.
#' @inheritParams continuous_continuous_loglik
#' @inheritParams fit_copula_submodel_OrdCont
#'
#' @inherit fit_copula_submodel_OrdCont return
#' @seealso [continuous_continuous_loglik()]
fit_copula_submodel_ContCont = function(X,
                                        Y,
                                        copula_family,
                                        marginal_X,
                                        marginal_Y,
                                        start_X,
                                        start_Y,
                                        start_copula,
                                        method = "BFGS",
                                        names_XY = c("Surr", "True"),
                                        twostep = FALSE,
                                        maxit,
                                        copula_transform = function(x) x,
                                        ...)
{
  # Number of parameters for X.
  p1 = marginal_X[[4]]
  # Number of parameters for Y.
  p2 = marginal_Y[[4]]

  # loglikelihood function to be maximized.
  log_lik_function = function(para) {
    para[p1 + p2 + 1] = copula_transform(para[p1 + p2 + 1])
    continuous_continuous_loglik(
      para = para,
      X = X,
      Y = Y,
      copula_family = copula_family,
      marginal_X = marginal_X,
      marginal_Y = marginal_Y,
      return_sum = FALSE
    )
  }
  # Estimate marginal distribution of X.
  param_X = estimate_marginal(X, marginal_X, start_X)
  names(param_X) = paste0(
    names_XY[1],
    "[",
    1:p1,
    "]"
  )
  # Estimate marginal distribution of Y.
  param_Y = estimate_marginal(Y, marginal_Y, start_Y)
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
    two_step_fit = fit_copula_submodel_ContCont(X,
                                                Y,
                                                copula_family,
                                                marginal_X,
                                                marginal_Y,
                                                start_X,
                                                start_Y,
                                                start_copula,
                                                twostep = TRUE,
                                                copula_transform = copula_transform)

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
    marginal_X = marginal_cont_constructor(marginal_X, est_X),
    marginal_Y = marginal_cont_constructor(marginal_Y, est_Y),
    copula_family = copula_family,
    data = data.frame(X = X, Y = Y),
    names_XY = names_XY
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
#' @param marginal_X,marginal_Y List with the following three elements (in order):
#' * Density function with first argument `x` and second argument `para` the parameter
#' vector for this distribution.
#' * Distribution function with first argument `x` and second argument `para` the parameter
#' vector for this distribution.
#' * Inverse distribution function with first argument `p` and second argument `para` the parameter
#' vector for this distribution.
#' * The number of elements in `para`.
#' * A vector of starting values for `para`.
#' @inheritParams loglik_copula_scale
#'
#' @return (numeric) loglikelihood value evaluated in `para`.
continuous_continuous_loglik <- function(para, X, Y, copula_family, marginal_X, marginal_Y, return_sum = TRUE){
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
      pdf_Y = pdf_Y,
      return_sum = return_sum
    )

  return(loglik)
}

conditional_mean_copula_ContCont = function(fitted_submodel, grid) {
  para = coef(fitted_submodel$ml_fit)
  pdf_x = fitted_submodel$marginal_X$pdf
  marginal_X = attr(fitted_submodel$marginal_X, "constructor")
  marginal_Y = attr(fitted_submodel$marginal_Y, "constructor")
  # Compute the expected value through numerical integration. First, define a
  # helper function that computes the bivariate estimated density.
  log_dens_joint = function(x, y) {
    log_dens = continuous_continuous_loglik(
      para = para,
      X = x,
      Y = y,
      copula_family = fitted_submodel$copula_family,
      marginal_X = marginal_X,
      marginal_Y = marginal_Y,
      return_sum = FALSE
    )
    log_dens[is.nan(log_dens)] = -Inf
    return(log_dens)
  }

  conditional_mean = sapply(grid,
         function(x) {
           constant = log(pdf_x(x))
           integrand = function(y) {
             y * exp(log_dens_joint(rep(x, length(y)), y) - constant)
           }

           cond_mean = stats::integrate(
             f = integrand,
             lower = -Inf,
             upper = +Inf
           )$value
           return(cond_mean)
         })
  return(conditional_mean)
}
