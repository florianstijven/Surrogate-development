#' Fit ordinal-continuous vine copula model
#'
#' [fit_copula_OrdCont()] fits the ordinal-continuous vine copula model. See
#' Details for more information about this model.
#'
#' @param data data frame with three columns in the following order: surrogate
#'   endpoint, true endpoint, and treatment indicator (0/1 coding). Ordinal endpoints
#'   should be integers starting from `1`.
#' @param K_T Number of categories in the true endpoint.
#' @param marginal_S0,marginal_S1 List with the
#'   following three elements (in order):
#' * Density function with first argument `x` and second argument `para` the parameter
#'   vector for this distribution.
#' * Distribution function with first argument `x` and second argument `para` the parameter
#'   vector for this distribution.
#' * Inverse distribution function with first argument `p` and second argument `para` the parameter
#'   vector for this distribution.
#' * The number of elements in `para`.
#' * A vector of starting values for `para`.
#' @inheritParams fit_model_SurvSurv
#' @inheritParams fit_copula_ContCont
#' @inheritParams fit_copula_submodel_OrdCont
#' @inheritDotParams fit_copula_submodel_OrdCont
#' @inherit ordinal_continuous_loglik details
#'
#' @return Returns an S3 object that can be used to perform the sensitivity
#'   analysis with [sensitivity_analysis_copula()].
#' @export
#'
#' @author Florian Stijven
#'
#' @seealso [sensitivity_analysis_copula()], [print.vine_copula_fit()],
#'   [plot.vine_copula_fit()]
fit_copula_OrdCont = function(data,
                              copula_family,
                              marginal_S0,
                              marginal_S1,
                              K_T,
                              start_copula,
                              method = "BFGS",
                              ...
) {
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

  submodel_0 = fit_copula_submodel_OrdCont(
    X = data0$true,
    Y = data0$surr,
    copula_family = copula_family[1],
    marginal_Y = marginal_S0,
    K = K_T,
    start_Y = start_S0,
    start_copula = start_copula,
    method = method,
    names_XY = c("True", "Surr"),
    ...
  )
  submodel_1 = fit_copula_submodel_OrdCont(
    X = data1$true,
    Y = data1$surr,
    copula_family = copula_family[2],
    marginal_Y = marginal_S1,
    K = K_T,
    start_Y = start_S1,
    start_copula = start_copula,
    method = method,
    names_XY = c("True", "Surr"),
    ...
  )

  return(
    new_vine_copula_fit(submodel_0, submodel_1, c("ordinal", "continuous"))
  )
}

#' Fit ordinal-continuous copula submodel
#'
#' The `fit_copula_submodel_OrdCont()` function fits the copula (sub)model for a
#' continuous surrogate and an ordinal true endpoint with maximum likelihood.
#'
#' @param names_XY Names for `X` and `Y`, respectively.
#' @param twostep (boolean) If `TRUE`, the starting values are fixed for the
#'   marginal distributions and only the copula parameter is estimated.
#' @param start_Y Starting values for the marginal distribution paramters for `Y`.
#' @param start_copula Starting value for the copula parameter.
#' @param ... Extra argument to pass onto maxLik::maxLik
#' @inheritParams ordinal_continuous_loglik
#' @inheritParams fit_model_SurvSurv
#'
#' @return A list with five elements:
#' * ml_fit: object of class `maxLik::maxLik` that contains the estimated copula
#'  model.
#' * marginal_X: list with the estimated cdf, pdf/pmf, and inverse cdf for X.
#' * marginal_Y: list with the estimated cdf, pdf/pmf, and inverse cdf for X.
#' * copula_family: string that indicates the copula family
#' * data: data frame containing `X` and `Y`
#' * names_XY: The names (i.e., `"Surr"` and `"True"`) for `X` and `Y`
#' @seealso [ordinal_continuous_loglik()]
fit_copula_submodel_OrdCont = function(X,
                                       Y,
                                       copula_family,
                                       marginal_Y,
                                       start_Y,
                                       start_copula,
                                       method = "BFGS",
                                       K,
                                       names_XY = c("Surr", "True"),
                                       twostep = FALSE,
                                       ...) {
  # Number of parameters for X.
  p1 = K - 1
  # Number of parameters for Y.
  p2 = marginal_Y[[4]]

  # loglikelihood function to be maximized.
  log_lik_function = function(para) {
    ordinal_continuous_loglik(
      para = para,
      X = X,
      Y = Y,
      copula_family = copula_family,
      marginal_Y = marginal_Y,
      K = K,
      return_sum = FALSE
    )
  }
  # Compute empirical proportions for X.
  props = sapply(1:K, function(category) {
    mean(X == category)
  })
  # Compute cumulative probabilities. These are transformed to the normal cdf
  # scale.
  param_X = qnorm(cumsum(props[-length(props)]))
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
    two_step_fit = fit_copula_submodel_OrdCont(X, Y, copula_family, marginal_Y, start_Y, start_copula, K = K, twostep = TRUE)
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
    marginal_Y = marginal_cont_constructor(marginal_Y, est_Y),
    copula_family = copula_family,
    data = data.frame(X = X, Y = Y),
    names_XY = names_XY
  )
  return(submodel_fit)
}


#' Loglikelihood function for ordinal-continuous copula model
#'
#' [ordinal_continuous_loglik()] computes the observed-data loglikelihood for a
#' bivariate copula model with a continuous and an ordinal endpoint. The model
#' is based on a latent variable representation of the ordinal endpoint.
#'
#' @details
#'
#' ## Vine Copula Model for Ordinal Endpoints
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
#' ## Observed-Data Likelihood
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
#' @param K Number of categories in `X`.
#' @param marginal_Y List with the following five elements (in order):
#' * Density function with first argument `x` and second argument `para` the parameter
#' vector for this distribution.
#' * Distribution function with first argument `x` and second argument `para`.
#' * Inverse distribution function with first argument `p` and second argument `para`.
#' * The number of elements in `para`.
#' * Starting values for `para`.
#' @inheritParams loglik_copula_scale
#'
#' @return (numeric) loglikelihood value evaluated in `para`.
ordinal_continuous_loglik <- function(para, X, Y, copula_family, marginal_Y, K, return_sum = TRUE){
  # Number of independent cut-points for the distribution of X.
  p1 = K - 1
  # Parameters for the distribution of X.
  para_X = para[1:p1]
  # Number of parameters for the distribution of Y.
  p2 = marginal_Y[[4]]
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
  loglik = log(lik)
  if (return_sum) loglik = sum(loglik)

  return(loglik)
}

conditional_mean_copula_OrdCont = function(fitted_submodel, grid) {
  para = coef(fitted_submodel$ml_fit)
  K = length(unique(fitted_submodel$data$X))
  pdf_y = fitted_submodel$marginal_Y$pdf
  marginal_Y = attr(fitted_submodel$marginal_Y, "constructor")
  # Compute the expected value through numerical integration. First, define a
  # helper function that computes the bivariate estimated density.
  dens_joint = function(x, y) {
    dens = exp(
      ordinal_continuous_loglik(
        para = para,
        X = x,
        Y = y,
        copula_family = fitted_submodel$copula_family,
        marginal_Y = marginal_Y,
        K = K,
        return_sum = FALSE
      )
    )
    dens[is.nan(dens)] = 0
    return(dens)
  }

  conditional_mean = sapply(grid,
                            function(y) {
                              sum(dens_joint(1:K, rep(y, K)) * 1:K) / pdf_y(y)
                            })
  return(conditional_mean)
}

