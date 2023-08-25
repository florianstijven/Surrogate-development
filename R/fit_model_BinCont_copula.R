#' Fit copula model for binary true endpoint and continuous surrogate endpoint
#'
#' The function [fit_copula_model_BinCont()] fits the copula model for a
#' continuous surrogate endpoint and binary true endpoint. Because the bivariate
#' distributions of the surrogate-true endpoint pairs are functionally
#' independent across treatment groups, a bivariate distribution is fitted in
#' each treatment group separately.
#'
#' @param twostep (boolean) if `TRUE`, the two step estimator implemented in
#' [twostep_BinCont()] is used for estimation.
#'
#' @inheritParams fit_model_SurvSurv
#' @inheritParams binary_continuous_loglik
#' @inheritParams twostep_BinCont
#'
#' @return WIP
#' @export
fit_copula_model_BinCont = function(data,
                             copula_family,
                             marginal_surrogate,
                             marginal_surrogate_estimator = NULL,
                             twostep = FALSE,
                             fitted_model = NULL,
                             maxit = 500,
                             method = "BFGS") {
  # Column names are added to make the interpretation of the further code
  # easier. surr refers to the surrogate, true refers to the true endpoint.
  colnames(data) = c("surr", "true", "Treat")

  # Split original dataset into two data sets, one for each treatment group.
  data0 = data[data$Treat == 0, ]
  data1 = data[data$Treat == 1, ]

  # If required, the full maximum likelihood estimator is used where for which
  # the twostep estimator provides starting values. Else, the twostep estimator
  # gives the final estimate.
  if (twostep) {
    fit0 = twostep_BinCont(
      X = data0$surr,
      Y = data0$true,
      copula_family = copula_family,
      marginal_surrogate = marginal_surrogate,
      marginal_surrogate_estimator = marginal_surrogate_estimator,
      method = method
    )
    fit1 = twostep_BinCont(
      X = data1$surr,
      Y = data1$true,
      copula_family = copula_family,
      marginal_surrogate = marginal_surrogate,
      marginal_surrogate_estimator = marginal_surrogate_estimator,
      method = method
    )
  }
  else {
    fit0 = fit_copula_submodel_BinCont(
      X = data0$surr,
      Y = data0$true,
      copula_family = copula_family,
      marginal_surrogate = marginal_surrogate,
      method = method
    )
    fit1 = fit_copula_submodel_BinCont(
      X = data1$surr,
      Y = data1$true,
      copula_family = copula_family,
      marginal_surrogate = marginal_surrogate,
      method = method
    )
  }

  #return an S3 object
  return(
    new_vine_copula_bc_fit(fit0,
                           fit1,
                           copula_family,
                           marginal_surrogate)
  )
}


#' Fit binary-continuous copula submodel
#'
#' The `fit_copula_submodel_BinCont()` function fits the copula (sub)model fir a
#' continuous surrogate and binary true endpoint with maximum likelihood.
#'
#' @inheritParams twostep_BinCont
#' @inherit twostep_BinCont params return
fit_copula_submodel_BinCont = function(X,
                                       Y,
                                       copula_family,
                                       marginal_surrogate,
                                       method = "BFGS") {
  # The starting values are determined by the two-step estimator.
  twostep_fit = twostep_BinCont(X, Y, copula_family, marginal_surrogate)
  starting_values = coef(twostep_fit$ml_fit)
  # Maximize likelihood
  log_lik_function = function(para) {
    temp_fun = binary_continuous_loglik(
      para = para,
      X = X,
      Y = Y,
      copula_family = copula_family,
      marginal_surrogate = marginal_surrogate
    )
    return(temp_fun)
  }
  suppressWarnings({
    ml_fit = maxLik::maxLik(logLik = log_lik_function,
                            start = starting_values,
                            method = method)
  })


  # Construct fitdistrplus::fitdist object from maximum likelihood estimates. To
  # do this, we first extract the estimated marginal surrogate distribution
  # parameters. The names of these parameters should correspond to the marginal
  # distribution arguments.
  fix.arg = coef(ml_fit)[c(-1, -length(coef(ml_fit)))]
  names(fix.arg) = stringr::str_sub(names(fix.arg),
                                    start = 1L, end = -5L)
  marginal_S_dist = marginal_distribution(
    x = X,
    distribution = marginal_surrogate,
    fix.arg = as.list(fix.arg)
  )

  # Return a list with the marginal surrogate distribution, fit object for
  # copula parameter, and element to indicate the copula family.
  submodel_fit = list(
    ml_fit = ml_fit,
    marginal_S_dist = marginal_S_dist,
    copula_family = copula_family
  )
  return(submodel_fit)
}


#' Fit binary-continuous copula submodel with two-step estimator
#'
#' The `twostep_BinCont()` function fits the copula (sub)model fir a continuous
#' surrogate and binary true endpoint with a two-step estimator. In the first
#' step, the marginal distribution parameters are estimated through maximum
#' likelihood. In the second step, the copula parameter is estimate while
#' holding the marginal distribution parameters fixed.
#'
#' @param X (numeric) Continuous surrogate variable
#' @param Y (integer) Binary true endpoint variable (\eqn{T_k \, \in \, \{0, 1\}})
#' @param marginal_surrogate_estimator Not yet implemented
#' @param method Optimization algorithm for maximizing the objective function.
#'   For all options, see `?maxLik::maxLik`. Defaults to `"BFGRS"`.
#'
#' @inheritParams loglik_copula_scale
#' @inheritParams binary_continuous_loglik
#'
#' @return A list with three elements:
#' * ml_fit: object of class `maxLik::maxLik` that contains the estimated copula
#'  model.
#' * marginal_S_dist: object of class `fitdistrplus::fitdist` that represents the
#'  marginal surrogate distribution.
#' * copula_family: string that indicates the copula family
#' @export
twostep_BinCont = function(X,
                           Y,
                           copula_family,
                           marginal_surrogate,
                           marginal_surrogate_estimator = NULL,
                           method = "BFGS") {
  # Estimate marginal distribution parameters. For clarity, the parameters are
  # named in the vector.
  mu_T = c("mu (T)" = -1 * qnorm(mean(Y)))
  # Estimate marginal surrogate distribution.
  marginal_S = c()
  if (is.null(marginal_surrogate_estimator)) {
    # Estimate parameters for parametric marginal distribution. Again for
    # clarity, the parameters are given a clear name. The
    # marginal_distirbution() function returns a named vector, but we still add
    # an indicator to its name.
    marginal_S_dist = marginal_distribution(X, marginal_surrogate)
    marginal_S = coef(marginal_S_dist)
    names(marginal_S) = paste(names(marginal_S), "(S)")
  }
  else {
    # Estimate marginal distribution non-parametrically.
  }

  # Vector of marginal parameter estimates.
  marginal_para = c(mu_T, marginal_S)



  # Maximize likelihood
  log_lik_function = function(para) {
    temp_fun = binary_continuous_loglik(
      para = para,
      X = X,
      Y = Y,
      copula_family = copula_family,
      marginal_surrogate = marginal_surrogate
    )
    return(temp_fun)
  }

  # Starting values come from estimated marginal distribution parameters and an
  # educated guess for the copula parameter.
  copula_start = BinCont_starting_values(X, Y, copula_family, marginal_surrogate)
  copula_start = copula_start[length(copula_start)]
  unname(copula_start)
  start = c(marginal_para,
            "theta (copula)" = copula_start)
  suppressWarnings({
    ml_fit = maxLik::maxLik(
      logLik = log_lik_function,
      start = start,
      method = method,
      fixed = 1:3
    )
  })


  # Return a list with the marginal surrogate distribution, fit object for
  # copula parameter, and element to indicate the copula family.
  submodel_fit = list(
    ml_fit = ml_fit,
    marginal_S_dist = marginal_S_dist,
    copula_family = copula_family
  )
  return(submodel_fit)
}

#' Fit marginal distribution
#'
#' The `marginal_distribution()` function is a wrapper for
#' `fitdistrplus::fitdist()` that fits a univariate distribution to a data
#' vector.
#'
#' @param x (numeric) data vector
#' @param distribution Distributional family. One of the follwing:
#' * `"normal"`: normal distribution
#' * `"logistic`: logistic distribution as parameterized in `dlogis()`
#' * `"t"`: student t distribution is parameterized in `dt()`
#' * `"lognormal"`: lognormal distribution as parameterized in `dlnorm()`
#' * `"gamma"`: gamma distribution as parameterized in `dgamma()`
#' * `"weibull"`: weibull distribution as parameterized in `dweibull()`
#' @param fix.arg An optional named list giving the values of fixed parameters
#'   of the named distribution or a function of data computing (fixed) parameter
#'   values and returning a named list. Parameters with fixed value are thus NOT
#'   estimated by this maximum likelihood procedure.
#'
#' @return Object of class `fitdistrplus::fitdist` that represents the marginal
#'   surrogate distribution.
marginal_distribution = function(x, distribution, fix.arg = NULL) {
  # No fixed arguments are provided.
  if (is.null(fix.arg)) {
    fitted_dist = switch(
      distribution,
      normal = fitdistrplus::fitdist(data = x, "norm", fix.arg = fix.arg),
      logistic = fitdistrplus::fitdist(data = x, "logis", fix.arg = fix.arg),
      t = fitdistrplus::fitdist(x, "t", fix.arg = fix.arg),
      lognormal = fitdistrplus::fitdist(x, "lnorm", fix.arg = fix.arg),
      gamma = fitdistrplus::fitdist(x, "gamma", start = list(shape = 1, rate = 1), fix.arg = fix.arg),
      weibull = fitdistrplus::fitdist(x, "weibull", start = list(shape = 1, scale = 1), fix.arg = fix.arg)
    )
  }
  # Fixed arguments are provided. In this case, we give the fixed arguments as
  # starting values and use a "custom" optimisation function that does not
  # iterate. Hence, the starting values are just returned.
  else {
    fitted_dist = switch(
      distribution,
      normal = fitdistrplus::fitdist(
        data = x,
        "norm",
        start = fix.arg,
        optim.method = "SANN",
        control = list(maxit = 0)
      ),
      logistic = fitdistrplus::fitdist(
        data = x,
        "logis",
        start = fix.arg,
        optim.method = "SANN",
        control = list(maxit = 0)
      ),
      t = fitdistrplus::fitdist(
        x,
        "t",
        start = fix.arg,
        optim.method = "SANN",
        control = list(maxit = 0)
      ),
      lognormal = fitdistrplus::fitdist(
        x,
        "lnorm",
        start = fix.arg,
        optim.method = "SANN",
        control = list(maxit = 0)
      ),
      gamma = fitdistrplus::fitdist(
        x,
        "gamma",
        start = fix.arg,
        optim.method = "SANN",
        control = list(maxit = 0)
      ),
      weibull = fitdistrplus::fitdist(
        x,
        "weibull",
        start = fix.arg,
        optim.method = "SANN",
        control = list(maxit = 0)
      )
    )
  }

  return(fitted_dist)
}


BinCont_starting_values = function(X, Y, copula_family, marginal_surrogate){
  # The starting value for the association parameter is obtained by estimating
  # the copula parameter through Kendall's tau, ignoring censoring. The
  # estimated Kendall's tau is then converted to the copula parameter scale.
  tau = cor(X, Y, method = "kendall")
  # tau = 0.05
  # Kendall's tau is converted to the copula parameter scale.
  if(copula_family == "gaussian"){
    inv_tau = copula::iTau(copula = copula::ellipCopula(family = "normal"),
                   tau = tau)
  }
  else if(copula_family == "clayton"){
    inv_tau = copula::iTau(copula = copula::claytonCopula(),
                     tau = tau)
  }
  else if(copula_family == "frank"){
    inv_tau = copula::iTau(copula = copula::frankCopula(),
                     tau = tau)
  }
  else if(copula_family == "gumbel"){
    inv_tau = copula::iTau(copula = copula::gumbelCopula(),
                     tau = tau)
  }

  # Compute mean of latent normal variable for the true endpoint.
  mu_T = -1 * qnorm(mean(Y))
  # Compute mean and standard deviation of surrogate endpoint.
  mu_S = mean(X)
  sd_S = sd(X)

  # Return vector of informed starting values.
  starting_values = c(mu_S = mu_S, mu_T = mu_T, add_param_S = sd_S, copula_param = inv_tau)
  return(starting_values)
}



#' Loglikelihood function for binary-continuous copula model
#'
#' @param para Parameter vector. The parameters are ordered as follows:
#' * `para[1]`: mean parameter for latent true endpoint distribution
#' * `para[2:p]`: Parameters for surrogate distribution, more details in
#'  `?Surrogate::cdf_fun` for the specific implementations.
#' * `para[p + 1]`: copula parameter
#' @param X First variable (continuous)
#' @param Y Second variable (binary, $0$ or $1$)
#' @param marginal_surrogate Marginal distribution for the surrogate. For all
#'   available options, see `?Surrogate::cdf_fun`.
#' @inheritParams loglik_copula_scale
#'
#' @return (numeric) loglikelihood value evaluated in `para`.
binary_continuous_loglik <- function(para, X, Y, copula_family, marginal_surrogate){
  # Parameter for true endpoint distribution.
  para_T = c(para[1], 1)
  # Parameter(s) for surrogate endpoint distribution.
  para_S = para[2:(length(para) - 1)]
  # Vector of copula parameters.
  theta = para[length(para)]

  # Construct marginal distribution and density functions.
  cdf_X = cdf_fun(para = para_S,
                  family = marginal_surrogate)
  pdf_X = pdf_fun(para = para_S,
                  family = marginal_surrogate)

  cdf_Y = cdf_fun(para = para_T,
                  family = "normal")
  pdf_Y = pdf_fun(para = para_T,
                  family = "normal")

  # Transform true endpoint variable to left/right-censoring indicator. d2 = 0
  # indicates right-censoring and d2 = -1 indicates left-censoring.
  d2 = Y - 1

  loglik = log_likelihood_copula_model(
    theta = theta,
    X = X,
    Y = rep(0, length(X)),
    d1 = rep(1, length(X)),
    d2 = d2,
    copula_family = copula_family,
    cdf_X = cdf_X,
    cdf_Y = cdf_Y,
    pdf_X = pdf_X,
    pdf_Y = pdf_Y
  )

  return(loglik)
}



# define S3 object for fitted model
new_vine_copula_bc_fit = function(fit0,
                                  fit1,
                                  copula_family,
                                  marginal_surrogate) {
  structure(
    .Data = list(
      submodel0 = fit0,
      submodel1 = fit1,
      copula_family = copula_family,
      marginal_surrogate = marginal_surrogate
    ),
    class = "vine_copula_bc_fit"
  )
}

model_fit_measures = function(fitted_model){
  #total number of parameters
  n_parameters = length(fitted_model$parameters0) + length(fitted_model$parameters1)
  #get fitted copula parameters
  copula_par0 = fitted_model$parameters0[length(fitted_model$parameters0)]
  copula_par1 = fitted_model$parameters1[length(fitted_model$parameters1)]
  #convert fitted copula parameters to kendall's tau scale
  tau_0 = conversion_copula_tau(copula_par = copula_par0,
                                copula_family = fitted_model$copula_family)
  tau_1 = conversion_copula_tau(copula_par = copula_par1,
                                copula_family = fitted_model$copula_family)

  #compute total maximized log likelihood
  log_lik = fitted_model$log_lik0 + fitted_model$log_lik1
  AIC = -2*log_lik + 2*n_parameters

  return(c(tau_0 = tau_0, tau_1 = tau_1, log_lik = log_lik, AIC = AIC))
}

conversion_copula_tau = function(copula_par, copula_family){
  if(copula_family == "frank"){
    return(
      tau(copula = frankCopula(copula_par))
    )
  }
  else if(copula_family == "gaussian"){
    #convert to correct scale
    correlation_scale = (exp(copula_par) - 1)/(exp(copula_par) + 1)
    return(
      tau(copula = normalCopula(correlation_scale))
    )
  }
  else if(copula_family == "clayton"){
    return(
      tau(copula = claytonCopula(copula_par))
    )
  }
  else if(copula_family == "gumbel"){
    return(
      tau(copula = gumbelCopula(copula_par))
    )
  }
}
