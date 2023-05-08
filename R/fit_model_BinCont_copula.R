#' Fit copula model for binary true endpoint and continuous surrogate endpoint
#'
#' The function `fit_copula_model_BinCont()` fits the copula model for a
#' continuous surrogate endpoint and binary true endpoint. Because the bivariate
#' distributions of the surrogate-true endpoint pairs are functionally
#' independent across treatment groups, a bivariate distribution is fitted in
#' each treatment group separately.
#'
#' @param data
#' @inheritParams fit_model_SurvSurv
#'
#' @return
#' @export
#'
#' @examples
fit_copula_model_BinCont = function(data,
                             copula_family,
                             marginal_surrogate,
                             fitted_model = NULL,
                             hessian = TRUE,
                             maxit = 500) {
  # Column names are added to make the interpretation of the further code
  # easier. surr refers to the surrogate, true refers to the true endpoint
  colnames(data) = c("surr", "true", "Treat")

  #choose correct log-likelihood function
  #starting value for the association parameter is obtained by
  #estimating the copula parameter through kendall's tau, ignoring censoring
  tau_0 = cor(data$surr[data$Treat == 0], data$true[data$Treat == 0],
              method = "kendall")
  tau_1 = cor(data$surr[data$Treat == 1], data$true[data$Treat == 1],
              method = "kendall")

  if(copula_family == "gaussian"){
    inv_tau_0 = iTau(copula = ellipCopula(family = "normal"),
                     tau = tau_0)
    inv_tau_0 = log(1 + inv_tau_0) - log(1 - inv_tau_0)
    inv_tau_1 = iTau(copula = ellipCopula(family = "normal"),
                     tau = tau_1)
    inv_tau_1 = log(1 + inv_tau_1) - log(1 - inv_tau_1)
  }
  else if(copula_family == "clayton"){
    inv_tau_0 = iTau(copula = claytonCopula(),
                     tau = tau_0)
    inv_tau_1 = iTau(copula = claytonCopula(),
                     tau = tau_1)
  }
  else if(copula_family == "frank"){
    inv_tau_0 = iTau(copula = frankCopula(),
                     tau = tau_0)
    inv_tau_1 = iTau(copula = frankCopula(),
                     tau = tau_1)
  }
  else if(copula_family == "gumbel"){
    inv_tau_0 = iTau(copula = gumbelCopula(),
                     tau = tau_0)
    inv_tau_1 = iTau(copula = gumbelCopula(),
                     tau = tau_1)
  }

  #use partly data based starting values
  inits_0 = c(mean(data$surr[data$Treat == 0]), 0,
              sd(data$surr[data$Treat == 0]),
              abs(inv_tau_0) + 0.1)
  inits_1 = c(mean(data$surr[data$Treat == 1]), 0,
              sd(data$surr[data$Treat == 1]),
              abs(inv_tau_1) + 0.1)

  fit_0 = optim(par = inits_0, fn = binary_continuous_loglik, method = "BFGS",
                X = data$surr[data$Treat == 0], Y = data$true[data$Treat == 0],
                marginal_true = marginal_true, marginal_surrogate = marginal_surrogate,
                copula_family = copula_family,
                control = list(maxit = maxit, fnscale = -1, reltol = 1e-8,
                               ndeps = rep(1e-5, 4)),
                hessian = TRUE)
  fit_1 = optim(par = inits_1, fn = binary_continuous_loglik, method = "BFGS",
                X = data$surr[data$Treat == 1], Y = data$true[data$Treat == 1],
                marginal_true = marginal_true, marginal_surrogate = marginal_surrogate,
                copula_family = copula_family,
                control = list(maxit = maxit, fnscale = -1, reltol = 1e-8,
                               ndeps = rep(1e-5, 4)),
                hessian = TRUE)

  #return an S3 object
  return(new_vine_copula_bc_fit(fit_0, fit_1, copula_family,
                                marginal_true, marginal_surrogate))
}

fit_copula_submodel_BinCont = function(X, Y, copula_family, marginal_surrogate) {
  # Determine starting values
  starting_values = BinCont_starting_values(X, Y, copula_family, marginal_surrogate)
  # Maximize likelihood
  ml_fit = maxLik::maxLik(
    logLik = binary_continuous_loglik,
    start = starting_values,
    X = X,
    Y = Y,
    copula_family = copula_family,
    marginal_surrogate = marginal_surrogate
  )
}



BinCont_starting_values = function(X, Y, copula_family, marginal_surrogate){
  # The starting value for the association parameter is obtained by estimating
  # the copula parameter through Kendall's tau, ignoring censoring. The
  # estimated Kendall's tau is then converted to the copula parameter scale.
  tau = cor(X, Y, method = "kendall")

  # Kendall's tau is converted to the copula parameter scale.
  if(copula_family == "gaussian"){
    inv_tau = iTau(copula = ellipCopula(family = "normal"),
                   tau = tau)
  }
  else if(copula_family == "clayton"){
    inv_tau = iTau(copula = claytonCopula(),
                     tau = tau)
  }
  else if(copula_family == "frank"){
    inv_tau = iTau(copula = frankCopula(),
                     tau = tau)
  }
  else if(copula_family == "gumbel"){
    inv_tau = iTau(copula = gumbelCopula(),
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
#' * `para[1]`: mean (location) parameter for surrogate distribution
#' * `para[2]`: mean parameter for latent true endpoint distribution
#' * `para[3]`: additional parameter for surrogate distribution
#' * `para[4]`: copula parameter
#' @param X First variable (continuous)
#' @param Y Second variable (binary, $0$ or $1$)
#' @param marginal_surrogate Marginal distribution for the surrogate
#' @inheritParams loglik_copula_scale
#'
#' @return
binary_continuous_loglik <- function(para, X, Y, copula_family, marginal_surrogate){
  # Vector with mean parameters of the marginal distribution.
  mean_S = para[1]
  mean_T = para[2]
  # Vector of additional parameters for marginal distributions of S.
  extra_par_S = para[3]
  # Vector of copula parameters.
  theta = para[4]

  # Construct marginal distribution and density functions.
  cdf_X = cdf_fun(mean = mean_S,
                  extra_par = extra_par_S,
                  family = marginal_surrogate)
  pdf_X = pdf_fun(mean = mean_S,
                  extra_par = extra_par_S,
                  family = marginal_surrogate)

  cdf_Y = cdf_fun(mean = mean_T,
                  extra_par = 1,
                  family = "normal")
  pdf_Y = pdf_fun(mean = mean_T,
                  extra_par = 1,
                  family = "normal")

  # Transform true endpoint variable to left/right-censoring indicator. d2 = 1
  # indicates right-censoring and d2 = -1 indicates left-censoring.
  d2 = (Y * 2) - 1

  loglik = log_likelihood_copula_model(
    theta = theta,
    X = X,
    Y = rep(0, length(x)),
    d1 = rep(0, length(x)),
    d2 = d2,
    copula = copula_family,
    cdf_X = cdf_X,
    cdf_Y = cdf_Y,
    pdf_X = pdf_X,
    pdf_Y = pdf_Y
  )

  return(loglik)
}



#' Title
#'
#' @param u
#' @param v
#' @param copula_par
#' @param family
#'
#' @return
#' @export
#'
#' @examples
partial_deriv_copula = function(u, v, copula_par, family){
  #This function return the value for the partial derivative of the copula
  #wrt v, evaluated in (u,v)

  if (family == "gaussian"){
    dC = ifelse(
      pmin(u,v) == 0,
      0,
      {
        correlation_scale = (exp(copula_par) - 1)/(exp(copula_par) + 1)
        numerator = psi_dot_v(u, v, correlation_scale)
        denominator = dnorm(qnorm(v))
        dC = numerator/denominator
      }
    )
  }
  else if (family == "clayton"){
    dC = ifelse(
      pmin(u,v) == 0,
      0,
      {
        C = (u^(-copula_par) + v^(-copula_par) - 1)^(-1/copula_par)
        dC = (C^(copula_par + 1)) / (v^(copula_par + 1))
      }
    )
  }
  else if (family == "frank"){
    dC = ifelse(
      pmin(u,v) == 0,
      0,
      {
        C = (-1/copula_par) *
          log(
            (1/(1 - exp(-copula_par))) *
              ((1 - exp(-copula_par)) - (1 - exp(-copula_par*u))*(1 - exp(-copula_par*v)))
          )
        (1 - exp(copula_par*C))/(1 - exp(copula_par*v))
      }
    )
  }
  else if (family == "gumbel"){
    dC = ifelse(
      pmin(u,v) == 0,
      0,
      {
        C = exp(-( ((-log(u))^copula_par) + ((-log(v))^copula_par) ))
        dC = C*(-log(v))^(copula_par - 1) / v*(-log(C))^(copula_par - 1)
      }
    )

  }
  return(ifelse(is.nan(dC), 0, dC))
}



# define S3 object for fitted model
new_vine_copula_bc_fit = function(fit_0, fit_1, copula_family,
                                  marginal_true, marginal_surrogate){
  structure(
    .Data = list(
      parameters0 = fit_0$par,
      parameters1 = fit_1$par,
      hessian0 = fit_0$hessian,
      hessian1 = fit_1$hessian,
      log_lik0 = fit_0$value,
      log_lik1 = fit_1$value,
      copula_family = copula_family,
      marginal_true = marginal_true,
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
