#' Title
#'
#' @param data
#' @inheritParams fit_model_SurvSurv
#'
#' @return
#' @export
#'
#' @examples
fit_model_BinCont = function(data,
                             copula_family,
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

#' Title
#'
#' @param para
#' @param X
#' @param Y
#' @param copula_family
#' @param marginal_surrogate
#' @param marginal_true
#'
#' @return
#' @export
#'
#' @examples
binary_continuous_loglik <- function(para, X, Y, copula_family,
                                     marginal_surrogate, marginal_true){
  # Vector with mean parameters of the marginal distribution.
  mean_S = para[1]
  mean_T = para[2]
  # Vector of additional parameters for marginal distributions of S.
  extra_par_S = para[3]
  # Vector of copula parameters.
  copula_par = para[4]

  log_likelihood_copula_model()
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

#helper function for the gaussian copula
psi_dot_v = function(u, v, rho){
  q_u = qnorm(p = u)
  q_v = qnorm(p = v)
  dens_q_v = dnorm(q_v)
  cond_mean = rho*q_v
  cond_var = 1 - rho^2
  p_q_u = pnorm(q = q_u, mean = cond_mean, sd = sqrt(cond_var), lower.tail = TRUE)

  return(dens_q_v*p_q_u)
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
