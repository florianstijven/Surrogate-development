#' Fit Survival-Survival model
#'
#' The function `fit_model_SurvSurv()` fits the copula model for time-to-event
#' surrogate and true endpoints (Stijven et al., 2022). Because the bivariate
#' distributions of the surrogate-true endpoint pairs are functionally
#' independent across treatment groups, a bivariate distribution is fitted in
#' each treatment group separately. The marginal distributions are based on the Royston-Parmar
#' survival model (Royston and Parmar, 2002).
#'
#' @details
#'
#' # Model
#'
#' In the causal-inference approach to evaluating surrogate endpoints, the first
#' step is to estimate the joint distribution of the relevant potential
#' outcomes. Let \eqn{(T_0, S_0, S_1, T_1)'} denote the vector of potential
#' outcomes where \eqn{(S_k, T_k)'} is the pair of potential outcomes under
#' treatment \eqn{Z = k}. \eqn{T} refers to the true endpoint, e.g., overall
#' survival. \eqn{S} refers to the composite surrogate endpoint, e.g.,
#' progression-free-survival. Because \eqn{S} is
#' usually a composite endpoint with death as possible event, modeling
#' difficulties arise because \eqn{Pr(S_k = T_k) > 0}.
#'
#' Due to difficulties in modeling the composite surrogate and the true endpoint
#' jointly, the time-to-surrogate event (\eqn{\tilde{S}}) is modeled instead of
#' the time-to-composite surrogate event (\eqn{S}). Using this new variable,
#' \eqn{\tilde{S}}, a D-vine copula model is proposed for \eqn{(T_0,
#' \tilde{S}_0, \tilde{S}_1, T_1)'} in Stijven et al. (2022). However, only the
#' following bivariate distributions are identifiable \eqn{(T_k, \tilde{S}_k)'}
#' for \eqn{k=0,1}. The margins in these bivariate distributions are based on
#' the Royston-Parmar survival model (Roystona and Parmar, 2002). The
#' association is modeled through two copulas of the same parametric form, but
#' with unique copula parameters.
#'
#' Two modelling choices are made before estimating the two bivariate
#' distributions described in the previous paragraph:
#' * The number of internal knots for the Royston-Parmar survival models. This
#' is specified through the `nknots` argument. The number of knots is assumed to
#' be equal across the four margins.
#' * The parametric family of the bivariate copulas. The parametric family is
#' assumed to be equal across treatment groups. This choice is specified through
#' the `copula_family` argument.
#'
#' # Data Format
#'
#' The data frame should have the semi-competing risks format. The columns must
#' be ordered as follows:
#'
#' * time to surrogate event, true event, or independent censoring; whichever
#' comes first
#' * time to true event, or independent censoring; whichever comes first
#' * treatment indicator: 0 or 1
#' * surrogate event indicator: 1 if surrogate event is observed, 0 otherwise
#' * true event indicator: 1 if true event is observed, 0 otherwise
#'
#' Note that according to the methodology in Stijven et al. (2022), the
#' surrogate event must not be the composite event. For example, when the
#' surrogacy of progression-free survival for overall survival is evaluated. The
#' surrogate event is progression, but not the composite event of progression or
#' death.
#'
#' @references Stijven, F., Alonso, a., Molenberghs, G., Van Der Elst, W., Van
#' Keilegom, I. (2022). An information-theoretic approach to the evaluation of
#' time-to-event surrogates for time-to-event true endpoints based on causal
#' inference.
#'
#' Royston, P., & Parmar, M. K. (2002). Flexible parametric proportional-hazards
#' and proportional-odds models for censored survival data, with application to
#' prognostic modelling and estimation of treatment effects. Statistics in
#' medicine, 21(15), 2175-2197.
#'
#'
#' @param data A data frame in the correct format (See details).
#' @param copula_family One of the following parametric copula families:
#'   `"clayton"`, `"frank"`, `"gaussian"`, or `"gumbel"`.
#' @param nknots Number of internal knots for the Royston-Parmar survival model.
#' @param fitted_model Fitted model from which initial values are extracted. If
#'   `NULL` (default), standard initial values are used. This option intended
#'   for when a model is repeatedly fitted, e.g., in a bootstrap.
#' @param hessian A boolean.
#' * `TRUE` (default): Hessian is computed
#' * `FALSE`: Hessian is not computed. This can save a small amount of time. This
#' can be useful when a model is repeatedly fitted, e.g., in a bootstrap.
#' @param maxit Maximum number of iterations for the numeric optimization,
#'   defaults to 500.
#'
#' @return Returns an S3 object that can be used to perform the sensitivity
#'   analysis with [ica_SurvSurv_sens()].
#' @export
#'
#' @author Florian Stijven
#'
#' @seealso [marginal_gof_scr()], [ica_SurvSurv_sens()]
#'
#' @examples
#' if(require(Surrogate)) {
#'   data("Ovarian")
#'   #For simplicity, data is not recoded to semi-competing risks format, but is
#'   #left in the composite event format.
#'   data = data.frame(Ovarian$Pfs,
#'                     Ovarian$Surv,
#'                     Ovarian$Treat,
#'                     Ovarian$PfsInd,
#'                     Ovarian$SurvInd)
#'   Surrogate::fit_model_SurvSurv(data = data,
#'                                 copula_family = "clayton",
#'                                 nknots = 1)
#' }
#'
#' @importFrom stats cor optim
fit_model_SurvSurv = function(data,
                              copula_family,
                              nknots = 2,
                              fitted_model = NULL,
                              method = "BFGS",
                              maxit = 500) {
  # Column names are added to make the intrepretation of the further code
  # easier. Pfs refers to the surrogate, Surv refers to the true endpoint.
  colnames(data) = c("Pfs", "Surv", "Treat", "PfsInd", "SurvInd")

  # Split original dataset into two data sets, one for each treatment group.
  data0 = data[data$Treat == 0, ]
  data1 = data[data$Treat == 1, ]

  # If required, the full maximum likelihood estimator is used where for which
  # the twostep estimator provides starting values. Else, the twostep estimator
  # gives the final estimate.
  twostep_fit0 = twostep_SurvSurv(
    X = data0$Pfs,
    delta_X = data0$PfsInd,
    Y = data0$Surv,
    delta_Y = data0$SurvInd,
    copula_family = copula_family,
    n_knots = nknots,
    method = method
  )
  twostep_fit1 = twostep_SurvSurv(
    X = data1$Pfs,
    delta_X = data1$PfsInd,
    Y = data1$Surv,
    delta_Y = data1$SurvInd,
    copula_family = copula_family,
    n_knots = nknots,
    method = method
  )

  # knots
  knots0 = twostep_fit0$marginal_S_dist$knots
  knott0 = twostep_fit0$marginal_T_dist$knots
  knots1 = twostep_fit1$marginal_S_dist$knots
  knott1 = twostep_fit1$marginal_T_dist$knots

  #allow for custom starting values
  #if inits_0/1 is provided as an argument, use that value
  if (is.null(fitted_model)) {
    inits_0 = coef(twostep_fit0$ml_fit)
    inits_1 = coef(twostep_fit1$ml_fit)
  }
  else{
    inits_0 = fitted_model$parameters0
    inits_1 = fitted_model$parameters1
  }

  # Maximize likelihood
  log_lik_function_0 = function(para) {
    temp_fun = survival_survival_loglik(
      para = para,
      X = data0$Pfs,
      delta_X = data0$PfsInd,
      Y = data0$Surv,
      delta_Y = data0$SurvInd,
      copula_family = copula_family,
      knotsx = knots0,
      knotsy = knott0
    )
    return(temp_fun)
  }

  # Starting values come from estimated marginal distribution parameters and an
  # educated guess for the copula parameter.
  suppressWarnings({
    ml_fit_0 = maxLik::maxLik(
      logLik = log_lik_function_0,
      start = inits_0,
      method = method,
    )
  })

  # Maximize likelihood
  log_lik_function_1 = function(para) {
    temp_fun = survival_survival_loglik(
      para = para,
      X = data1$Pfs,
      delta_X = data1$PfsInd,
      Y = data1$Surv,
      delta_Y = data1$SurvInd,
      copula_family = copula_family,
      knotsx = knots1,
      knotsy = knott1
    )
    return(temp_fun)
  }

  # Starting values come from estimated marginal distribution parameters and an
  # educated guess for the copula parameter.
  suppressWarnings({
    ml_fit_1 = maxLik::maxLik(
      logLik = log_lik_function_1,
      start = inits_1,
      method = method,
    )
  })


  return(
    new_vine_copula_ss_fit(
      fit_0 = ml_fit_0,
      fit_1 = ml_fit_1,
      copula_family = copula_family,
      knots0 = knots0,
      knots1 = knots1,
      knott0 = knott0,
      knott1 = knott1
    )
  )
}


#' Constructor for vine copula model
#'
#' @param fit_0 Estimated parameters in the control group.
#' @param fit_1 Estimated parameters in the experimental group
#' @param copula_family Parametric copula family
#' @param knots0 placement of knots for Royston-Parmar model
#' @param knots1 placement of knots for Royston-Parmar model
#' @param knott0 placement of knots for Royston-Parmar model
#' @param knott1 placement of knots for Royston-Parmar model
#'
#' @return S3 object
#'
#' @examples
#' #should not be used be the user
new_vine_copula_ss_fit = function(fit_0, fit_1, copula_family,
                                  knots0, knots1, knott0, knott1){
  structure(
    .Data = list(
      fit_0 = fit_0,
      fit_1 = fit_1,
      copula_family = copula_family,
      knots0 = knots0,
      knots1 = knots1,
      knott0 = knott0,
      knott1 = knott1
    ),
    class = "vine_copula_bc_fit"
  )
}


#' Goodness of fit information for survival-survival model
#'
#' This function returns several goodness-of-fit measures for a model fitted by
#' [fit_model_SurvSurv()]. These are primarily intended for model selection.
#'
#' The following goodness-of-fit measures are returned in a named vector:
#'
#' * `tau_0` and `tau_1`: (latent) value for Kendall's tau in the estimated
#' model.
#' * `log_lik`: the maximized log-likelihood value.
#' * `AIC`: the Aikaike information criterion of the fitted model.
#'
#' @param fitted_model returned value from [fit_model_SurvSurv()].
#'
#' @return a named vector containing the goodness-of-fit measures
#' @export
#'
#' @examples
#' library(Surrogate)
#' data("Ovarian")
#' #For simplicity, data is not recoded to semi-competing risks format, but is
#' #left in the composite event format.
#' data = data.frame(
#'   Ovarian$Pfs,
#'   Ovarian$Surv,
#'   Ovarian$Treat,
#'   Ovarian$PfsInd,
#'   Ovarian$SurvInd
#' )
#' ovarian_fitted =
#'     fit_model_SurvSurv(data = data,
#'                        copula_family = "clayton",
#'                        nknots = 1)
#' model_fit_measures(ovarian_fitted)
model_fit_measures = function(fitted_model){
  #number of internal knots
  nknots = length(fitted_model$knots0) - 2
  #total number of parameters
  n_parameters = length(fitted_model$parameters0) + length(fitted_model$parameters1)

  #get fitted copula parameters
  copula_par0 = fitted_model$parameters0[2*(nknots + 2) + 1]
  copula_par1 = fitted_model$parameters1[2*(nknots + 2) + 1]
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

#' Fit survival-survival copula submodel with two-step estimator
#'
#' The `twostep_SurvSurv()` function fits the copula (sub)model for a
#' time-to-event surrogate and true endpoint with a two-step estimator. In the
#' first step, the marginal distribution parameters are estimated through
#' maximum likelihood. In the second step, the copula parameter is estimate
#' while holding the marginal distribution parameters fixed.
#'
#' @param X (numeric) Possibly right-censored time-to-surrogate event
#' @param Y (numeric) Possibly right-censored time-to-true endpoint event
#' @param marginal_surrogate_estimator Not yet implemented
#' @param method Optimization algorithm for maximizing the objective function.
#'   For all options, see `?maxLik::maxLik`. Defaults to `"BFGRS"`.
#'
#' @inheritParams loglik_copula_scale
#' @inheritParams survival_survival_loglik
#'
#' @return A list with three elements:
#' * ml_fit: object of class `maxLik::maxLik` that contains the estimated copula
#'  model.
#' * marginal_S_dist: object of class `fitdistrplus::fitdist` that represents the
#'  marginal surrogate distribution.
#' * copula_family: string that indicates the copula family
#' @export
twostep_SurvSurv = function(X,
                            delta_X,
                            Y,
                            delta_Y,
                            copula_family,
                            n_knots,
                            method = "BFGS") {
  # Compute starting values; This function already fits the marginal survival
  # models by maximum likelihood.
  starting_values_list = SurvSurv_starting_values(X, delta_X, Y, delta_Y, copula_family, n_knots)
  marg_coef_x = coef(starting_values_list$fit_s)
  marg_coef_y = coef(starting_values_list$fit_t)
  # The placement of knots of the above fitted survival models are further used.
  knotsx = starting_values_list$fit_s$knots
  knotsy = starting_values_list$fit_t$knots

  data = data.frame(
    X = X,
    delta_X = delta_X,
    Y = Y,
    delta_Y = delta_Y
  )
  # Maximize likelihood
  log_lik_function = function(para) {
    temp_fun = survival_survival_loglik(
      para = para,
      X = X,
      delta_X = delta_X,
      Y = Y,
      delta_Y = delta_Y,
      copula_family = copula_family,
      knotsx = knotsx,
      knotsy = knotsy
    )
    return(temp_fun)
  }

  # Starting values come from estimated marginal distribution parameters and an
  # educated guess for the copula parameter.
  copula_start = starting_values_list$inv_tau
  start = c(marg_coef_x, marg_coef_y, copula_start)
  suppressWarnings({
    ml_fit = maxLik::maxLik(
      logLik = log_lik_function,
      start = start,
      method = method,
      fixed = 1:(2 * n_knots + 4)
    )
  })


  # Return a list with the marginal surrogate distribution, fit object for
  # copula parameter, and element to indicate the copula family.
  submodel_fit = list(
    ml_fit = ml_fit,
    marginal_S_dist = starting_values_list$fit_s,
    marginal_T_dist = starting_values_list$fit_t,
    copula_family = copula_family
  )
  return(submodel_fit)
}

#' Loglikelihood function for survival-continuous copula model
#'
#'
#' @return
survival_survival_loglik =  function(para,
                                     X,
                                     delta_X,
                                     Y,
                                     delta_Y,
                                     copula_family,
                                     knotsx,
                                     knotsy)
{
  k = length(knotsx) - 2
  switch(
    copula_family,
    "clayton" = clayton_loglik(para, X, Y, delta_X, delta_Y, k, knotsx, knotsy),
    "frank" = frank_loglik(para, X, Y, delta_X, delta_Y, k, knotsx, knotsy),
    "gumbel" = gumbel_loglik(para, X, Y, delta_X, delta_Y, k, knotsx, knotsy),
    "gaussian" = normal_loglik(para, X, Y, delta_X, delta_Y, k, knotsx, knotsy),
  )
}

SurvSurv_starting_values = function(X, delta_X, Y, delta_Y, copula_family, nknots){
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

  data = data.frame(
    X = X,
    delta_X = delta_X,
    Y = Y,
    delta_Y = delta_Y
  )
  fit_s = flexsurv::flexsurvspline(
    formula = survival::Surv(X, delta_X) ~ 1,
    data = data,
    k = nknots,
    scale = "hazard"
  )
  fit_t = flexsurv::flexsurvspline(
    formula = survival::Surv(Y, delta_Y) ~ 1,
    data = data,
    k = nknots,
    scale = "hazard"
  )

  # Return vector of informed starting values.
  starting_values = list(fit_s = fit_s, fit_t = fit_t, inv_tau = inv_tau)
  return(starting_values)
}

