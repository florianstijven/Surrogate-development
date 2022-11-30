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
fit_model_SurvSurv = function(data, copula_family, nknots = 2,
                     fitted_model = NULL, hessian = TRUE,
                     maxit = 500){
  #colnames are added to make the intrepretation of the furthe code easier
  #Pfs refers to the surrogate, Surv refers to the true endpoint
  colnames(data) = c("Pfs", "Surv", "Treat", "PfsInd", "SurvInd")
  #choose correct log-likelihood function
  #starting value for the association parameter is obtained by
  #estimating the copula parameter through kendall's tau, ignoring censoring
  tau_0 = cor(data$Pfs[data$Treat == 0], data$Surv[data$Treat == 0],
              method = "kendall")
  tau_1 = cor(data$Pfs[data$Treat == 1], data$Surv[data$Treat == 1],
              method = "kendall")
  if (copula_family == "gaussian") {
    log_lik_fn = normal_loglik
    inv_tau_0 = copula::iTau(copula = copula::ellipCopula(family = "normal"),
                             tau = tau_0)
    inv_tau_0 = log(1 + inv_tau_0) - log(1 - inv_tau_0)
    inv_tau_1 = copula::iTau(copula = copula::ellipCopula(family = "normal"),
                             tau = tau_1)
    inv_tau_1 = log(1 + inv_tau_1) - log(1 - inv_tau_1)
  }
  else if (copula_family == "clayton") {
    log_lik_fn = clayton_loglik
    inv_tau_0 = copula::iTau(copula = copula::claytonCopula(),
                             tau = tau_0)
    inv_tau_1 = copula::iTau(copula = copula::claytonCopula(),
                             tau = tau_1)
  }
  else if (copula_family == "frank") {
    log_lik_fn = frank_loglik
    inv_tau_0 = copula::iTau(copula = copula::frankCopula(),
                             tau = tau_0)
    inv_tau_1 = copula::iTau(copula = copula::frankCopula(),
                             tau = tau_1)
  }
  else if (copula_family == "gumbel") {
    log_lik_fn = gumbel_loglik
    inv_tau_0 = copula::iTau(copula = copula::gumbelCopula(),
                             tau = tau_0)
    inv_tau_1 = copula::iTau(copula = copula::gumbelCopula(),
                             tau = tau_1)
  }
  #fit univariate models to provide starting values
  #if fitted model is provided to provide starting values
  #the knots from that fitted model are used
  if (is.null(fitted_model)) {
    fit_s0 = flexsurv::flexsurvspline(
      formula = survival::Surv(Pfs, PfsInd) ~ 1,
      data = data,
      subset = data$Treat == 0,
      k = nknots,
      scale = "hazard"
    )
    fit_t0 = flexsurv::flexsurvspline(
      formula = survival::Surv(Surv, SurvInd) ~ 1,
      data = data,
      subset = data$Treat == 0,
      k = nknots,
      scale = "hazard"
    )
    fit_s1 = flexsurv::flexsurvspline(
      formula = survival::Surv(Pfs, PfsInd) ~ 1,
      data = data,
      subset = data$Treat == 1,
      k = nknots,
      scale = "hazard"
    )
    fit_t1 = flexsurv::flexsurvspline(
      formula = survival::Surv(Surv, SurvInd) ~ 1,
      data = data,
      subset = data$Treat == 1,
      k = nknots,
      scale = "hazard"
    )
  }
  else{
    fit_s0 = flexsurv::flexsurvspline(
      formula = survival::Surv(Pfs, PfsInd) ~ 1,
      data = data,
      subset = data$Treat == 0,
      k = nknots,
      scale = "hazard",
      knots = fitted_model$knots0[2:(nknots + 1)],
      bknots = fitted_model$knots0[c(1, nknots + 2)]
    )
    fit_t0 = flexsurv::flexsurvspline(
      formula = survival::Surv(Surv, SurvInd) ~ 1,
      data = data,
      subset = data$Treat == 0,
      k = nknots,
      scale = "hazard",
      knots = fitted_model$knott0[2:(nknots + 1)],
      bknots = fitted_model$knott0[c(1, nknots + 2)]
    )
    fit_s1 = flexsurv::flexsurvspline(
      formula = survival::Surv(Pfs, PfsInd) ~ 1,
      data = data,
      subset = data$Treat == 1,
      k = nknots,
      scale = "hazard",
      knots = fitted_model$knots1[2:(nknots + 1)],
      bknots = fitted_model$knots1[c(1, nknots + 2)]
    )
    fit_t1 = flexsurv::flexsurvspline(
      formula = survival::Surv(Surv, SurvInd) ~ 1,
      data = data,
      subset = data$Treat == 1,
      k = nknots,
      scale = "hazard",
      knots = fitted_model$knott1[2:(nknots + 1)],
      bknots = fitted_model$knott1[c(1, nknots + 2)]
    )
  }


  #allow for custom starting values
  #if inits_0/1 is provided as an argument, use that value
  if (is.null(fitted_model)) {
    inits_0 = c(fit_s0$coefficients, fit_t0$coefficients, inv_tau_0)
    inits_1 = c(fit_s1$coefficients, fit_t1$coefficients, inv_tau_1)
  }
  else{
    inits_0 = fitted_model$parameters0
    inits_1 = fitted_model$parameters1
  }

  fit_0 = optim(
    par = inits_0,
    fn = log_lik_fn,
    method = "BFGS",
    X = data$Pfs[data$Treat == 0],
    Y = data$Surv[data$Treat == 0],
    d1 = data$PfsInd[data$Treat == 0],
    d2 = data$SurvInd[data$Treat == 0],
    k = nknots,
    knotsx = fit_s0$knots,
    knotsy = fit_t0$knots,
    control = list(
      maxit = maxit,
      fnscale = -1,
      reltol = 1e-8,
      ndeps = rep(1e-5, 2 * (nknots + 2) + 1)
    ),
    hessian = hessian
  )
  fit_1 = optim(
    par = inits_1,
    fn = log_lik_fn,
    method = "BFGS",
    X = data$Pfs[data$Treat == 1],
    Y = data$Surv[data$Treat == 1],
    d1 = data$PfsInd[data$Treat == 1],
    d2 = data$SurvInd[data$Treat == 1],
    k = nknots,
    knotsx = fit_s1$knots,
    knotsy = fit_t1$knots,
    control = list(
      maxit = maxit,
      fnscale = -1,
      reltol = 1e-8,
      ndeps = rep(1e-5, 2 * (nknots + 2) + 1)
    ),
    hessian = hessian
  )

  return(
    new_vine_copula_ss_fit(
      fit_0 = fit_0,
      fit_1 = fit_1,
      copula_family = copula_family,
      knots0 = fit_s0$knots,
      knots1 = fit_s1$knots,
      knott0 = fit_t0$knots,
      knott1 = fit_t1$knots
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
      parameters0 = fit_0$par,
      parameters1 = fit_1$par,
      hessian0 = fit_0$hessian,
      hessian1 = fit_1$hessian,
      log_lik0 = fit_0$value,
      log_lik1 = fit_1$value,
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
