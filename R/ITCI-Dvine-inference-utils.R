#' Confidence interval for the ICA given the unidentifiable parameters
#'
#' [Dvine_ICA_confint()] computes the confidence interval for the ICA
#' in the D-vine copula model. The unidentifiable parameters are fixed at the user
#' supplied values.
#'
#' @param alpha (numeric) `1 - alpha` is the level of the confidence interval
#' @inheritParams ICA_given_model_constructor_SurvSurv
#' @inheritParams summary_level_bootstrap_ICA
#'
#' @return (numeric) Vector with the limits of the two-sided `1 - alpha`
#'   confidence interval.
Dvine_ICA_confint = function(fitted_model,
                             alpha,
                             copula_par_unid,
                             copula_family2,
                             rotation_par_unid,
                             n_prec,
                             mutinfo_estimator = NULL,
                             composite,
                             B,
                             seed)
{
  bootstrap_replications = summary_level_bootstrap_ICA(
    fitted_model = fitted_model,
    copula_par_unid = copula_par_unid,
    copula_family2 = copula_family2,
    rotation_par_unid = rotation_par_unid,
    n_prec = n_prec,
    mutinfo_estimator = mutinfo_estimator,
    composite = composite,
    B = B,
    seed = seed
  )

  limits = stats::quantile(bootstrap_replications, prob = c(alpha / 2, 1 - (alpha / 2)))
  return(limits)
}

#' Variance of log-mutual information based on the delta method
#'
#' [delta_method_log_mutinfo()] computes the variance of the estimated log
#' mutual information, given the unidentifiable parameters.
#'
#' @param eps (numeric) Step size for finite difference in numeric
#'   differentiation
#' @inheritParams ICA_given_model_constructor_SurvSurv
#'
#' @details
#'  This function should not be used. The ICA is computed through numerical
#'  methods with a considerable error. This error is negligible in individual estimates
#'  of the ICA; however, this error easily breaks the numeric differentiation because
#'  finite differences are inflated by this error.
#'
#' @return (numeric) Variance for the estimated ICA based on the delta method,
#'   holding the unidentifiable parameters fixed at the user supplied values.
delta_method_log_mutinfo = function(fitted_model,
                                    copula_par_unid,
                                    copula_family2,
                                    rotation_par_unid,
                                    n_prec,
                                    mutinfo_estimator = NULL,
                                    composite,
                                    seed,
                                    eps = 1e-3) {
  # Number of knots
  k = length(fitted_model$knots1)

  ICA_given_model = ICA_given_model_constructor_SurvSurv(
    fitted_model = fitted_model,
    copula_par_unid = copula_par_unid,
    copula_family2 = copula_family2,
    rotation_par_unid = rotation_par_unid,
    n_prec = n_prec,
    mutinfo_estimator = mutinfo_estimator,
    composite = composite,
    seed = seed
  )

  # Compute gradient of the ICA as a function of the model parameters, evaluated
  # at the maximum likelihood estimates.
  gradient_vec = maxLik::numericGradient(
    f = function(theta)
      - 0.5 * log(1 - ICA_given_model(theta)),
    t0 = c(coef(fitted_model$fit_0), coef(fitted_model$fit_1)),
    eps = eps
  )
  # Convert gradient to a column vector.
  gradient_vec = matrix(gradient_vec, ncol = 1)

  # Compute variance following the delta method.
  zeros_matrix = matrix(rep(0, (2 * k + 1) ** 2), ncol = (2 * k + 1))
  vcov_matrix = rbind(cbind(vcov(fitted_model$fit_0), zeros_matrix),
                      cbind(zeros_matrix, vcov(fitted_model$fit_1)))

  variance_log_mutinfo = t(gradient_vec) %*% vcov_matrix %*% gradient_vec

  return(variance_log_mutinfo)
}

#' Bootstrap based on the multivariate normal sampling distribution
#'
#' [summary_level_bootstrap_ICA()] performs a parametric type of bootstrap based
#' on the estimated multivariate normal sampling distribution of the maximum
#' likelihood estimator for the (observable) D-vine copula model parameters.
#'
#' @details
#' Let \eqn{\hat{\boldsymbol{\beta}}} be the estimated identifiable parameter
#' vector, \eqn{\hat{\Sigma}} the corresponding estimated covariance matrix, and
#' \eqn{\boldsymbol{\nu}} a fixed value for the sensitivity parameter. The
#' bootstrap is then performed in the following steps
#'
#' 1. Resample the identifiable parameters from the estimated sampling distribution,
#' \deqn{\hat{\boldsymbol{\beta}}^{(b)} \sim N(\hat{\boldsymbol{\beta}}, \hat{\Sigma}).}
#' 2. For each resampled parameter vector and the fixed sensitivty parameter,
#' compute the ICA as \eqn{ICA(\hat{\boldsymbol{\beta}}^{(b)},
#' \boldsymbol{\nu})}.
#'
#'
#' @param B Number of bootstrap replications
#' @inheritParams ICA_given_model_constructor_SurvSurv
#' @inheritParams sensitivity_analysis_SurvSurv_copula
#'
#' @return (numeric) Vector of bootstrap replications for the estimated ICA.
summary_level_bootstrap_ICA = function(fitted_model,
                                       copula_par_unid,
                                       copula_family2,
                                       rotation_par_unid,
                                       n_prec,
                                       B,
                                       measure = "ICA",
                                       mutinfo_estimator = NULL,
                                       ICA_estimator = NULL,
                                       composite = FALSE,
                                       seed,
                                       restr_time = +Inf,
                                       ncores = 1) {
  withr::local_seed(seed)
  # The structure of the fitted model objects is different (for now) between the
  # survival-survival setting and the other settings.
  if (!is.null(fitted_model$endpoint_types)) {
    # Not the survival-survival setting
    # Parameter estimates
    theta_hat = c(coef(fitted_model$fit_0$ml_fit), coef(fitted_model$fit_1$ml_fit))

    # Compute covariance matrix of the sampling distribution.
    zeros_matrix = matrix(
      rep(0, length(coef(fitted_model$fit_0$ml_fit)) * length(coef(fitted_model$fit_1$ml_fit))),
      ncol = length(coef(fitted_model$fit_1$ml_fit))
    )
    vcov_matrix = rbind(
      cbind(vcov(fitted_model$fit_0$ml_fit), zeros_matrix),
      cbind(t(zeros_matrix), vcov(fitted_model$fit_1$ml_fit))
      )
  }
  else {
    # Survival-survival setting
    # Parameter estimates
    theta_hat = c(coef(fitted_model$fit_0), coef(fitted_model$fit_1))

    # Compute covariance matrix of the sampling distribution.
    zeros_matrix = matrix(rep(0, length(coef(fitted_model$fit_0)) * length(coef(fitted_model$fit_1))),
                          ncol = length(coef(fitted_model$fit_1)))
    vcov_matrix = rbind(cbind(vcov(fitted_model$fit_0), zeros_matrix),
                        cbind(t(zeros_matrix), vcov(fitted_model$fit_1)))
  }


  ICA_given_model_original = ICA_given_model_constructor(
    fitted_model = force(fitted_model),
    copula_par_unid = force(copula_par_unid),
    copula_family2 = force(copula_family2),
    rotation_par_unid = force(rotation_par_unid),
    n_prec = force(n_prec),
    measure = force(measure),
    mutinfo_estimator = force(mutinfo_estimator),
    ICA_estimator = force(ICA_estimator),
    composite = force(composite),
    seed = force(seed),
    restr_time = force(restr_time)
  )
  ICA_given_model = function(theta) {
    tryCatch(
      ICA_given_model_original(theta),
      error = function(e) {
        print(e)
        return(NA)
      }
    )
  }


  # Resample parameter estimates from the estimated multivariate normal sampling
  # distribution.
  requireNamespace("mvtnorm")
  theta_resampled = mvtnorm::rmvnorm(n = B,
                                     mean = theta_hat,
                                     sigma = vcov_matrix)

  # For the Gaussian copula, fisher's Z transformation was applied. We have to
  # backtransform to the correlation scale in that case.
  if (!is.null(fitted_model$endpoint_types)) {
    # Not the survival-survival setting
    a = length(coef(fitted_model$fit_0$ml_fit))
    b = length(coef(fitted_model$fit_1$ml_fit))

    if (fitted_model$fit_0$copula_family == "gaussian") {
      theta_resampled[, a] = (exp(2 * theta_resampled[, a]) - 1) / (exp(2 * theta_resampled[, a]) + 1)
    }
    if (fitted_model$fit_1$copula_family == "gaussian") {
      theta_resampled[, a + b] = (exp(2 * theta_resampled[, a + b]) - 1) / (exp(2 * theta_resampled[, a + b]) + 1)
    }
  }
  else {
    # Survival-survival setting
    a = length(coef(fitted_model$fit_0))
    b = length(coef(fitted_model$fit_1))

    if (fitted_model$copula_family[1] == "gaussian") {
      theta_resampled[, a] = (exp(2 * theta_resampled[, a]) - 1) / (exp(2 * theta_resampled[, a]) + 1)
    }
    if (fitted_model$copula_family[1] == "gaussian") {
      theta_resampled[, a + b] = (exp(2 * theta_resampled[, a + b]) - 1) / (exp(2 * theta_resampled[, a + b]) + 1)
    }
  }


  if (ncores > 1 & requireNamespace("parallel")) {
    cl  <- parallel::makeCluster(ncores)
    # helper function
    # surrogacy_sample_sens <- surrogacy_sample_sens
    print("Starting parallel simulations")

    # Get current search path and set the same search path in the new instances
    # of R. Usually, this would not be necessary, but if the user changed the
    # search path before running this function, there could be an issue.
    search_path = .libPaths()
    force(search_path)
    parallel::clusterExport(
      cl = cl,
      varlist = c("search_path", "ICA_given_model_original"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl = cl, expr = .libPaths(new = search_path))
    ICA_hats = parallel::parApply(
      cl = cl,
      X = theta_resampled,
      FUN = ICA_given_model,
      MARGIN = 1
    )
    print("Finishing parallel simulations")
    on.exit(expr = {
      parallel::stopCluster(cl)
      rm("cl")
    })
  }
  else {
    # Compute the ICA for the resampled parameter estimates.
    ICA_hats = apply(
      X = theta_resampled,
      FUN = ICA_given_model,
      MARGIN = 1
    )
  }

  return(ICA_hats)
}

#' Constructor for the function that returns that ICA as a function of the
#' identifiable parameters for survival-survival
#'
#' [ICA_given_model_constructor_SurvSurv()] returns a function fixes the unidentifiable
#' parameters at user-specified values and takes the identifiable parameters as
#' argument.
#'
#' @param copula_par_unid Parameter vector for the sequence of *unidentifiable*
#'   bivariate copulas that define the D-vine copula. The elements of
#'   `copula_par` correspond to \eqn{(c_{23}, c_{13;2}, c_{24;3}, c_{14;23})}.
#' @param rotation_par_unid Vector of rotation parameters for the sequence of
#'   *unidentifiable* bivariate copulas that define the D-vine copula. The elements of
#'   `rotation_par` correspond to \eqn{(c_{23}, c_{13;2},
#'   c_{24;3}, c_{14;23})}.
#' @inheritParams compute_ICA_SurvSurv
#' @inheritParams sensitivity_analysis_SurvSurv_copula
#' @inheritParams sensitivity_intervals_Dvine
#'
#' @return A function that computes the ICA as a function of the identifiable
#'   parameters. In this computation, the unidentifiable parameters are fixed at
#'   the values supplied as arguments to [ICA_given_model_constructor_SurvSurv()] or
#'   [ICA_given_model_constructor()].
ICA_given_model_constructor_SurvSurv = function(fitted_model,
                                       copula_par_unid,
                                       copula_family2,
                                       rotation_par_unid,
                                       n_prec,
                                       measure = "ICA",
                                       mutinfo_estimator,
                                       composite,
                                       seed,
                                       restr_time = +Inf) {
  # Number of knots
  ks0 = length(fitted_model$knots0)
  ks1 = length(fitted_model$knots1)
  kt0 = length(fitted_model$knott0)
  kt1 = length(fitted_model$knott1)
  n_par0 = length(coef(fitted_model$fit_0))
  n_par1 = length(coef(fitted_model$fit_1))
  # Location of the knots
  knots0 = fitted_model$knots0
  knots1 = fitted_model$knots1
  knott0 = fitted_model$knott0
  knott1 = fitted_model$knott1
  # Copula Families
  copula_family1 = fitted_model$copula_family
  # Copula rotations for the identifiable copulas
  r_12 = fitted_model$copula_rotations[1]
  r_34 = fitted_model$copula_rotations[2]
  marginal_sp_rho = TRUE
  if (measure %in% c("ICA", "sp_rho")) marginal_sp_rho = FALSE

  force(n_prec); force(mutinfo_estimator); force(copula_par_unid)
  force(copula_family2); force(rotation_par_unid); force(composite)
  force(seed); force(restr_time)

  ICA_given_model = function(theta) {
    # The first k + 1 elements of theta correspond to the parameters in fit_0,
    # the next k + 1 elements correspond to the parameters in fit_1. In a second
    # step, the parameters for each marginal distribution are extracted.
    theta0 = theta[1:n_par0]
    theta1 = theta[(n_par0 + 1):(n_par0 + n_par1)]

    gammas0 = theta0[1:ks0]
    gammas1 = theta1[1:ks1]
    gammat0 = theta0[(ks0 + 1):(ks0 + kt0)]
    gammat1 = theta1[(ks1 + 1):(ks1 + kt1)]

    # The last elements of theta_0 and theta_1 contain the copula parameters.
    c_12 = theta0[n_par0]
    c_34 = theta1[n_par1]

    # Quantile functions for the marginal distributions
    q_S0 = function(p) {
      flexsurv::qsurvspline(p = p,
                            gamma = gammas0,
                            knots = knots0)
    }
    q_S1 = function(p) {
      flexsurv::qsurvspline(p = p,
                            gamma = gammas1,
                            knots = knots1)
    }
    q_T0 = function(p) {
      flexsurv::qsurvspline(p = p,
                            gamma = gammat0,
                            knots = knott0)
    }
    q_T1 = function(p) {
      flexsurv::qsurvspline(p = p,
                            gamma = gammat1,
                            knots = knott1)
    }

    # Compute the log of the mutual information
    ICA = compute_ICA_SurvSurv(
      copula_par = c(c_12,
                     copula_par_unid[1],
                     c_34,
                     copula_par_unid[2:4]),
      rotation_par = c(r_12,
                       rotation_par_unid[1],
                       r_34,
                       rotation_par_unid[2:4]),
      copula_family1 = copula_family1,
      copula_family2 = copula_family2,
      n_prec = n_prec,
      q_S0 = q_S0,
      q_T0 = q_T0,
      q_S1 = q_S1,
      q_T1 = q_T1,
      composite = composite,
      mutinfo_estimator = mutinfo_estimator,
      seed = seed,
      marginal_sp_rho = marginal_sp_rho,
      restr_time = restr_time
    )[measure]
    return(ICA)
  }
  return(ICA_given_model)
}


#' Constructor for the function that returns that ICA as a function of the
#' identifiable parameters
#'
#' [ICA_given_model_constructor()] returns a function fixes the unidentifiable
#' parameters at user-specified values and takes the identifiable parameters as
#' argument.
#'
#' @param mutinfo_estimator Function that estimates the mutual information
#'  between the first two arguments which are numeric vectors. Defaults to
#'  `FNN::mutinfo()` with default arguments in the survival-survival setting. This
#'  argument is not used for non-survival-survival settings.
#' @param ICA_estimator Function that estimates the ICA between the first two
#'   arguments which are numeric vectors. Defaults to `NULL` which corresponds
#'   to using [estimate_ICA_ContCont()], [estimate_ICA_OrdCont()], or
#'   [estimate_ICA_OrdOrd()] (depending on the endpoint types). This argument is
#'   not used in the survival-survival setting.
#'
#' @inheritParams ICA_given_model_constructor_SurvSurv
#' @inherit ICA_given_model_constructor_SurvSurv return
ICA_given_model_constructor = function(fitted_model,
                                       copula_par_unid,
                                       copula_family2,
                                       rotation_par_unid,
                                       n_prec,
                                       measure = "ICA",
                                       mutinfo_estimator = NULL,
                                       ICA_estimator = NULL,
                                       seed,
                                       composite = NULL,
                                       restr_time = +Inf) {
  if (is.null(fitted_model$endpoint_types)) {
    return(
      ICA_given_model_constructor_SurvSurv(fitted_model,
                                           copula_par_unid,
                                           copula_family2,
                                           rotation_par_unid,
                                           n_prec,
                                           measure,
                                           mutinfo_estimator,
                                           seed,
                                           composite,
                                           restr_time)
    )
  }
  else {
    # Number of parameters corresponding to the first and second submodel.
    n_par0 = fitted_model$fit_0$marginal_X$n_param + fitted_model$fit_0$marginal_Y$n_param + 1
    n_par1 = fitted_model$fit_1$marginal_X$n_param + fitted_model$fit_1$marginal_Y$n_param + 1

    # Number of parameters corresponding to the distribution of X in the first
    # and second submodel.
    n_parX0 = fitted_model$fit_0$marginal_X$n_param
    n_parX1 = fitted_model$fit_1$marginal_X$n_param

    # Number of parameters corresponding to the distribution of Y in the first
    # and second submodel.
    n_parY0 = fitted_model$fit_0$marginal_Y$n_param
    n_parY1 = fitted_model$fit_1$marginal_Y$n_param

    # Marginal models
    marginal_X0 = attr(fitted_model$fit_0$marginal_X, "constructor")
    marginal_Y0 = attr(fitted_model$fit_0$marginal_Y, "constructor")
    marginal_X1 = attr(fitted_model$fit_1$marginal_X, "constructor")
    marginal_Y1 = attr(fitted_model$fit_1$marginal_Y, "constructor")

    copula_family1 = c(fitted_model$fit_0$copula_family,
                       fitted_model$fit_1$copula_family)
    # Copula rotations for the identifiable copulas
    r_12 = 0
    r_34 = 0
    marginal_sp_rho = TRUE
    if (measure %in% c("ICA", "sp_rho")) marginal_sp_rho = FALSE

    force(n_prec); force(ICA_estimator); force(copula_par_unid)
    force(copula_family2); force(rotation_par_unid); force(composite)
    force(seed); force(restr_time)

    # Determine which endoints X and Y correspond in fitted_model. We
    # start from X = S and Y = T as default. Otherwise, everything is switched.
    switch_XY = fitted_model$fit_0$names_XY[1] == "True"
    force(switch_XY)
    ICA_given_model = function(theta) {
      # Extract parameters for each submodel.
      theta0 = theta[1:n_par0]
      theta1 = theta[(n_par0 + 1):(n_par0 + n_par1)]

      # The last elements of theta_0 and theta_1 contain the copula parameters.
      c_12 = theta0[n_par0]
      c_34 = theta1[n_par1]

      # Quantile functions for the marginal distributions
      if (fitted_model$endpoint_types[1] == "continuous") {
        q_S0 = marginal_cont_constructor(marginal_X0, theta0[1:n_parX0])$inv_cdf
        q_S1 = marginal_cont_constructor(marginal_X1, theta1[1:n_parX1])$inv_cdf
      }
      else {
        q_S0 = marginal_ord_constructor(theta0[1:n_parX0])$inv_cdf
        q_S1 = marginal_ord_constructor(theta1[1:n_parX1])$inv_cdf
      }

      if (fitted_model$endpoint_types[2] == "continuous") {
        q_T0 = marginal_cont_constructor(marginal_Y0,
                                         theta0[(n_parX0 + 1):(n_parX0 + n_parY0)])$inv_cdf
        q_T1 = marginal_cont_constructor(marginal_Y1,
                                         theta1[(n_parX1 + 1):(n_parX1 + n_parY1)])$inv_cdf
      }
      else {
        q_T0 = marginal_ord_constructor(theta0[(n_parX0 + 1):(n_parX0 + n_parY0)])$inv_cdf
        q_T1 = marginal_ord_constructor(theta1[(n_parX0 + 1):(n_parX0 + n_parY0)])$inv_cdf
      }

      if (switch_XY) {
        # Switch S and T if necessary
        q_S0_old = q_S0
        q_S1_old = q_S1

        q_S0 = q_T0
        q_S1 = q_T1
        q_T0 = q_S0_old
        q_T1 = q_S1_old
      }

      if (all(fitted_model$endpoint_types == c("continuous", "continuous"))) {
        compute_ICA = compute_ICA_ContCont
      }
      else if (all(fitted_model$endpoint_types == c("ordinal", "continuous")) |
               all(fitted_model$endpoint_types == c("continuous", "ordinal"))) {
        compute_ICA = compute_ICA_OrdCont
      }
      else if (all(fitted_model$endpoint_types == c("ordinal", "ordinal"))) {
        compute_ICA = compute_ICA_OrdOrd
      }

      ICA = compute_ICA(
        copula_par = c(c_12,
                       copula_par_unid[1],
                       c_34,
                       copula_par_unid[2:4]),
        rotation_par = c(r_12,
                         rotation_par_unid[1],
                         r_34,
                         rotation_par_unid[2:4]),
        copula_family1 = copula_family1,
        copula_family2 = copula_family2,
        n_prec = n_prec,
        q_S0 = q_S0,
        q_T0 = q_T0,
        q_S1 = q_S1,
        q_T1 = q_T1,
        ICA_estimator = ICA_estimator,
        seed = seed,
        marginal_sp_rho = marginal_sp_rho
      )[measure]

      return(ICA)
    }
  }


  return(ICA_given_model)
}
