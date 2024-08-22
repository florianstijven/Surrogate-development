#' Perform Sensitivity Analysis for the Individual Causal Association based on a
#' D-vine copula model
#'
#' @inheritSection sensitivity_analysis_SurvSurv_copula Information-Theoretic Causal Inference Framework
#' @inheritSection estimate_ICA_OrdCont Individual Causal Association
#' @inheritSection sensitivity_analysis_SurvSurv_copula Sensitivity Analysis
#' @inheritSection sensitivity_analysis_SurvSurv_copula Additional Assumptions
#'
#' @param fitted_model Returned value from [fit_copula_model_OrdOrd()],
#'   [fit_copula_model_OrdCont()], or [fit_copula_model_ContCont()]. This object
#'   contains the estimated identifiable part of the joint distribution for the
#'   potential outcomes.
#' @param mutinfo_estimator Function that estimates the mutual information
#'  between the first two arguments which are numeric vectors. Defaults to
#'  `FNN::mutinfo()` for the continuous-continuous setting.
#' @inheritParams sensitivity_analysis_SurvSurv_copula
#' @inheritParams sample_copula_parameters
#'
#' @return A data frame is returned. Each row represents one replication in the
#'   sensitivity analysis. The returned data frame always contains the following
#'   columns:
#' * `R2H`, `sp_rho`: ICA as quantified by \eqn{R^2_H} and Spearman's rho, respectively.
#' * `c12`, `c34`: estimated copula parameters.
#' * `c23`, `c13_2`, `c24_3`, `c14_23`: sampled copula parameters of the
#'   unidentifiable copulas in the D-vine copula. The parameters correspond to
#'   the parameterization of the `copula_family2` copula as in the `copula`
#'   R-package.
#' * `r12`, `r34`: Fixed rotation parameters for the two identifiable copulas.
#' * `r23`, `r13_2`, `r24_3`, `r14_23`: Sampled rotation parameters of the
#'   unidentifiable copulas in the D-vine copula. These values are constant for
#'   the Gaussian copula family since that copula is invariant to rotations.
#'
#' The returned data frame also contains the following columns when
#' `marg_association` is `TRUE`:
#' * `sp_s0s1`, `sp_s0t0`, `sp_s0t1`, `sp_s1t0`, `sp_s1t1`, `sp_t0t1`:
#' Spearman's rho between the corresponding potential outcomes. Note that these
#' associations refer to the observable potential outcomes. In contrast, the
#' estimated association parameters from [fit_copula_model_OrdOrd()] and
#' [fit_copula_model_OrdCont] refer to associations on a latent scale.
#' @export
#'
#' @references
#' Alonso, A. (2018). An information-theoretic approach for the evaluation of
#' surrogate endpoints. In Wiley StatsRef: Statistics Reference Online. John
#' Wiley & Sons, Ltd.
#'
#' Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., and Burzykowski, T.
#' (2015). On the relationship between the causal-inference and meta-analytic
#' paradigms for the validation of surrogate endpoints. Biometrics 71, 15–24.
sensitivity_analysis_copula = function(fitted_model,
                                       n_sim,
                                       eq_cond_association = TRUE,
                                       lower = c(-1, -1, -1, -1),
                                       upper = c(1, 1, 1, 1),
                                       degrees = c(0, 90, 180, 270),
                                       marg_association = TRUE,
                                       copula_family2 = fitted_model$copula_family[1],
                                       n_prec = 1e4,
                                       ncores = 1,
                                       mutinfo_estimator = NULL
) {
  # Extract relevant estimated parameters/objects for the fitted copula model.

  # Identifiable copula families.
  copula_family1 = c(fitted_model$fit_0$copula_family,
                     fitted_model$fit_1$copula_family)
  # If copula_family2 contains only 1 element, this vector is appended to
  # the correct length.
  if(length(copula_family2) == 1) copula_family2 = rep(copula_family2, 4)

  # The "X" variable in fitted_model can correspond to the surrogate or the true
  # endpoint. This is checked manually.
  if(fitted_model$fit_0$names_XY[1] == "Surr") {
    q_S0 = fitted_model$fit_0$marginal_X$inv_cdf
    q_S1 = fitted_model$fit_1$marginal_X$inv_cdf
    q_T0 = fitted_model$fit_0$marginal_Y$inv_cdf
    q_T1 = fitted_model$fit_1$marginal_Y$inv_cdf
  }
  if(fitted_model$fit_0$names_XY[1] == "True") {
    q_S0 = fitted_model$fit_0$marginal_Y$inv_cdf
    q_S1 = fitted_model$fit_1$marginal_Y$inv_cdf
    q_T0 = fitted_model$fit_0$marginal_X$inv_cdf
    q_T1 = fitted_model$fit_1$marginal_X$inv_cdf
  }

  c12 = coef(fitted_model$fit_0$ml_fit)[length(coef(fitted_model$fit_0$ml_fit))]
  c34 = coef(fitted_model$fit_1$ml_fit)[length(coef(fitted_model$fit_1$ml_fit))]

  # For the Gaussian copula, fisher's Z transformation was applied. We have to
  # backtransform to the correlation scale in that case.
  if (copula_family1[1] == "gaussian") {
    c12 = (exp(2 * c12) - 1) / (exp(2 * c12) + 1)
  }
  if (copula_family1[2] == "gaussian") {
    c34 = (exp(2 * c34) - 1) / (exp(2 * c34) + 1)
  }
  # Sample unidentifiable copula parameters.
  c = sample_copula_parameters(
    copula_family2 = copula_family2,
    n_sim = n_sim,
    eq_cond_association = eq_cond_association,
    lower = lower,
    upper = upper
  )
  c = cbind(rep(c12, n_sim), c[, 1], rep(c34, n_sim), c[, 2:4])
  c_list = purrr::map(.x = split(c, seq(nrow(c))), .f = as.double)

  # Sample rotation parameters. If the copula family does allows for negative
  # associations, the rotation parameter is set to zero.
  r = sample_rotation_parameters(n_sim, degrees = degrees)
  if (copula_family2[1] %in% c("gaussian", "frank")) r[, 1] = rep(0, n_sim)
  if (copula_family2[2] %in% c("gaussian", "frank")) r[, 2] = rep(0, n_sim)
  if (copula_family2[3] %in% c("gaussian", "frank")) r[, 3] = rep(0, n_sim)
  if (copula_family2[4] %in% c("gaussian", "frank")) r[, 4] = rep(0, n_sim)

  r = cbind(rep(0, n_sim), r[, 1], rep(0, n_sim), r[, 2:4])
  r_list = purrr::map(.x = split(r, seq(nrow(r))), .f = as.double)
  # For every set of sampled unidentifiable parameters, compute the
  # required quantities.

  # Put all other arguments in a list for the apply function.
  MoreArgs = list(
    copula_family1 = copula_family1,
    copula_family2 = copula_family2,
    n_prec = n_prec,
    q_S0 = q_S0,
    q_S1 = q_S1,
    q_T0 = q_T0,
    q_T1 = q_T1,
    marginal_sp_rho = marg_association,
    seed = 1,
    mutinfo_estimator = mutinfo_estimator
  )

  # Determine which ICA to compute.
  if (all(fitted_model$endpoint_types == c("continuous", "continuous"))) {
    compute_ICA = compute_ICA_ContCont
  }
  if (all(fitted_model$endpoint_types == c("ordinal", "continuous"))) {
    compute_ICA = compute_ICA_OrdCont
  }
  if (all(fitted_model$endpoint_types == c("ordinal", "ordinal"))) {
    compute_ICA = compute_ICA_OrdOrd
  }

  if (ncores > 1 & requireNamespace("parallel")) {
    cl  <- parallel::makeCluster(ncores)
    print("Starting parallel simulations")
    # Get current search path and set the same search path in the new instances
    # of R. Usually, this would not be necessary, but if the user changed the
    # search path before running this function, there could be an issue.
    search_path = .libPaths()
    force(search_path)
    parallel::clusterExport(
      cl = cl,
      varlist = c("fitted_model", "search_path"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl = cl, expr = .libPaths(new = search_path))
    temp = parallel::clusterMap(
      cl = cl,
      fun = compute_ICA,
      copula_par = c_list,
      rotation_par = r_list,
      MoreArgs = MoreArgs
    )
    print("Finishing parallel simulations")
    on.exit(expr = {
      parallel::stopCluster(cl)
      rm("cl")
    })
  }
  else if (ncores == 1) {
    temp = mapply(
      FUN = compute_ICA,
      copula_par = c_list,
      rotation_par = r_list,
      MoreArgs = MoreArgs,
      SIMPLIFY = FALSE
    )
  }

  measures_df = t(as.data.frame(temp))
  rownames(measures_df) = NULL

  colnames(c) = c("c12", "c23", "c34", "c13_2", "c24_3", "c14_23")
  colnames(r) = c("r12", "r23", "r34", "r13_2", "r24_3", "r14_23")

  sens_results = cbind(as.data.frame(measures_df), c, r)

  attr(sens_results, which = "copula_family1") = copula_family1
  attr(sens_results, which = "copula_family2") = copula_family2
  attr(sens_results, which = "composite") = FALSE
  # Return a data frame with the results of the sensitivity analysis.
  return(sens_results)

}
