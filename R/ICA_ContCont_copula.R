#' Compute Individual Causal Association for a given D-vine copula model in the
#' Continuous-Continuous Setting
#'
#' The [compute_ICA_ContCont()] function computes the individual causal
#' association (and associated quantities) for a fully identified D-vine copula
#' model in the continuous-continuous setting.
#'
#' @inheritParams compute_ICA_SurvSurv
#'
#' @return (numeric) A Named vector with the following elements:
#'  * ICA
#'  * Spearman's rho, \eqn{\rho_s (\Delta S, \Delta T)} (if asked)
#'  * Marginal association parameters in terms of Spearman's rho (if asked):
#'  \deqn{\rho_{s}(T_0, S_0), \rho_{s}(T_0, S_1), \rho_{s}(T_0, T_1),
#'  \rho_{s}(S_0, S_1), \rho_{s}(S_0, T_1),
#'  \rho_{s}(S_1, T_1)}
compute_ICA_ContCont = function(copula_par,
                                rotation_par,
                                copula_family1,
                                copula_family2,
                                n_prec,
                                q_S0,
                                q_T0,
                                q_S1,
                                q_T1,
                                marginal_sp_rho = TRUE,
                                seed = 1,
                                mutinfo_estimator = NULL,
                                plot_deltas = FALSE) {
  computed_ICA = compute_ICA_SurvSurv(
    copula_par = copula_par,
    rotation_par = rotation_par,
    copula_family1 = copula_family1,
    copula_family2 = copula_family2,
    n_prec = n_prec,
    q_S0 = q_S0,
    q_T0 = q_T0,
    q_S1 = q_S1,
    q_T1 = q_T1,
    marginal_sp_rho = marginal_sp_rho,
    seed = seed,
    composite = FALSE,
    mutinfo_estimator = mutinfo_estimator,
    plot_deltas = plot_deltas,
    restr_time = +Inf
  )
  # Remove the computed population strata for survival which don't make sense in
  # the continuous-continuous setting.
  return(computed_ICA[1:8])
}
