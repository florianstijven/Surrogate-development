#' Compute Individual Causal Association for a given D-vine copula model in the
#' Survival-Survival Setting
#'
#' The [compute_ICA_SurvSurv()] function computes the individual causal
#' association (and associated quantities) for a fully identified D-vine copula
#' model in the survival-survival setting.
#'
#' @param q_T0
#' @param q_T1
#' @param seed
#'
#' @inheritParams compute_ICA_BinCont
#' @inheritParams estimate_mutual_information_SurvSurv
#'
#' @return
#' @export
#'
#' @examples
compute_ICA_SurvSurv = function(copula_par,
                                rotation_par,
                                copula_family1,
                                copula_family2 = copula_family1,
                                n_prec,
                                minfo_prec,
                                q_S0,
                                q_T0,
                                q_S1,
                                q_T1,
                                composite,
                                marginal_sp_rho = TRUE,
                                seed = 1) {
  withr::local_seed(seed)
  # Sample individual causal treatment effects from the given model. If
  # marginal_sp_rho = TRUE, then the Spearman's correlation matrix is also
  # computed for the 4d joint distribution of potential outcomes.
  delta_list = sample_deltas_BinCont(
    copula_par = copula_par,
    rotation_par = rotation_par,
    copula_family1 = copula_family1,
    copula_family2 = copula_family1,
    n = n_prec,
    q_S0 = q_S0,
    q_S1 = q_S1,
    q_T0 = q_T0,
    q_T1 = q_T1,
    marginal_sp_rho = marginal_sp_rho,
    setting = "SurvSurv",
    composite = composite
  )

  delta_df = delta_list$Delta_dataframe
  sp_rho_matrix = delta_list$marginal_sp_rho_matrix
  # Compute mutual information between Delta S and Delta T.
  mutual_information = estimate_mutual_information_SurvSurv(delta_df$DeltaS, delta_df$DeltaT, minfo_prec)
  # Compute ICA
  ICA = 1 - exp(-2 * mutual_information)
  sp_rho = stats::cor(delta_df$DeltaS, delta_df$DeltaT, method = "spearman")

  return(
    c(
      ICA = ICA,
      sp_rho = sp_rho,
      sp_rho_t0s0 = sp_rho_matrix[1, 2],
      sp_rho_t0s1 = sp_rho_matrix[1, 3],
      sp_rho_t0t1 = sp_rho_matrix[1, 4],
      sp_rho_s0s1 = sp_rho_matrix[2, 3],
      sp_rho_s0t1 = sp_rho_matrix[2, 3],
      sp_rho_s1t1 = sp_rho_matrix[3, 4],
      prop_harmed = delta_list$survival_classification["prop_harmed"],
      prop_protected = delta_list$survival_classification["prop_protected"],
      prop_always = delta_list$survival_classification["prop_always"],
      prop_never = delta_list$survival_classification["prop_never"]
    )
  )
}

#' Estimate the Mutual Information in the Survival-Survival Setting
#'
#' [estimate_mutual_information_SurvSurv()] estimates the mutual information for
#' a sample of individual causal treatment effects with a time-to-event
#' surrogate and a time-to-event true endpoint. The mutual information is
#' estimated by first estimating the bivariate density and then computing the
#' mutual information for the estimated density.
#'
#' @inheritParams estimate_ICA_BinCont
#' @param minfo_prec Number of quasi Monte-Carlo samples for the numerical
#'   integration to obtain the mutual information. If this value is 0 (default),
#'   the mutual information is not computed and `NA` is returned for the mutual
#'   information and derived quantities.
#'
#' @return (numeric) estimated mutual information.
#' @importFrom kdecopula kdecop dep_measures
estimate_mutual_information_SurvSurv = function(delta_S, delta_T, minfo_prec) {
  # When minfo = 0, we do not want to compute the mutual information. Therefore,
  # a NA is returned in that case.
  if (minfo_prec == 0) {
    return(NA)
  }
  # If minfo > 0, we go proceed with computing the mutual information.
  a = rank(delta_S) / (length(delta_S) + 1)
  b = rank(delta_T) / (length(delta_T) + 1)
  temp_matrix = matrix(data = c(a, b), ncol = 2)
  tryCatch(expr = {
    t_kde = kdecopula::kdecop(udata = temp_matrix, method = "TLL2nn")
    minfo = kdecopula::dep_measures(
      object = t_kde,
      n_qmc = minfo_prec + 1,
      measures = "minfo",
      seed = 1
    )
  })

  return(minfo)
}
