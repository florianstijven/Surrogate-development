#' Compute Sensitivity Intervals
#'
#' [sensitivity_intervals_Dvine()] computes the estimated intervals of ignorance
#' and uncertainty within the information-theoretic causal inference framework
#' when the data are modeled with a D-vine copula model.
#'
#' @details
#' # Intervals of Ignorance and Uncertainty
#'
#' Vansteelandt et al. (2006) formalized sensitivity analysis for partly
#' identifiable parameters in the context of missing data and MNAR. These
#' concepts can be applied to the estimation of the ICA. Indeed, the ICA is also
#' partly identifiable because 50% if the potential outcomes are missing.
#'
#' Vansteelandt et al. (2006) replace a point estimate with a interval estimate:
#' the estimated interval of ignorance. In addition, they proposed several
#' extension of the classic confidence interval together with appropriate
#' definitions of coverage; these are termed intervals of uncertainty.
#'
#' [sensitivity_intervals_Dvine()] implements the estimated interval of
#' ignorance and the pointwise and strong intervals of uncertainty. Let \eqn{\boldsymbol{\nu}_l}
#' and \eqn{\boldsymbol{\nu}_u} be the values for the sensitivity parameter that
#' lead to the lowest and largest ICA, respectively, while fixing the identifiable
#' parameter at its estimated value \eqn{\hat{\boldsymbol{\beta}}}. See also
#' [summary_level_bootstrap_ICA()]. The following intervals are implemented:
#'
#' 1. *Estimated interval of ignorance*. This interval is defined as
#' \eqn{[ICA(\hat{\boldsymbol{\beta}}, \boldsymbol{\nu}_l), ICA(\hat{\boldsymbol{\beta}}, \boldsymbol{\nu}_u)]}.
#' 2. *Pointiwse interval of uncertainty*. Let \eqn{C_l} (and \eqn{C_u}) be the
#' lower (and upper) limit of a one-sided \eqn{1 - \alpha} CI for
#' \eqn{ICA(\boldsymbol{\beta_0}, \boldsymbol{\nu}_l)} (and
#' \eqn{ICA(\boldsymbol{\beta_0}, \boldsymbol{\nu}_l)}). This interval is then
#' defined as \eqn{[C_l, C_u]} when the ignorance is much larger than the
#' statistical imprecision.
#' 3. *Strong interval of uncertainty*. Let \eqn{C_l} (and \eqn{C_u}) be the
#' lower (and upper) limit of a two-sided \eqn{1 - \alpha} CI for
#' \eqn{ICA(\boldsymbol{\beta_0}, \boldsymbol{\nu}_l)} (and
#' \eqn{ICA(\boldsymbol{\beta_0}, \boldsymbol{\nu}_l)}). This interval is then
#' defined as \eqn{[C_l, C_u]}.
#'
#' The CIs, which are need for the intervals of uncertainty, are based on
#' percentile bootstrap confidence intervals, as documented in
#' [summary_level_bootstrap_ICA()]. In addition, \eqn{\boldsymbol{\nu}_l} is not
#' known. Therefore, it is estimated as
#' \deqn{\arg \min_{\boldsymbol{\nu} \in \Gamma} ICA(\hat{\boldsymbol{\beta}}, \boldsymbol{\nu}),}
#' and similarly for \eqn{\boldsymbol{\nu}_u}.
#'
#' @param sens_results Dataframe returned by
#'   [sensitivity_analysis_SurvSurv_copula()]. If additional assumptions need to
#'   be incorporated, this dataframe can first be filtered.
#' @param measure Compute intervals for which measure of surrogacy? Defaults to
#'   `"ICA"`. See first column names of `sens_results` for other possibilities.
#' @inheritParams Dvine_ICA_confint
#' @inheritParams summary_level_bootstrap_ICA
#' @inheritParams sensitivity_analysis_SurvSurv_copula
#'
#' @return An S3 object of the class `sensitivity_intervals_Dvine` which can be
#' printed.
#' @export
#'
#' @inherit sensitivity_analysis_SurvSurv_copula examples
#' @references
#' Vansteelandt, Stijn, et al. "Ignorance and uncertainty regions as inferential
#' tools in a sensitivity analysis." Statistica Sinica (2006): 953-979.
sensitivity_intervals_Dvine = function(fitted_model,
                                       sens_results,
                                       measure = "ICA",
                                       B = 2e2,
                                       alpha = 0.05,
                                       n_prec = 5e3,
                                       mutinfo_estimator = NULL,
                                       restr_time = +Inf) {
  # Select row with the smallest value for the ICA.
  upper_limit_row = sens_results[which.max(sens_results[[measure]]), ]
  # Select row with the largest value for the ICA.
  lower_limit_row = sens_results[which.min(sens_results[[measure]]), ]

  copula_family2 = attr(sens_results, "copula_family2")
  composite = attr(sens_results, "composite")

  est_interval_of_ignorance = c(
    as.numeric(lower_limit_row[measure]),
    as.numeric(upper_limit_row[measure])
  )

  # Perform parametric bootstrap for setting that yielded smallest ICA.
  bootstrap_replications_lower = summary_level_bootstrap_ICA(
    fitted_model = fitted_model,
    copula_par_unid = as.numeric(lower_limit_row[1, c("c23", "c13_2", "c24_3", "c14_23")]),
    copula_family2 = copula_family2,
    rotation_par_unid = as.numeric(lower_limit_row[1, c("r23", "r13_2", "r24_3", "r14_23")]),
    n_prec = n_prec,
    composite = composite,
    B = B,
    seed = 1,
    measure = measure,
    mutinfo_estimator = mutinfo_estimator,
    restr_time = restr_time
  )

  # Perform parametric bootstrap for setting that yielded largest ICA.
  bootstrap_replications_upper = summary_level_bootstrap_ICA(
    fitted_model = fitted_model,
    copula_par_unid = as.numeric(upper_limit_row[1, c("c23", "c13_2", "c24_3", "c14_23")]),
    copula_family2 = copula_family2,
    rotation_par_unid = as.numeric(upper_limit_row[1, c("r23", "r13_2", "r24_3", "r14_23")]),
    n_prec = n_prec,
    composite = composite,
    B = B,
    seed = 1,
    measure = measure,
    mutinfo_estimator = mutinfo_estimator,
    restr_time = restr_time
  )

  # Compute confidence intervals for the two scenarios selected above. For the
  # uncertainty intervals with pointwise 1 - alpha coverage, we can consider the
  # lower (upper) limit of a one-sided 1 - alpha confidence interval for the smallest
  # (largest) estimated ICA.
  interval_of_uncertainty_pointwise_coverage = c(
    stats::quantile(bootstrap_replications_lower, prob = alpha),
    stats::quantile(bootstrap_replications_upper, prob = 1 - alpha)
  )
  # For the uncertainty interval with strong 1 - alpha coverage, we need to
  # consider the lower (upper) limit of two-sided 1 - alpha confidence intervals
  # for the smallest (largest) estimated ICA.
  interval_of_uncertainty_strong_coverage = c(
    stats::quantile(bootstrap_replications_lower, prob = alpha / 2),
    stats::quantile(bootstrap_replications_upper, prob = 1 - alpha / 2)
  )

  return(
    new_sensitivity_intervals_Dvine(
      fitted_model,
      upper_limit = list(
        scenario_information = upper_limit_row
        # confint_alpha = upper_confint_alpha,
        # confint_2alpha = upper_confint_2alpha
      ),
      lower_limit = list(
        scenario_information = lower_limit_row
        # confint_alpha = lower_confint_alpha,
        # confint_2alpha = lower_confint_2alpha
      ),
      est_interval_of_ignorance = est_interval_of_ignorance,
      interval_of_uncertainty_pointwise_coverage = interval_of_uncertainty_pointwise_coverage,
      interval_of_uncertainty_strong_coverage = interval_of_uncertainty_strong_coverage,
      alpha = alpha,
      copula_family2 = copula_family2
    )
  )
}

new_sensitivity_intervals_Dvine = function(fitted_model,
                                           upper_limit,
                                           lower_limit,
                                           est_interval_of_ignorance,
                                           interval_of_uncertainty_pointwise_coverage,
                                           interval_of_uncertainty_strong_coverage,
                                           alpha,
                                           copula_family2) {
  names(est_interval_of_ignorance) = NULL
  names(interval_of_uncertainty_pointwise_coverage) = NULL
  names(interval_of_uncertainty_strong_coverage) = NULL
  structure(
    .Data = list(
      fitted_model = fitted_model,
      upper_limit = upper_limit,
      lower_limit = lower_limit,
      est_interval_of_ignorance = est_interval_of_ignorance,
      interval_of_uncertainty_pointwise_coverage = interval_of_uncertainty_pointwise_coverage,
      interval_of_uncertainty_strong_coverage = interval_of_uncertainty_strong_coverage,
      alpha = alpha,
      copula_family2 = copula_family2
    ),
    class = "sensitivity_intervals_Dvine"
  )
}

print.sensitivity_intervals_Dvine = function(x, ...) {
  cat("Sensitivity analysis based on the D-vine copula model\n\n")
  cat("Identifiable Copula Family (c12 and c34): "); cat(x$fitted_model$copula_family); cat("\n")
  cat("Unidentifiable Copula Family (c23, c13_2, c24_3, and c14_23): "); cat(x$copula_family2); cat("\n")
  cat("Number of Internal Knots: "); cat(length(x$fitted_model$knots0) - 2)
  cat("\n\n")
  cat("Estimated interval of ignorance: "); print_interval(x$est_interval_of_ignorance); cat("\n")
  cat("Interval of uncertainty (pointwise): "); print_interval(x$interval_of_uncertainty_pointwise_coverage); cat("\n")
  cat("Interval of uncertainty (strong): "); print_interval(x$interval_of_uncertainty_strong_coverage); cat("\n")
  cat("alpha = "); cat(x$alpha)
  cat("\n")
}

print_interval = function(x, digits = 3) {
  cat("(")
  cat(round(x[1], digits))
  cat(", ")
  cat(round(x[2], digits))
  cat(")")
}
