sensitivity_intervals_Dvine = function(fitted_model,
                                       sens_results,
                                       measure = "ICA",
                                       B = 2e2,
                                       alpha = 0.05,
                                       n_prec = 5e3) {
  # Select row with the smallest value for the ICA.
  upper_limit_row = sens_results[which.max(sens_results[[measure]]), ]
  # Select row with the largest value for the ICA.
  lower_limit_row = sens_results[which.min(sens_results[[measure]]), ]

  copula_family2 = attr(sens_results, "copula_family2")
  composite = attr(sens_results, "composite")

  # Compute confidence interval for the two scenarios selected above.
  upper_confint = Dvine_ICA_confint(
    fitted_model = fitted_model,
    alpha = alpha,
    copula_par_unid = as.numeric(upper_limit_row[1, c("c23", "c13_2", "c24_3", "c14_23")]),
    copula_family2 = copula_family2,
    rotation_par_unid = as.numeric(upper_limit_row[1, c("r23", "r13_2", "r24_3", "r14_23")]),
    n_prec = n_prec,
    composite = composite,
    B = B,
    seed = 1
  )
  lower_confint = Dvine_ICA_confint(
    fitted_model = fitted_model,
    alpha = alpha,
    copula_par_unid = as.numeric(lower_limit_row[1, c("c23", "c13_2", "c24_3", "c14_23")]),
    copula_family2 = copula_family2,
    rotation_par_unid = as.numeric(lower_limit_row[1, c("r23", "r13_2", "r24_3", "r14_23")]),
    n_prec = n_prec,
    composite = composite,
    B = B,
    seed = 1
  )

  interval_of_ignorance = c(
    as.numeric(lower_limit_row[measure]),
    as.numeric(upper_limit_row[measure])
  )
  interval_of_uncertainty = c(
    lower_confint[1],
    upper_confint[2]
  )

  return(
    new_sensitivity_intervals_Dvine(
      fitted_model,
      upper_limit = list(
        scenario_information = upper_limit_row,
        confint = upper_confint
      ),
      lower_limit = list(
        scenario_information = lower_limit_row,
        confint = lower_confint
      ),
      interval_of_ignorance = interval_of_ignorance,
      interval_of_uncertainty = interval_of_uncertainty,
      alpha = alpha,
      copula_family2 = copula_family2
    )
  )
}

new_sensitivity_intervals_Dvine = function(fitted_model,
                                           upper_limit,
                                           interval_of_ignorance,
                                           interval_of_uncertainty,
                                           lower_limit,
                                           alpha,
                                           copula_family2) {
  names(interval_of_ignorance) = NULL
  names(interval_of_uncertainty) = NULL
  structure(
    .Data = list(
      fitted_model = fitted_model,
      upper_limit = upper_limit,
      lower_limit = lower_limit,
      interval_of_ignorance = interval_of_ignorance,
      interval_of_uncertainty = interval_of_uncertainty,
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
  cat("Estimated interval of ignorance: "); print_interval(x$interval_of_ignorance); cat("\n")
  cat("Estimated interval of uncertainty: "); print_interval(x$interval_of_uncertainty); cat("\n")
  cat("alpha = "); cat(x$alpha)
}

print_interval = function(x, digits = 3) {
  cat("(")
  cat(round(x[1], digits))
  cat(", ")
  cat(round(x[2], digits))
  cat(")")
}
