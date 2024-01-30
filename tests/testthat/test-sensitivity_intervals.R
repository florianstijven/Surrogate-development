test_that("Intervals of ignorance and uncertainty are correctly computed", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  fitted_model_comp = readRDS(test_path("fixtures", "ovarian-dvine-gaussian-scr.rds"))
  # Load precomputed data set that is outputted by the sensitivity analysis
  # functions.
  sens_results = readRDS(test_path("fixtures", "sens-results-ovarian-clayton.rds"))
  sens_results_comp = readRDS(test_path("fixtures", "sens-results-ovarian-gaussian-comp.rds"))
  # Create sensitivity intervals object
  set.seed(1)
  sensitivity_intervals_estimated = sensitivity_intervals_Dvine(
    fitted_model = fitted_model,
    sens_results = sens_results,
    B = 20,
    n_prec = 1e3
  )
  sensitivity_intervals_estimated_comp = sensitivity_intervals_Dvine(
    fitted_model = fitted_model_comp,
    sens_results = sens_results_comp,
    B = 20,
    n_prec = 1e3
  )
  est_interval_of_ignorance = sensitivity_intervals_estimated$est_interval_of_ignorance
  interval_of_uncertainty_strong = sensitivity_intervals_estimated$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise = sensitivity_intervals_estimated$interval_of_uncertainty_pointwise_coverage

  est_interval_of_ignorance_comp = sensitivity_intervals_estimated_comp$est_interval_of_ignorance
  interval_of_uncertainty_strong_comp = sensitivity_intervals_estimated_comp$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise_comp = sensitivity_intervals_estimated_comp$interval_of_uncertainty_pointwise_coverage

  expect_equal(
    c(est_interval_of_ignorance, interval_of_uncertainty_strong, interval_of_uncertainty_pointwise),
    c(0.835910421200, 0.997707840148, 0.748958557291, 0.997230581629, 0.754652670814, 0.997228210639)
  )
  expect_equal(
    c(est_interval_of_ignorance_comp, interval_of_uncertainty_strong_comp, interval_of_uncertainty_pointwise_comp),
    c(0.730641382393, 0.988266353657, 0.680891211087, 0.988232184888, 0.695156138941, 0.98815982158)
  )
})


