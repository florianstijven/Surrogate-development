test_that("Intervals of ignorance and uncertainty are correctly computed", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Load precomputed data set that is outputted by the sensitivity analysis
  # functions.
  sens_results = readRDS(test_path("fixtures", "sens-results-ovarian-clayton.rds"))
  # Create sensitivity intervals object
  set.seed(1)
  sensitivity_intervals_estimated = sensitivity_intervals_Dvine(
    fitted_model = fitted_model,
    sens_results = sens_results,
    B = 20,
    n_prec = 1e3
  )
  est_interval_of_ignorance = sensitivity_intervals_estimated$est_interval_of_ignorance
  interval_of_uncertainty_strong = sensitivity_intervals_estimated$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise = sensitivity_intervals_estimated$interval_of_uncertainty_pointwise_coverage

  expect_equal(
    c(est_interval_of_ignorance, interval_of_uncertainty_strong, interval_of_uncertainty_pointwise),
    c(0.930258399130, 0.998591189134, 0.848113716973, 0.997954708498, 0.851955788189, 0.997932000240)
  )
})


