test_that("Intervals of ignorance and uncertainty are correctly computed", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Load precomputed data set that is outputted by the sensitivity analysis
  # functions.
  sens_results = readRDS(test_path("fixtures", "sens-results-ovarian-clayton.rds"))
  # Create sensitivity intervals object
  sensitivity_intervals_estimated = sensitivity_intervals_Dvine(
    fitted_model = fitted_model,
    sens_results = sens_results,
    B = 5,
    n_prec = 1e3
  )
  interval_of_ignorance = sensitivity_intervals_estimated$interval_of_ignorance
  interval_of_uncertainty = sensitivity_intervals_estimated$interval_of_uncertainty

  expect_equal(
    c(interval_of_ignorance, interval_of_uncertainty),
    c(0.930258399130, 0.998591189134, 0.872962115841, 0.997528944411)
  )
})
