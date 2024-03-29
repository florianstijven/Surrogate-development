test_that("Intervals of ignorance and uncertainty are correctly computed", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  fitted_model_comp = readRDS(test_path("fixtures", "ovarian-dvine-gaussian-scr.rds"))
  # Load precomputed data set that is outputted by the sensitivity analysis
  # functions.
  sens_results = readRDS(test_path("fixtures", "sens-results-ovarian-clayton.rds"))
  sens_results_comp = readRDS(test_path("fixtures", "sens-results-ovarian-gaussian-comp.rds"))
  sens_results_comp_sprho_restr = readRDS(test_path("fixtures", "sens-results-ovarian-gaussian-comp-sprho-restr.rds"))
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
  sensitivity_intervals_estimated_comp_sprho_restr = sensitivity_intervals_Dvine(
    fitted_model = fitted_model_comp,
    sens_results = sens_results_comp_sprho_restr,
    B = 20,
    n_prec = 1e3,
    mutinfo_estimator = function(x, y) 1 - exp(-2 * stats::cor(x, y, method = "spearman")),
    restr_time = 1
  )
  est_interval_of_ignorance = sensitivity_intervals_estimated$est_interval_of_ignorance
  interval_of_uncertainty_strong = sensitivity_intervals_estimated$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise = sensitivity_intervals_estimated$interval_of_uncertainty_pointwise_coverage

  est_interval_of_ignorance_comp = sensitivity_intervals_estimated_comp$est_interval_of_ignorance
  interval_of_uncertainty_strong_comp = sensitivity_intervals_estimated_comp$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise_comp = sensitivity_intervals_estimated_comp$interval_of_uncertainty_pointwise_coverage

  est_interval_of_ignorance_comp_sprho_restr = sensitivity_intervals_estimated_comp_sprho_restr$est_interval_of_ignorance
  interval_of_uncertainty_strong_comp_sprho_restr = sensitivity_intervals_estimated_comp_sprho_restr$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise_comp_sprho_restr = sensitivity_intervals_estimated_comp_sprho_restr$interval_of_uncertainty_pointwise_coverage


  expect_equal(
    c(
      est_interval_of_ignorance,
      interval_of_uncertainty_strong,
      interval_of_uncertainty_pointwise
    ),
    c(
      0.835910421200,
      0.997707840148,
      0.748958557291,
      0.997291085878,
      0.754652670814,
      0.997185131382
    )
  )
  expect_equal(
    c(
      est_interval_of_ignorance_comp,
      interval_of_uncertainty_strong_comp,
      interval_of_uncertainty_pointwise_comp
    ),
    c(
      0.642646510145,
      0.991589287727,
      0.504365537108,
      0.991103510500,
      0.521834842359,
      0.990989205985
    )
  )
  expect_equal(
    c(
      est_interval_of_ignorance_comp_sprho_restr,
      interval_of_uncertainty_strong_comp_sprho_restr,
      interval_of_uncertainty_pointwise_comp_sprho_restr
    ),
    c(
      0.751294689091,
      0.820984966266,
      0.711347430929,
      0.820793874360,
      0.712601776077,
      0.820715994919
    )
  )
})


test_that("Intervals of ignorance and uncertainty work on multiple cores", {
  testthat::skip_on_cran()
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Load precomputed data set that is outputted by the sensitivity analysis
  # functions.
  sens_results = readRDS(test_path("fixtures", "sens-results-ovarian-clayton.rds"))
  # Create sensitivity intervals object
  set.seed(1)
  n_prec = 1e3; B = 20
  sensitivity_intervals_estimated = sensitivity_intervals_Dvine(
    fitted_model = fitted_model,
    sens_results = sens_results,
    B = B,
    n_prec = n_prec,
    ncores = 2
  )
  est_interval_of_ignorance = sensitivity_intervals_estimated$est_interval_of_ignorance
  interval_of_uncertainty_strong = sensitivity_intervals_estimated$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise = sensitivity_intervals_estimated$interval_of_uncertainty_pointwise_coverage

  expect_equal(
    c(
      est_interval_of_ignorance,
      interval_of_uncertainty_strong,
      interval_of_uncertainty_pointwise
    ),
    c(
      0.835910421200,
      0.997707840148,
      0.748958557291,
      0.997291085878,
      0.754652670814,
      0.997185131382
    )
  )
})
