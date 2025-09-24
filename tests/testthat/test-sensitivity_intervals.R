test_that("Intervals of ignorance and uncertainty are correctly computed for survival endpoints", {
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
  suppressWarnings({
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
      mutinfo_estimator = function(x, y)
        - 0.5 *  log(1 - stats::cor(x, y, method = "spearman")),
      restr_time = 1
    )
  })

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
      0.903859303709,
      0.997214965810,
      0.890634385971,
      0.995473758827,
      0.892237833093,
      0.995439312912
    ),
    tolerance = 0.1
  )
  expect_equal(
    c(
      est_interval_of_ignorance_comp,
      interval_of_uncertainty_strong_comp,
      interval_of_uncertainty_pointwise_comp
    ),
    c(
      0.728857500144,
      0.990287672120,
      0.581820825288,
      0.991316728016,
      0.590704158826,
      0.991155587546
    ),
    tolerance = 0.1
  )
  expect_equal(
    c(
      est_interval_of_ignorance_comp_sprho_restr,
      interval_of_uncertainty_strong_comp_sprho_restr,
      interval_of_uncertainty_pointwise_comp_sprho_restr
    ),
    c(
      0.613278117755,
      0.967098312700,
      0.520161800751,
      0.973583298994,
      0.521745786609,
      0.973544737359
    ),
    tolerance = 0.1
  )
})

test_that("Intervals of ignorance and uncertainty are correctly computed for ordinal/continuous endpoints", {
  # Load fitted copula model.
  fitted_model_contcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-ContCont.rds"))
  fitted_model_ordcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-OrdCont.rds"))
  fitted_model_ordord = readRDS(test_path("fixtures", "schizo-dvine-gaussian-OrdOrd.rds"))
  # Load precomputed data set that is outputted by the sensitivity analysis
  # functions.
  sens_results_contcont = readRDS(test_path("fixtures", "sens-results-schizo-clayton-ContCont.rds"))
  sens_results_ordcont = readRDS(test_path("fixtures", "sens-results-schizo-clayton-OrdCont.rds"))
  sens_results_ordord = readRDS(test_path("fixtures", "sens-results-schizo-gaussian-OrdOrd.rds"))
  # Create sensitivity intervals object
  set.seed(1)
  sensitivity_intervals_estimated_contcont = sensitivity_intervals_Dvine(
    fitted_model = fitted_model_contcont,
    sens_results = sens_results_contcont,
    B = 20,
    n_prec = 1e3
  )
  set.seed(1)
  sensitivity_intervals_estimated_ordcont = sensitivity_intervals_Dvine(
    fitted_model = fitted_model_ordcont,
    sens_results = sens_results_ordcont,
    B = 20,
    n_prec = 1e3
  )
  set.seed(1)
  sensitivity_intervals_estimated_ordord = sensitivity_intervals_Dvine(
    fitted_model = fitted_model_ordord,
    sens_results = sens_results_ordord,
    B = 20,
    n_prec = 1e3
  )
  est_interval_of_ignorance_contcont = sensitivity_intervals_estimated_contcont$est_interval_of_ignorance
  interval_of_uncertainty_strong_contcont = sensitivity_intervals_estimated_contcont$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise_contcont = sensitivity_intervals_estimated_contcont$interval_of_uncertainty_pointwise_coverage

  est_interval_of_ignorance_ordcont = sensitivity_intervals_estimated_ordcont$est_interval_of_ignorance
  interval_of_uncertainty_strong_ordcont = sensitivity_intervals_estimated_ordcont$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise_ordcont = sensitivity_intervals_estimated_ordcont$interval_of_uncertainty_pointwise_coverage

  est_interval_of_ignorance_ordord = sensitivity_intervals_estimated_ordord$est_interval_of_ignorance
  interval_of_uncertainty_strong_ordord = sensitivity_intervals_estimated_ordord$interval_of_uncertainty_strong_coverage
  interval_of_uncertainty_pointwise_ordord = sensitivity_intervals_estimated_ordord$interval_of_uncertainty_pointwise_coverage

  expect_equal(
    c(
      est_interval_of_ignorance_contcont,
      interval_of_uncertainty_strong_contcont,
      interval_of_uncertainty_pointwise_contcont
    ),
    c(
      0.764378683162,
      0.937104257947,
      0.751912123271,
      0.951424114168,
      0.752176460121,
      0.951197679875
    ),
    tolerance = 0.1
  )
  expect_equal(
    c(
      est_interval_of_ignorance_ordcont,
      interval_of_uncertainty_strong_ordcont,
      interval_of_uncertainty_pointwise_ordcont
    ),
    c(
      0.218531467403,
      0.370530021884,
      0.177987805110,
      0.378222005120,
      0.182317520585,
      0.375366736593
    ),
    tolerance = 0.1
  )
  expect_equal(
    c(
      est_interval_of_ignorance_ordord,
      interval_of_uncertainty_strong_ordord,
      interval_of_uncertainty_pointwise_ordord
    ),
    c(
      0.0837791500542,
      0.1932191584169,
      0.0708307634417,
      0.2094639344559,
      0.0725304694619,
      0.2072494890744
    ),
    tolerance = 0.1
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
  suppressWarnings({
    sensitivity_intervals_estimated = sensitivity_intervals_Dvine(
      fitted_model = fitted_model,
      sens_results = sens_results,
      B = B,
      n_prec = n_prec,
      ncores = 2
    )
  })

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
      0.903859303709,
      0.997214965810,
      0.890634385971,
      0.995473758827,
      0.892237833093,
      0.995439312912
    ),
    tolerance = 0.1
  )
})
