test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Clayton copula", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    lower = c(-1, 0, 0, -1),
    upper = c(1, 0, 0, 1),
    n_sim = 5,
    n_prec = 500
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1])
  check_vector = c(0.988990899030, 0.993308293233)
  expect_equal(output_vector, check_vector, tolerance = 1)
  # We're only checking whether the functions run without errors;
  # hence the large tolerance.
})

test_that("sensitivity_analysis_SurvSurv_copula() works with variable number of knots", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-variable.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    lower = c(-1, 0, 0, -1),
    upper = c(1, 0, 0, 1),
    n_sim = 5,
    n_prec = 500
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1])
  check_vector = c(0.979074242038, 0.971166764667)
  expect_equal(output_vector, check_vector, tolerance = 1)
  # We're only checking whether the functions run without errors;
  # hence the large tolerance.
})

test_that("sensitivity_analysis_SurvSurv_copula() works on 2 cores with Clayton copula", {
  # Do not run this test on CRAN because multiple cores may not be available on
  # the CRAN computer.
  skip_on_cran()
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    lower = c(-1, 0, 0, -1),
    upper = c(1, 0, 0, 1),
    n_sim = 5,
    n_prec = 500,
    ncores = 2
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1])
  check_vector = c(0.988990899030, 0.993308293233)
  expect_equal(output_vector, check_vector, tolerance = 1)
  # We're only checking whether the functions run without errors;
  # hence the large tolerance.
})

test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with variable copula", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-variable.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    n_sim = 5,
    n_prec = 500
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1])
  check_vector = c(0.973975831732, 0.972747026988)
  expect_equal(output_vector, check_vector, tolerance = 1)
  # We're only checking whether the functions run without errors;
  # hence the large tolerance.
})

test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Clayton copula and four different unidentifiable copulas", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    n_sim = 1,
    n_prec = 2e3,
    copula_family2 = c("clayton", "frank", "gaussian", "frank")
  )
  set.seed(1)
  sens_results_cond1 = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    lower = c(-1, 0, 0, -1),
    upper = c(1, 0, 0, 1),
    n_sim = 1,
    n_prec = 2e3,
    copula_family2 = c("clayton", "frank", "gaussian", "frank")
  )
  set.seed(1)
  sens_results_cond2 = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    lower = c(-1, 0, 0, -1),
    upper = c(1, 0, 0, 1),
    n_sim = 1,
    n_prec = 2e3,
    copula_family2 = c("clayton", "frank", "frank", "frank")
  )

  # Check results for setting without conditional independence.
  expect_equal(
    c(sens_results$ICA[1],
      sens_results$sp_rho[1]),
    c(0.996814266963, 0.995511262378),
    tolerance = 1
  )
  # We're only checking whether the functions run without errors;
  # hence the large tolerance.

  # Check that results for two conditional independence settings are identical.
  # The only difference is in the copula for which we assume conditional
  # independence.
  expect_equal(
    sens_results_cond1,
    sens_results_cond2, ignore_attr = "copula_family2"
  )
})

test_that("sensitivity_analysis_SurvSurv_copula() wth restricted survival times and Spearman's rho as ICA", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    lower = rep(0.5, 4),
    composite = TRUE,
    n_sim = 1,
    n_prec = 2e3,
    mutinfo_estimator = function(x, y) 1 - exp(-2 * stats::cor(x, y, method = "spearman")),
    restr_time = 2
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1])
  check_vector = c(0.822136179803, 0.995515244445)
  expect_equal(output_vector, check_vector, tolerance = 1)
  # We're only checking whether the functions run without errors;
  # hence the large tolerance.
})
