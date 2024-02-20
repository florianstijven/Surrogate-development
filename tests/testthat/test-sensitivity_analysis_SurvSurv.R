test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Clayton copula", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 5,
    n_prec = 500
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1],
                    sens_results$c23[3])
  check_vector = c(0.983848446048, 0.984428257713, 1.374917942896)
  expect_equal(output_vector, check_vector)
})

test_that("sensitivity_analysis_SurvSurv_copula() works with variable number of knots", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-variable.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 5,
    n_prec = 500
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1],
                    sens_results$c23[3])
  check_vector = c(0.973170806302, 0.971849999400, 1.374917942896)
  expect_equal(output_vector, check_vector)
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
    cond_ind = TRUE,
    n_sim = 5,
    n_prec = 500,
    ncores = 2
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1],
                    sens_results$c23[3])
  check_vector = c(0.983848446048, 0.984428257713, 1.374917942896)
  expect_equal(output_vector, check_vector)
})

test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with variable copula", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-variable.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 5,
    n_prec = 500
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1],
                    sens_results$c23[3])
  check_vector = c(0.973170806302, 0.971849999400, 1.374917942896)
  expect_equal(output_vector, check_vector)
})

test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Clayton copula and four different unidentifiable copulas", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    cond_ind = FALSE,
    n_sim = 1,
    n_prec = 2e3,
    copula_family2 = c("clayton", "frank", "gaussian", "frank")
  )
  set.seed(1)
  sens_results_cond1 = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 1,
    n_prec = 2e3,
    copula_family2 = c("clayton", "frank", "gaussian", "frank")
  )
  set.seed(1)
  sens_results_cond2 = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 1,
    n_prec = 2e3,
    copula_family2 = c("clayton", "frank", "frank", "frank")
  )

  # Check results for setting without conditional independence.
  expect_equal(
    c(sens_results$ICA[1],
      sens_results$sp_rho[1],
      sens_results$c23[1]),
    c(0.995992037557, 0.995217038804, 0.436161760666)
  )

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
    cond_ind = TRUE,
    n_sim = 1,
    n_prec = 2e3,
    mutinfo_estimator = function(x, y) 1 - exp(-2 * stats::cor(x, y, method = "spearman")),
    restr_time = 2
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1],
                    sens_results$c23[1])
  check_vector = c(0.822189801988, 0.995997149542, 1.682036932130)
  expect_equal(output_vector, check_vector, tolerance = 1e-5)
})
