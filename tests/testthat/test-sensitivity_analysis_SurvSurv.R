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
  check_vector = c(0.987833035047, 0.984428257713, 1.374917942896)
  expect_equal(output_vector, check_vector)
})

test_that("sensitivity_analysis_SurvSurv_copula() works on 2 cores with Clayton copula", {
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
  check_vector = c(0.987833035047, 0.984428257713, 1.374917942896)
  expect_equal(output_vector, check_vector)
})

test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Gaussian copula", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-gaussian.rds"))
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
  check_vector = c(0.957439366743, 0.962676202705, 0.152435752847)
  expect_equal(output_vector, check_vector)
})

test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Frank copula", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-frank.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 1,
    n_prec = 2e3
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1],
                    sens_results$c23[1])
  check_vector = c(0.990380494599, 0.986181510545, -3.171396109738)
  expect_equal(output_vector, check_vector, tolerance = 1e-5)
})

test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Frank copula and four different unidentifiable copulas", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-frank.rds"))
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    fitted_model,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 1,
    n_prec = 2e3,
    copula_family2 = c("clayton", "frank", "gaussian", "gumbel")
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$sp_rho[1],
                    sens_results$c23[1])
  check_vector = c(0.991389516772, 0.988913341228, 0.436161760666)
  expect_equal(output_vector, check_vector, tolerance = 1e-5)
})
