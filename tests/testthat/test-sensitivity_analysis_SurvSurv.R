test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Clayton copula", {
  data("Ovarian")
  # For simplicity, data is not recoded to semi-competing risks format, but the
  # data are left in the composite event format.
  data = data.frame(
    Ovarian$Pfs,
    Ovarian$Surv,
    Ovarian$Treat,
    Ovarian$PfsInd,
    Ovarian$SurvInd
  )
  ovarian_fitted =
      fit_model_SurvSurv(data = data,
                         copula_family = "clayton",
                         nknots = 1)
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    ovarian_fitted,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 5,
    n_prec = 500,
    minfo_prec = 2e3
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$c23[3])
  check_vector = c(0.95596036, 1.37491794289595)
  expect_equal(output_vector, check_vector)
})

test_that("sensitivity_analysis_SurvSurv_copula() works on 2 cores with Clayton copula", {
  data("Ovarian")
  # For simplicity, data is not recoded to semi-competing risks format, but the
  # data are left in the composite event format.
  data = data.frame(
    Ovarian$Pfs,
    Ovarian$Surv,
    Ovarian$Treat,
    Ovarian$PfsInd,
    Ovarian$SurvInd
  )
  ovarian_fitted =
    fit_model_SurvSurv(data = data,
                       copula_family = "clayton",
                       nknots = 1)
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    ovarian_fitted,
    composite = TRUE,
    cond_ind = FALSE,
    n_sim = 5,
    n_prec = 500,
    minfo_prec = 500,
    marg_association = TRUE,
    ncores = 2
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$c23[3])
  check_vector = c(0.96872392, 1.37491794289595)
  expect_equal(output_vector, check_vector)
})

test_that("sensitivity_analysis_SurvSurv_copula() works on a single core with Gaussian copula", {
  data("Ovarian")
  # For simplicity, data is not recoded to semi-competing risks format, but the
  # data are left in the composite event format.
  data = data.frame(
    Ovarian$Pfs,
    Ovarian$Surv,
    Ovarian$Treat,
    Ovarian$PfsInd,
    Ovarian$SurvInd
  )
  ovarian_fitted =
    fit_model_SurvSurv(data = data,
                       copula_family = "gaussian",
                       nknots = 1)
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = sensitivity_analysis_SurvSurv_copula(
    ovarian_fitted,
    composite = TRUE,
    cond_ind = TRUE,
    n_sim = 5,
    n_prec = 500,
    minfo_prec = 2e3
  )
  output_vector = c(sens_results$ICA[1],
                    sens_results$c23[3])
  check_vector = c(0.93654753, 0.15243575)
  expect_equal(output_vector, check_vector)
})
