test_that("Delta method based on numerical derivatives works", {
  data("Ovarian")
  #For simplicity, data is not recoded to semi-competing risks format, but is
  #left in the composite event format.
  data = data.frame(Ovarian$Pfs,
                    Ovarian$Surv,
                    Ovarian$Treat,
                    Ovarian$PfsInd,
                    Ovarian$SurvInd)
  fitted_model = fit_model_SurvSurv(data = data,
                                    copula_family = "clayton",
                                    n_knots = 1)
  # The computation of the ICA is way to inaccurate here to combine it with
  # numerical differentiation. This test is added to check breaking changes.
  estimated_variance = delta_method_log_mutinfo(
    fitted_model = fitted_model,
    rotation_par_id = 0,
    rotation_par_unid = c(0,0,0,0),
    copula_par_unid = 1:4,
    n_prec = 1e3,
    composite = TRUE,
    seed = 1,
    eps = 1e-2
  )
  expect_equal(estimated_variance, matrix(0.215427099938, nrow = 1, ncol = 1))
})

test_that("Resampling method for the ICA based on fitted model summary information works", {
  data("Ovarian")
  #For simplicity, data is not recoded to semi-competing risks format, but is
  #left in the composite event format.
  data = data.frame(Ovarian$Pfs,
                    Ovarian$Surv,
                    Ovarian$Treat,
                    Ovarian$PfsInd,
                    Ovarian$SurvInd)
  fitted_model = fit_model_SurvSurv(data = data,
                                    copula_family = "clayton",
                                    n_knots = 1)
  # The computation of the ICA is way to inaccurate here to combine it with
  # numerical differentiation. This test is added to check breaking changes.
  estimated_ICAs = summary_level_bootstrap_ICA(
    fitted_model = fitted_model,
    rotation_par_id = 0,
    rotation_par_unid = c(0,0,0,0),
    copula_par_unid = 1:4,
    n_prec = 1e3,
    composite = TRUE,
    seed = 1,
    B = 10L
  )
  mean_ICA = mean(estimated_ICAs)
  sd_ICA = sd(estimated_ICAs)
  expect_equal(c(mean_ICA, sd_ICA), c(0.98869172852483, 0.00220452715409))
})
