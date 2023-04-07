test_that("SurvivalSurvival model works with Clayton copula for Ovarian data", {
  if (require(Surrogate)) {
    data("Ovarian")
    #For simplicity, data is not recoded to semi-competing risks format, but is
    #left in the composite event format.
    data = data.frame(Ovarian$Pfs,
                      Ovarian$Surv,
                      Ovarian$Treat,
                      Ovarian$PfsInd,
                      Ovarian$SurvInd)
    fitted_model = Surrogate::fit_model_SurvSurv(data = data,
                                  copula_family = "clayton",
                                  nknots = 1)
  }
  log_lik_fitted = c(fitted_model$log_lik0, fitted_model$log_lik1)
  expect_equal(c(362.79316, 283.69633), log_lik_fitted)
})

test_that("SurvivalSurvival model works with Gaussian copula for Ovarian data", {
  if (require(Surrogate)) {
    data("Ovarian")
    #For simplicity, data is not recoded to semi-competing risks format, but is
    #left in the composite event format.
    data = data.frame(Ovarian$Pfs,
                      Ovarian$Surv,
                      Ovarian$Treat,
                      Ovarian$PfsInd,
                      Ovarian$SurvInd)
    fitted_model = Surrogate::fit_model_SurvSurv(data = data,
                                                 copula_family = "gaussian",
                                                 nknots = 1)
  }
  log_lik_fitted = c(fitted_model$log_lik0, fitted_model$log_lik1)
  expect_equal(c(319.42767, 256.69588), log_lik_fitted)
})

test_that("SurvivalSurvival model works with Frank copula for Ovarian data", {
  if (require(Surrogate)) {
    data("Ovarian")
    #For simplicity, data is not recoded to semi-competing risks format, but is
    #left in the composite event format.
    data = data.frame(Ovarian$Pfs,
                      Ovarian$Surv,
                      Ovarian$Treat,
                      Ovarian$PfsInd,
                      Ovarian$SurvInd)
    fitted_model = Surrogate::fit_model_SurvSurv(data = data,
                                                 copula_family = "frank",
                                                 nknots = 1)
  }
  log_lik_fitted = c(fitted_model$log_lik0, fitted_model$log_lik1)
  expect_equal(c(355.93429, 291.65151), log_lik_fitted)
})

test_that("SurvivalSurvival model works with Gumbel copula for Ovarian data", {
  if (require(Surrogate)) {
    data("Ovarian")
    #For simplicity, data is not recoded to semi-competing risks format, but is
    #left in the composite event format.
    data = data.frame(Ovarian$Pfs,
                      Ovarian$Surv,
                      Ovarian$Treat,
                      Ovarian$PfsInd,
                      Ovarian$SurvInd)
    fitted_model = Surrogate::fit_model_SurvSurv(data = data,
                                                 copula_family = "gumbel",
                                                 nknots = 1)
  }
  log_lik_fitted = c(fitted_model$log_lik0, fitted_model$log_lik1)
  expect_equal(c(338.562039, 274.257801), log_lik_fitted)
})
