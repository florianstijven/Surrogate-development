test_that("SurvivalSurvival model works with Clayton copula for Ovarian data",
          {
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
            log_lik_fitted = c(logLik(fitted_model$fit_0), logLik(fitted_model$fit_1))
            expect_equal(log_lik_fitted, c(362.79316, 283.69633), ignore_attr = "df")
          })

test_that("SurvivalSurvival model works with Gaussian copula for Ovarian data",
          {
            data("Ovarian")
            #For simplicity, data is not recoded to semi-competing risks format, but is
            #left in the composite event format.
            data = data.frame(Ovarian$Pfs,
                              Ovarian$Surv,
                              Ovarian$Treat,
                              Ovarian$PfsInd,
                              Ovarian$SurvInd)
            fitted_model = fit_model_SurvSurv(data = data,
                                              copula_family = "gaussian",
                                              n_knots = 2)
            log_lik_fitted = c(logLik(fitted_model$fit_0), logLik(fitted_model$fit_1))
            expect_equal(log_lik_fitted, c(336.84319, 274.04324), ignore_attr = "df")
          })

test_that("SurvivalSurvival model works with Frank copula for Ovarian data",
          {
            data("Ovarian")
            #For simplicity, data is not recoded to semi-competing risks format, but is
            #left in the composite event format.
            data = data.frame(Ovarian$Pfs,
                              Ovarian$Surv,
                              Ovarian$Treat,
                              Ovarian$PfsInd,
                              Ovarian$SurvInd)
            fitted_model = fit_model_SurvSurv(data = data,
                                              copula_family = "frank",
                                              n_knots = 1,
                                              method = "BFGS")
            log_lik_fitted = c(logLik(fitted_model$fit_0), logLik(fitted_model$fit_1))
            expect_equal(log_lik_fitted, c(355.93244, 288.48185),
                         ignore_attr = "df",
                         # Tolerance is very high at the moment. This is because
                         # there are non-trivial numerical differences depending
                         # on in which environment the tests are run.
                         tolerance = 1)
          })

test_that("SurvivalSurvival model works with Gumbel copula for Ovarian data",
          {
            data("Ovarian")
            #For simplicity, data is not recoded to semi-competing risks format, but is
            #left in the composite event format.
            data = data.frame(Ovarian$Pfs,
                              Ovarian$Surv,
                              Ovarian$Treat,
                              Ovarian$PfsInd,
                              Ovarian$SurvInd)
            fitted_model = fit_model_SurvSurv(data = data,
                                              copula_family = "gumbel",
                                              n_knots = 1)
            log_lik_fitted = c(logLik(fitted_model$fit_0), logLik(fitted_model$fit_1))
            expect_equal(log_lik_fitted, c(338.562039, 274.257801), ignore_attr = "df")
          })
