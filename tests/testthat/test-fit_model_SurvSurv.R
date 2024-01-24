test_that("SurvivalSurvival model works with Clayton and Gaussian copula for Ovarian data",
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
                                              copula_family = c("clayton", "gaussian"),
                                              n_knots = 1)
            log_lik_fitted = c(logLik(fitted_model$fit_0), logLik(fitted_model$fit_1))
            expect_equal(log_lik_fitted, c(362.793159888, 256.695876425), ignore_attr = "df")
          })

test_that("SurvivalSurvival model works with Frank and Gumbel copula for Ovarian data in semi-competing risks format",
          {
            # Load Ovarian data in semi-competing risks format.
            data = readRDS(test_path("fixtures", "ovarian-data-scr.rds"))

            fitted_model = fit_model_SurvSurv(data = data,
                                              copula_family = c("frank", "gumbel"),
                                              n_knots = 2)
            log_lik_fitted = c(logLik(fitted_model$fit_0), logLik(fitted_model$fit_1))
            expect_equal(log_lik_fitted, c(123.9074566805, 21.7039951207), ignore_attr = "df")
          })


