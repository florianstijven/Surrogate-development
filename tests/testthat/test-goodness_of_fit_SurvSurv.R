test_that("mean_S_before_T() works well with Gaussian, Frank, Clayton, and Gumbel copulas",
          {
            data("Ovarian")
            data_scr = data.frame(
              ttp = Ovarian$Pfs,
              os = Ovarian$Surv,
              treat = Ovarian$Treat,
              ttp_ind = ifelse(
                Ovarian$Pfs == Ovarian$Surv &
                  Ovarian$SurvInd == 1,
                0,
                Ovarian$PfsInd
              ),
              os_ind = Ovarian$SurvInd
            )
            # Load fitted copula models.
            fitted_clayton = readRDS(test_path("fixtures", "ovarian-dvine-clayton-scr.rds"))
            fitted_gaussian = readRDS(test_path("fixtures", "ovarian-dvine-gaussian-scr.rds"))
            fitted_frank = readRDS(test_path("fixtures", "ovarian-dvine-frank-scr.rds"))
            fitted_gumbel = readRDS(test_path("fixtures", "ovarian-dvine-gumbel-scr.rds"))


            clayton0 = mean_S_before_T(1.1, fitted_model = fitted_clayton, treated = 0)
            gaussian0 = mean_S_before_T(1.1, fitted_model = fitted_gaussian, treated = 0)
            frank1 = mean_S_before_T(1.1, fitted_model = fitted_frank, treated = 1)
            gumbel1 = mean_S_before_T(1.1, fitted_model = fitted_gumbel, treated = 1)

            output_vec = c(clayton0,
                           gaussian0,
                           frank1,
                           gumbel1)
            expect_equal(output_vec, c(0.724653499621, 0.606863379335, 0.711095534584, 0.704344059670))
          })
