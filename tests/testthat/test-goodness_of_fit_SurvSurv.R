test_that("mean_S_before_T() works well with Gaussian, Frank, Clayton, and Gumbel copulas",
          {
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
            expect_equal(output_vec,
                         c(
                           0.724653499621,
                           0.606863379335,
                           0.711095534584,
                           0.704344059670
                         ))
          })

test_that("plotting functions for GoF work as expected.",
          {
            # Load fitted copula model.
            fitted_gumbel = readRDS(test_path("fixtures", "ovarian-dvine-gumbel-scr.rds"))

            # Compare the plotted functions with a snapshot. Note that this
            # comparison is fragile, and unimportant changes could lead to a
            # failure of this test.
            vdiffr::expect_doppelganger(
              title = "mean surrogate before dying - gof",
              fig = function() {
                mean_S_before_T_plots_scr(fitted_model = fitted_gumbel,
                                          grid = seq(
                                            from = 1e-5,
                                            to = 2,
                                            length.out = 20
                                          ))
              }
            )
            vdiffr::expect_doppelganger(
              title = "marginal survival functions - gof",
              fig = function() {
                marginal_gof_plots_scr(fitted_model = fitted_gumbel,
                                       grid = seq(
                                         from = 1e-5,
                                         to = 2,
                                         length.out = 20
                                       ))
              }
            )
            vdiffr::expect_doppelganger(
              title = "probability of dying without progression - gof",
              fig = function() {
                prob_dying_without_progression_plots(fitted_model = fitted_gumbel,
                                                     grid = seq(
                                                       from = 1e-5,
                                                       to = 2,
                                                       length.out = 20
                                                     ))
              }
            )
          })
