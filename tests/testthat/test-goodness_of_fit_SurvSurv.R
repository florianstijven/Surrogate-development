test_that("mean_S_before_T() works well with Gaussian, Frank, Clayton, and Gumbel copulas",
          {
            # Load fitted copula models.
            fitted_gaussian = readRDS(test_path("fixtures", "ovarian-dvine-gaussian-scr.rds"))

            gaussian0 = mean_S_before_T(1.1, fitted_model = fitted_gaussian, treated = 0)

            expect_equal(gaussian0, 0.604292841017)
          })

test_that("plotting functions for GoF work as expected.",
          {
            # Load fitted copula model.
            fitted_gumbel = readRDS(test_path("fixtures", "ovarian-dvine-gaussian-scr.rds"))

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
              title = "prob of dying wo progression - gof",
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
