test_that("sensitivity_analysis_BinCont_copula() works in a simplified setting",
          {
            # Load fitted copula model.
            fitted_model = readRDS(test_path("fixtures", "schizo-dvine-clayton.rds"))
            n_sim = 2
            set.seed(1)
            sens_results = sensitivity_analysis_BinCont_copula(
              fitted_model,
              n_sim,
              cond_ind = TRUE,
              lower = c(-1,-1,-1,-1),
              upper = c(1, 1, 1, 1),
              n_prec = 1e3
            )
            expect_equal(
              sens_results[,1],
              c(0.43244243, 0.30876395),
              tolerance = 1e-3
            )
          })

test_that("sensitivity_analysis_BinCont_copula() works in a simplified setting with limits in the sensitivity analysis",
          {
            # Load fitted copula model.
            fitted_model = readRDS(test_path("fixtures", "schizo-dvine-clayton.rds"))

            n_sim = 2
            set.seed(1)
            sens_results = sensitivity_analysis_BinCont_copula(
              fitted_model,
              n_sim,
              cond_ind = TRUE,
              lower = c(0, 0, 0, 0),
              upper = c(0.95, 0.95, 0.95, 0.95),
              n_prec = 1e3
            )
            expect_equal(
              sens_results[,1],
              c(0.41079208, 0.3107483),
              tolerance = 1e-3
            )
          })
