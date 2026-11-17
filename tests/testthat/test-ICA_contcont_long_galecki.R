library(testthat)

test_that("ICA_contcont_long_galecki() works with a sensitivity grid", {
  V0 <- matrix(c(0.04513345, 0.04375355,
                 0.04375355, 0.04866820),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  V1 <- matrix(c(0.05231045, 0.05177398,
                 0.05177398, 0.05702684),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  result <- ICA_contcont_long_galecki(
    V0 = V0,
    V1 = V1,
    p = 5,
    grid_V = seq(-1, 1, by = 1),
    return_candidates = FALSE,
    tol = 1e-12
  )

  expect_s3_class(result, "ICA_contcont_long_galecki")
  expect_equal(result$Total.Num.Matrices, 81)
  expect_equal(result$Num.Pos.Def.V, 1)
  expect_equal(length(result$R2_Lambda), 1)
  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))

  result_summary <- summary_longitudinal_ica(result)

  expect_s3_class(result_summary, "data.frame")
  expect_equal(nrow(result_summary), 1)
})
