library(testthat)

test_that("ICA_contcont_long_ri() works with exponential correlation", {
  V0 <- matrix(c(0.02579325, 0.02514841,
                 0.02514841, 0.02872256),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  V1 <- matrix(c(0.02907611, 0.02962312,
                 0.02962312, 0.03400179),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  D0 <- matrix(c(0.02587092, 0.02458997,
                 0.02458997, 0.02405581),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  D1 <- matrix(c(0.03649246, 0.03225987,
                 0.03225987, 0.02897932),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  result <- ICA_contcont_long_ri(
    V0 = V0,
    V1 = V1,
    D0 = D0,
    D1 = D1,
    times = c(1, 2, 4, 6, 8),
    correlation = "exponential",
    temporal_rho = 0.7117544,
    grid_V = seq(-1, 1, by = 1),
    grid_D = seq(-1, 1, by = 1),
    return_candidates = FALSE,
    tol = 1e-12
  )

  expect_s3_class(result, "ICA_contcont_long_ri")
  expect_equal(result$correlation, "exponential")
  expect_equal(result$temporal_rho, 0.7117544)
  expect_equal(result$Num.Pos.Def.V, 1)
  expect_equal(result$Num.Pos.Def.D, 1)
  expect_equal(result$Num.Pos.Def.Pairs, 1)
  expect_equal(length(result$R2_Lambda), 1)
  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))

  result_summary <- summary_longitudinal_ica(result)
  expect_s3_class(result_summary, "data.frame")
  expect_equal(nrow(result_summary), 1)
})


test_that("ICA_contcont_long_ri() works with independence correlation", {
  V0 <- matrix(c(0.02579325, 0.02514841,
                 0.02514841, 0.02872256),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  V1 <- matrix(c(0.02907611, 0.02962312,
                 0.02962312, 0.03400179),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  D0 <- matrix(c(0.02587092, 0.02458997,
                 0.02458997, 0.02405581),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  D1 <- matrix(c(0.03649246, 0.03225987,
                 0.03225987, 0.02897932),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  result <- ICA_contcont_long_ri(
    V0 = V0,
    V1 = V1,
    D0 = D0,
    D1 = D1,
    times = c(1, 2, 4, 6, 8),
    correlation = "independence",
    grid_V = seq(-1, 1, by = 1),
    grid_D = seq(-1, 1, by = 1)
  )

  expect_s3_class(result, "ICA_contcont_long_ri")
  expect_equal(result$correlation, "independence")
  expect_true(is.na(result$temporal_rho))
  expect_equal(length(result$R2_Lambda), result$Num.Pos.Def.Pairs)
  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))
})

test_that("ICA_contcont_long_ri() returns candidate grids when requested", {
  V0 <- matrix(c(0.02579325, 0.02514841,
                 0.02514841, 0.02872256),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  V1 <- matrix(c(0.02907611, 0.02962312,
                 0.02962312, 0.03400179),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  D0 <- matrix(c(0.02587092, 0.02458997,
                 0.02458997, 0.02405581),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  D1 <- matrix(c(0.03649246, 0.03225987,
                 0.03225987, 0.02897932),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  result <- ICA_contcont_long_ri(
    V0 = V0,
    V1 = V1,
    D0 = D0,
    D1 = D1,
    times = c(1, 2, 4, 6, 8),
    correlation = "exponential",
    temporal_rho = 0.7117544,
    grid_V = seq(-1, 1, by = 1),
    grid_D = seq(-1, 1, by = 1),
    return_candidates = TRUE
  )

  expect_s3_class(result, "ICA_contcont_long_ri")
  expect_equal(result$Num.Pos.Def.V, nrow(result$V_candidates))
  expect_equal(result$Num.Pos.Def.D, nrow(result$D_candidates))
  expect_equal(result$Num.Pos.Def.Pairs, length(result$R2_Lambda))

  expect_s3_class(result$V_candidates, "data.frame")
  expect_s3_class(result$D_candidates, "data.frame")

  expect_true(all(c("c11", "c12", "c22", "rho2") %in% names(result$V_candidates)))
  expect_true(all(c("c11", "c12", "c22", "rho2") %in% names(result$D_candidates)))

  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))
})
