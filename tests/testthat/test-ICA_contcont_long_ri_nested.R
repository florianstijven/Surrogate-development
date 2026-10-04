library(testthat)

test_that("ICA_contcont_long_ri_nested() works with exponential correlation", {
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

  times_S <- list(
    week1 = c(1),
    week1_2 = c(1, 2),
    week1_2_4 = c(1, 2, 4)
  )

  result <- ICA_contcont_long_ri_nested(
    V0 = V0,
    V1 = V1,
    D0 = D0,
    D1 = D1,
    times_T = c(1, 2, 4, 6, 8),
    times_S = times_S,
    correlation = "exponential",
    temporal_rho = 0.7117544,
    grid_V = seq(-1, 1, by = 0.5),
    grid_D = seq(-1, 1, by = 0.5),
    return_candidates = FALSE,
    tol = 1e-12
  )

  expect_s3_class(result, "ICA_contcont_long_ri_nested")
  expect_equal(result$correlation, "exponential")
  expect_equal(result$temporal_rho, 0.7117544)

  expect_true(is.matrix(result$R2_Lambda))
  expect_equal(dim(result$R2_Lambda), c(9, 3))
  expect_equal(colnames(result$R2_Lambda), names(times_S))

  expect_equal(result$Num.Pos.Def.V, 3)
  expect_equal(result$Num.Pos.Def.D, 3)
  expect_equal(result$Num.Pos.Def.Pairs, 9)

  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))

  result_summary <- summary_longitudinal_ica(
    list(R2_Lambda = result$R2_Lambda[, "week1_2_4"])
  )

  expect_s3_class(result_summary, "data.frame")
  expect_equal(nrow(result_summary), 1)
})

test_that("ICA_contcont_long_ri_nested() works with compound symmetry correlation", {
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

  times_S <- list(
    week1 = c(1),
    week1_2 = c(1, 2),
    week1_2_4 = c(1, 2, 4)
  )

  result <- ICA_contcont_long_ri_nested(
    V0 = V0,
    V1 = V1,
    D0 = D0,
    D1 = D1,
    times_T = c(1, 2, 4, 6, 8),
    times_S = times_S,
    correlation = "compound_symmetry",
    temporal_rho = 0.2,
    grid_V = seq(-1, 1, by = 0.5),
    grid_D = seq(-1, 1, by = 0.5)
  )

  expect_s3_class(result, "ICA_contcont_long_ri_nested")
  expect_equal(result$correlation, "compound_symmetry")
  expect_equal(result$temporal_rho, 0.2)

  expect_true(is.matrix(result$R2_Lambda))
  expect_equal(dim(result$R2_Lambda), c(result$Num.Pos.Def.Pairs, length(times_S)))
  expect_equal(colnames(result$R2_Lambda), names(times_S))

  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))
})

test_that("ICA_contcont_long_ri_nested() works with independence correlation", {
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

  times_S <- list(
    week1 = c(1),
    week1_2 = c(1, 2),
    week1_2_4 = c(1, 2, 4)
  )

  result <- ICA_contcont_long_ri_nested(
    V0 = V0,
    V1 = V1,
    D0 = D0,
    D1 = D1,
    times_T = c(1, 2, 4, 6, 8),
    times_S = times_S,
    correlation = "independence",
    grid_V = seq(-1, 1, by = 0.5),
    grid_D = seq(-1, 1, by = 0.5)
  )

  expect_s3_class(result, "ICA_contcont_long_ri_nested")
  expect_equal(result$correlation, "independence")
  expect_true(is.na(result$temporal_rho))

  expect_true(is.matrix(result$R2_Lambda))
  expect_equal(dim(result$R2_Lambda), c(result$Num.Pos.Def.Pairs, length(times_S)))
  expect_equal(colnames(result$R2_Lambda), names(times_S))

  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))
})

test_that("ICA_contcont_long_ri_nested() returns candidate grids when requested", {
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

  times_S <- list(
    week1 = c(1),
    week1_2 = c(1, 2),
    week1_2_4 = c(1, 2, 4)
  )

  result <- ICA_contcont_long_ri_nested(
    V0 = V0,
    V1 = V1,
    D0 = D0,
    D1 = D1,
    times_T = c(1, 2, 4, 6, 8),
    times_S = times_S,
    correlation = "exponential",
    temporal_rho = 0.7117544,
    grid_V = seq(-1, 1, by = 0.5),
    grid_D = seq(-1, 1, by = 0.5),
    return_candidates = TRUE
  )

  expect_s3_class(result, "ICA_contcont_long_ri_nested")
  expect_equal(result$Num.Pos.Def.V, nrow(result$V_candidates))
  expect_equal(result$Num.Pos.Def.D, nrow(result$D_candidates))
  expect_equal(result$Num.Pos.Def.Pairs, nrow(result$R2_Lambda))

  expect_s3_class(result$V_candidates, "data.frame")
  expect_s3_class(result$D_candidates, "data.frame")

  expect_true(all(c("c11", "c12", "c22", "rho2") %in% names(result$V_candidates)))
  expect_true(all(c("c11", "c12", "c22", "rho2") %in% names(result$D_candidates)))

  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))
})
