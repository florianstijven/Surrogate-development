library(testthat)

test_that("ICA_contcont_long_shared_ri() works with exponential correlation", {
  V0 <- matrix(c(0.03025730, 0.02925077,
                 0.02925077, 0.03327464),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  V1 <- matrix(c(0.03397819, 0.03390776,
                 0.03390776, 0.03849923),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  result <- ICA_contcont_long_shared_ri(
    V0 = V0,
    V1 = V1,
    D0_variance = 0.02160313,
    D1_variance = 0.03115772,
    times = c(1, 2, 4, 6, 8),
    correlation = "exponential",
    temporal_rho = 0.7629443,
    grid_V = seq(-1, 1, by = 0.5),
    grid_D = seq(-1, 1, by = 0.5),
    return_candidates = FALSE,
    tol = 1e-12
  )

  expect_s3_class(result, "ICA_contcont_long_shared_ri")
  expect_equal(result$correlation, "exponential")
  expect_equal(result$temporal_rho, 0.7629443)
  expect_equal(result$Num.Pos.Def.V, 3)
  expect_equal(result$Num.Pos.Def.D, 3)
  expect_equal(result$Num.Pos.Def.Pairs, 9)
  expect_equal(length(result$R2_Lambda), 9)
  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))

  result_summary <- summary_longitudinal_ica(result)
  expect_s3_class(result_summary, "data.frame")
  expect_equal(nrow(result_summary), 1)
})

test_that("ICA_contcont_long_shared_ri() works with compound symmetry correlation", {
  V0 <- matrix(c(0.03025730, 0.02925077,
                 0.02925077, 0.03327464),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  V1 <- matrix(c(0.03397819, 0.03390776,
                 0.03390776, 0.03849923),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  result <- ICA_contcont_long_shared_ri(
    V0 = V0,
    V1 = V1,
    D0_variance = 0.02160313,
    D1_variance = 0.03115772,
    times = c(1, 2, 4, 6, 8),
    correlation = "compound_symmetry",
    temporal_rho = 0.2,
    grid_V = seq(-1, 1, by = 0.5),
    grid_D = seq(-1, 1, by = 0.5)
  )

  expect_s3_class(result, "ICA_contcont_long_shared_ri")
  expect_equal(result$correlation, "compound_symmetry")
  expect_equal(result$temporal_rho, 0.2)
  expect_equal(length(result$R2_Lambda), result$Num.Pos.Def.Pairs)
  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))
})

test_that("ICA_contcont_long_shared_ri() works with independence correlation", {
  V0 <- matrix(c(0.03025730, 0.02925077,
                 0.02925077, 0.03327464),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  V1 <- matrix(c(0.03397819, 0.03390776,
                 0.03390776, 0.03849923),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  result <- ICA_contcont_long_shared_ri(
    V0 = V0,
    V1 = V1,
    D0_variance = 0.02160313,
    D1_variance = 0.03115772,
    times = c(1, 2, 4, 6, 8),
    correlation = "independence",
    grid_V = seq(-1, 1, by = 0.5),
    grid_D = seq(-1, 1, by = 0.5)
  )

  expect_s3_class(result, "ICA_contcont_long_shared_ri")
  expect_equal(result$correlation, "independence")
  expect_true(is.na(result$temporal_rho))
  expect_equal(length(result$R2_Lambda), result$Num.Pos.Def.Pairs)
  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))
})

test_that("ICA_contcont_long_shared_ri() returns candidate grids when requested", {
  V0 <- matrix(c(0.03025730, 0.02925077,
                 0.02925077, 0.03327464),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  V1 <- matrix(c(0.03397819, 0.03390776,
                 0.03390776, 0.03849923),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("T", "S"), c("T", "S")))

  result <- ICA_contcont_long_shared_ri(
    V0 = V0,
    V1 = V1,
    D0_variance = 0.02160313,
    D1_variance = 0.03115772,
    times = c(1, 2, 4, 6, 8),
    correlation = "exponential",
    temporal_rho = 0.7629443,
    grid_V = seq(-1, 1, by = 0.5),
    grid_D = seq(-1, 1, by = 0.5),
    return_candidates = TRUE
  )

  expect_s3_class(result, "ICA_contcont_long_shared_ri")
  expect_equal(result$Num.Pos.Def.V, nrow(result$V_candidates))
  expect_equal(result$Num.Pos.Def.D, nrow(result$D_candidates))
  expect_equal(result$Num.Pos.Def.Pairs, length(result$R2_Lambda))

  expect_s3_class(result$V_candidates, "data.frame")
  expect_s3_class(result$D_candidates, "data.frame")

  expect_true(all(c("c11", "c12", "c22", "rho2") %in% names(result$V_candidates)))
  expect_true(all(c("rho_b0b1", "d_delta") %in% names(result$D_candidates)))
  expect_true(all(is.finite(result$D_candidates$d_delta)))
  expect_true(all(result$D_candidates$d_delta >= 0))

  expect_true(all(result$R2_Lambda >= 0))
  expect_true(all(result$R2_Lambda <= 1))
})
