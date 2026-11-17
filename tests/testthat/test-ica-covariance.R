test_that("determinant ICA agrees with squared correlation in dimension one", {
  out <- R2Lambda_from_covariance(
    Sigma_TT = matrix(1, 1, 1),
    Sigma_SS = matrix(1, 1, 1),
    Sigma_TS = matrix(0.5, 1, 1)
  )
  expect_equal(out, 0.25, tolerance = 1e-12)
})

test_that("independent covariance blocks give zero ICA", {
  out <- R2Lambda_from_covariance(
    Sigma_TT = diag(2),
    Sigma_SS = diag(3),
    Sigma_TS = matrix(0, 2, 3)
  )
  expect_equal(out, 0, tolerance = 1e-12)
})

test_that("endpoint-RI update agrees with direct calculation", {
  times <- list(T = c(1, 2), S = c(1, 2), common = TRUE)
  R0 <- 0.5 ^ abs(outer(c(1, 2), c(1, 2), "-"))
  R <- list(TT = R0, SS = R0, TS = R0)
  U <- data.frame(c11 = 1, c12 = 0.3, c22 = 1, rho2 = 0.09)
  DD <- data.frame(c11 = 1, c12 = 0.7, c22 = 1)

  fast <- Surrogate:::.surr_endpoint_ri_ica(
    U, DD, R, times, tol = 1e-12
  )

  J <- matrix(1, 2, 2)
  full <- rbind(
    cbind(R0 + J, 0.3 * R0 + 0.7 * J),
    cbind(0.3 * R0 + 0.7 * J, R0 + J)
  )
  direct <- R2Lambda_from_covariance(
    full[1:2, 1:2],
    full[3:4, 3:4],
    full[1:2, 3:4]
  )
  expect_equal(fast, direct, tolerance = 1e-10)
})
