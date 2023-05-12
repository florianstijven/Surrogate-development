test_that("sample_dvine() works with clayton copula", {
  copula_par = 1:6
  rotation_par = rep(0, 6)
  copula_family = "clayton"
  set.seed(1)
  U = sample_dvine(copula_par, rotation_par, copula_family, n = 1e2)
  expected_values = c(0.98379818, 0.02441270, 0.44578131, 0.06633612)
  output_values = c(U[1, 1], U[2, 2], U[3, 3], U[4, 4])
  expect_equal(output_values,
               expected_values)
})

test_that("sample_deltas_BinCont() works with clayton copula", {
  copula_par = 1:6
  rotation_par = rep(0, 6)
  copula_family = "clayton"
  # Quantile functions.
  q_S0 = function(x) {
    qnorm(x)
  }
  q_S1 = function(x) {
    qnorm(x, mean = 1)
  }
  set.seed(1)
  Delta = sample_deltas_BinCont(
    copula_par,
    rotation_par,
    copula_family,
    n = 1e2,
    q_S0 = q_S0,
    q_S1 = q_S1
  )
  expected_values = c(1.03557703, 1.07562675, 1L, 0L)
  output_values = c(Delta$DeltaS[1:2], Delta$DeltaT[3:4])
  expect_equal(output_values,
               expected_values)
})


test_that("estimate_ICA_BinCont() works in a sample example", {
  copula_par = 1:6
  rotation_par = rep(0, 6)
  copula_family = "clayton"
  # Quantile functions.
  q_S0 = function(x) {
    qnorm(x)
  }
  q_S1 = function(x) {
    qnorm(x, mean = 1)
  }
  set.seed(1)
  Delta = sample_deltas_BinCont(
    copula_par,
    rotation_par,
    copula_family,
    n = 1e2,
    q_S0 = q_S0,
    q_S1 = q_S1
  )
  ICA = estimate_ICA_BinCont(Delta$DeltaS, Delta$DeltaT)
  expected_values = c(1.03557703, 1.07562675, 1L, 0L)
  output_values = c(Delta$DeltaS[1:2], Delta$DeltaT[3:4])
  expect_equal(output_values,
               expected_values)
})
