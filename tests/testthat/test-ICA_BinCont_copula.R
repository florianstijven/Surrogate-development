# The test below is no longer used because the VineCopula package will not
# produce consistent samples across package versions and platforms.

# test_that("sample_dvine() works with clayton copula", {
#   copula_par = 1:6
#   rotation_par = rep(0, 6)
#   copula_family = "clayton"
#   set.seed(1)
#   U = sample_dvine(copula_par, rotation_par, copula_family, n = 1e2)
#   expected_values = c(0.99152565, 0.02532641, 0.52018781, 0.06633612)
#   output_values = c(U[1, 1], U[2, 2], U[3, 3], U[4, 4])
#   expect_equal(output_values,
#                expected_values)
# })

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
  )$Delta_dataframe
  # The VineCopula package cannot generate consistent samples across package
  # versions and platforms. So, we cannot construct robust tests for sampling.

  # expected_values = c(1.20168104, 0.99663129, 0L, 0L)
  # output_values = c(Delta$DeltaS[1:2], Delta$DeltaT[3:4])
  # expect_equal(output_values,
  #              expected_values)
})


# test_that("estimate_ICA_BinCont() works in a sample example", {
#   copula_par = 1:6
#   rotation_par = rep(0, 6)
#   copula_family = "clayton"
#   # Quantile functions.
#   q_S0 = function(x) {
#     qnorm(x)
#   }
#   q_S1 = function(x) {
#     qnorm(x, mean = 1)
#   }
#   set.seed(1)
#   Delta = sample_deltas_BinCont(
#     copula_par,
#     rotation_par,
#     copula_family,
#     n = 1e2,
#     q_S0 = q_S0,
#     q_S1 = q_S1
#   )$Delta_dataframe
#   ICA = estimate_ICA_BinCont(Delta$DeltaS, Delta$DeltaT)
#   expected_values = c(1.20168104, 0.99663129, 0L, 0L)
#   output_values = c(Delta$DeltaS[1:2], Delta$DeltaT[3:4])
#   expect_equal(output_values,
#                expected_values)
# })

test_that("compute_ICA_BinCont() works in a sample example", {
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
  ICA = compute_ICA_BinCont(
    copula_par = copula_par,
    rotation_par = rotation_par,
    copula_family1 = copula_family,
    copula_family2 = copula_family,
    n_prec = 1e3,
    q_S0 = q_S0,
    q_S1 = q_S1,
    seed = 1
  )
  expected_values = c(0.25048788, -0.53843674, 0.3952069, 0.74461876)
  output_values = ICA[1:4]
  expect_equal(output_values,
               expected_values,
               ignore_attr = "names",
               tolerance = 0.05)
  # Tolerance should be sufficiently large because the VineCopula package, which
  # is used for sampling from D-vine copulas, will not produce consistent
  # samples across package versions and platforms.
})

