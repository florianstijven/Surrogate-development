test_that("compute_ICA_SurvSurv works", {
  copula_par = c(5, 2, 5, 2, 2, 2)
  rotation_par = rep(0, 6)
  copula_family1 = "clayton"
  copula_family2 = copula_family1
  n_prec = 5e3
  q_S0 = qnorm
  q_T0 = function(p) qnorm(p, mean = 1)
  q_S1 = qnorm
  q_T1 = function(p) qnorm(p, mean = 1)
  composite = TRUE
  marginal_sp_rho = TRUE
  seed = 1
  # BUG TO FIX REGARDING MUTUAL INFORMATION ESTIMATOR
  output_vector = compute_ICA_SurvSurv(copula_par,
                       rotation_par,
                       copula_family1,
                       copula_family2 = copula_family1,
                       n_prec,
                       q_S0,
                       q_T0,
                       q_S1,
                       q_T1,
                       composite,
                       marginal_sp_rho,
                       seed)
  output_vector = unname(output_vector[1:4])
  check_vector = c(0.329931193410, 0.385933009133, 0.890199861304, 0.835027399097)
  expect_equal(output_vector, check_vector, tolerance = 0.05)
})

test_that("compute_ICA_SurvSurv works for different unidentifiable copulas", {
  copula_par = 1:6
  copula_par[5] = 0.8
  rotation_par = rep(0, 6)
  copula_family1 = "frank"
  copula_family2 = c("clayton", "frank", "gaussian", "gumbel")
  n_prec = 1e4
  q_S0 = qnorm
  q_T0 = qnorm
  q_S1 = qnorm
  q_T1 = qnorm
  composite = TRUE
  marginal_sp_rho = TRUE
  seed = 1
  output_vector = compute_ICA_SurvSurv(copula_par,
                                       rotation_par,
                                       copula_family1,
                                       copula_family2,
                                       n_prec,
                                       q_S0,
                                       q_T0,
                                       q_S1,
                                       q_T1,
                                       composite,
                                       marginal_sp_rho,
                                       seed)
  output_vector = unname(output_vector[1:4])
  check_vector = c(0.440951593781, 0.310282184611, 0.676431933056, 0.511494812175)
  expect_equal(output_vector, check_vector, tolerance = 0.05)
})
