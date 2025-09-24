test_that("compute_ICA_ContCont() works", {
  copula_par = c(5, 2, 5, 2, 2, 2)
  rotation_par = rep(0, 6)
  copula_family1 = "clayton"
  copula_family2 = copula_family1
  n_prec = 5e3
  q_S0 = qnorm
  q_T0 = function(p) qnorm(p, mean = 1)
  q_S1 = qnorm
  q_T1 = function(p) qnorm(p, mean = 1)
  marginal_sp_rho = TRUE
  seed = 1
  # BUG TO FIX REGARDING MUTUAL INFORMATION ESTIMATOR
  output_vector = compute_ICA_ContCont(copula_par = copula_par,
                                       rotation_par = rotation_par,
                                       copula_family1 = copula_family1,
                                       copula_family2 = copula_family1,
                                       n_prec = n_prec,
                                       q_S0 = q_S0,
                                       q_T0 = q_T0,
                                       q_S1 = q_S1,
                                       q_T1 = q_T1,
                                       marginal_sp_rho = marginal_sp_rho,
                                       seed = seed)
  output_vector = unname(output_vector[1:4])
  check_vector = c(0.318517383927, 0.384746331870, 0.885863084555, 0.826571557479)

  expect_equal(output_vector, check_vector, tolerance = 0.05)
  # Tolerance should be sufficiently large because the VineCopula package, which
  # is used for sampling from D-vine copulas, will produce consistent samples
  # across package versions and platforms.
})
