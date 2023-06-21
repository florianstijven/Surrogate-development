test_that("compute_ICA_SurvSurv works", {
  copula_par = 1:6
  rotation_par = rep(0, 6)
  copula_family1 = "clayton"
  copula_family2 = copula_family1
  n_prec = 1e3
  minfo_prec = 1e3
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
                       copula_family2 = copula_family1,
                       n_prec,
                       minfo_prec,
                       q_S0,
                       q_T0,
                       q_S1,
                       q_T1,
                       composite,
                       marginal_sp_rho,
                       seed)
  output_vector = unname(output_vector[1:4])
  check_vector = c(0.46076498, 0.57390746, 0.80508611, 0.69431985)
  expect_equal(output_vector, check_vector)
})
