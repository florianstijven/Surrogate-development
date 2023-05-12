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
