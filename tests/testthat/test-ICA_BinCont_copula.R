test_that("sample_dvine() works with clayton copula", {
  copula_par = 1:6
  rotation_par = rep(0, 6)
  copula_family = "clayton"
  set.seed(1)
  U = sample_dvine(copula_par, rotation_par, copula_family, n = 1e2)
})
