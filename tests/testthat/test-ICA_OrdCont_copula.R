test_that("estimate_mutual_information_OrdCont() and estimate_ICA_OrdCont() work in a sample example", {
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
  mut_info = estimate_mutual_information_OrdCont(Delta$DeltaS, Delta$DeltaT)
  ICA = estimate_ICA_OrdCont(Delta$DeltaS, Delta$DeltaT)
  # Computed mutual information and ICA equal expected values.
  expect_equal(c(mut_info, ICA),
               c(0.138087089646, 0.224333940442),
               tolerance = 0.3)
  # Computed ICA equals the value returned by
  # estimate_mutual_information_BinCont() and estimate_ICA_BinCont().
  expect_equal(c(mut_info, ICA),
               c(
                 estimate_mutual_information_BinCont(Delta$DeltaS, Delta$DeltaT),
                 estimate_ICA_BinCont(Delta$DeltaS, Delta$DeltaT)
               ))
})

test_that("compute_ICA_OrdCont() works in with ordinal endpoint with 3 categories", {
  copula_par = 1:6
  rotation_par = rep(0, 6)
  n_prec = 1e4
  copula_family = "clayton"
  # Quantile functions.
  q_S0 = function(x) {
    qnorm(x)
  }
  q_S1 = function(x) {
    qnorm(x, mean = 1)
  }
  q_T0 = function(p) {
    ifelse(p < 0.3, 1, ifelse(p < 0.6, 2, 3))
  }
  q_T1 = function(p) {
    ifelse(p < 0.4, 1, ifelse(p < 0.7, 2, 3))
  }

  ICA_marg_true = compute_ICA_OrdCont(
    copula_par = copula_par,
    rotation_par = rotation_par,
    copula_family1 = copula_family,
    copula_family2 = copula_family,
    n_prec = n_prec,
    q_S0 = q_S0,
    q_T0 = q_T0,
    q_S1 = q_S1,
    q_T1 = q_T1,
    marginal_sp_rho = TRUE,
    seed = 1
  )

  ICA_marg_false = compute_ICA_OrdCont(
    copula_par = copula_par,
    rotation_par = rotation_par,
    copula_family1 = copula_family,
    copula_family2 = copula_family,
    n_prec = n_prec,
    q_S0 = q_S0,
    q_T0 = q_T0,
    q_S1 = q_S1,
    q_T1 = q_T1,
    marginal_sp_rho = FALSE,
    seed = 1
  )
  # Computed ICA equals expected value.
  expect_equal(
    ICA_marg_true,
    c(
      0.300182238020, -0.689622488197, 0.437928990939, 0.815037420952, 0.582869469012, 0.670079357825, 0.864715817871, 0.717112698175
    ),
    ignore_attr = "names",
    tolerance = 0.1
  )
  # Tolerance should be sufficiently large because the VineCopula package, which
  # is used for sampling from D-vine copulas, will produce consistent samples
  # across package versions and platforms.

  # Computed ICA and Spearman's rho values do not depend on marginal_sp_rho
  # option.
  expect_equal(
    ICA_marg_true[1:2],
    ICA_marg_false[1:2]
  )
})


