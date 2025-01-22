test_that(
  "estimate_mutual_information_OrdOrd() and compute_ICA_OrdOrd() work with ordinal endpoints with 3 categories",
  {
    copula_par = c(10, 2, 10, 2, 2, 2)
    rotation_par = rep(0, 6)
    n_prec = 1e4
    copula_family = "clayton"
    # Quantile functions.
    q_S0 = function(p) {
      ifelse(p < 0.4, 1, ifelse(p < 0.7, 2, 3))
    }
    q_S1 = function(p) {
      ifelse(p < 0.4, 1, ifelse(p < 0.7, 2, 3))
    }
    q_T0 = function(p) {
      ifelse(p < 0.3, 1, ifelse(p < 0.6, 2, 3))
    }
    q_T1 = function(p) {
      ifelse(p < 0.4, 1, ifelse(p < 0.7, 2, 3))
    }

    set.seed(1)
    Delta = sample_deltas_BinCont(
      copula_par,
      rotation_par,
      copula_family,
      n = 1e4,
      q_S0 = q_S0,
      q_S1 = q_S1,
      q_T0 = q_T0,
      q_T1 = q_T1,
      marginal_sp_rho = FALSE,
      setting = "SurvSurv",
      composite = FALSE
    )$Delta_dataframe

    mut_info = estimate_mutual_information_OrdOrd(Delta$DeltaS, Delta$DeltaT)
    ICA_estimate = estimate_ICA_OrdOrd(Delta$DeltaS, Delta$DeltaT)
    ICA_computed  = compute_ICA_OrdOrd(
      copula_par = copula_par,
      rotation_par = rotation_par,
      copula_family1 = copula_family,
      copula_family2 = copula_family,
      n_prec = 1e4,
      q_S0 = q_S0,
      q_S1 = q_S1,
      q_T0 = q_T0,
      q_T1 = q_T1,
      marginal_sp_rho = TRUE,
      seed = 1
    )
    # Computed mutual information is correct.
    expect_equal(mut_info, 0.287093723436, tolerance = 0.05)
    # Computed ICA is correct.
    expect_equal(ICA_estimate, 0.355527102471, tolerance = 0.05)
    # Tolerance should be sufficiently large because the VineCopula package, which
    # is used for sampling from D-vine copulas, will produce consistent samples
    # across package versions and platforms.

    # Estimated ICA and computed ICA are equal.
    expect_equal(ICA_estimate, ICA_computed[1], ignore_attr = "names")
  }
)
