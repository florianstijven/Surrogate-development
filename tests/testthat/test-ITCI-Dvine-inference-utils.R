test_that("Delta method based on numerical derivatives works", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # The computation of the ICA is way to inaccurate here to combine it with
  # numerical differentiation. This test is added to check breaking changes.
  estimated_variance = delta_method_log_mutinfo(
    fitted_model = fitted_model,
    rotation_par_unid = c(0,0,0,0),
    copula_par_unid = 1:4,
    copula_family2 = "clayton",
    n_prec = 1e3,
    composite = TRUE,
    seed = 1,
    eps = 1e-2
  )
  expect_equal(estimated_variance, matrix(0.309301300833, nrow = 1, ncol = 1))
})

test_that("Resampling method for the ICA based on fitted model summary information works", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # The computation of the ICA is way to inaccurate here to combine it with
  # numerical differentiation. This test is added to check breaking changes.
  set.seed(1)
  estimated_ICAs = summary_level_bootstrap_ICA(
    fitted_model = fitted_model,
    rotation_par_unid = c(0,0,0,0),
    copula_par_unid = 1:4,
    copula_family2 = "clayton",
    n_prec = 1e3,
    composite = TRUE,
    seed = 1,
    B = 10L
  )
  set.seed(1)
  estimated_ICAs_p = summary_level_bootstrap_ICA(
    fitted_model = fitted_model,
    rotation_par_unid = c(0,0,0,0),
    copula_par_unid = 1:4,
    copula_family2 = "clayton",
    n_prec = 1e3,
    composite = TRUE,
    seed = 1,
    B = 10L,
    ncores = 2
  )
  mean_ICA = mean(estimated_ICAs); sd_ICA = sd(estimated_ICAs)
  mean_ICA_p = mean(estimated_ICAs_p); sd_ICA_p = sd(estimated_ICAs_p)
  expect_equal(c(mean_ICA, sd_ICA), c(0.99096747638744, 0.00145160206983))
  expect_equal(c(mean_ICA_p, sd_ICA_p), c(0.99096747638744, 0.00145160206983))
})
