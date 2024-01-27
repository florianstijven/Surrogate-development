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
  expect_equal(estimated_variance, matrix(0.158789300236, nrow = 1, ncol = 1))
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
  mean_ICA = mean(estimated_ICAs)
  sd_ICA = sd(estimated_ICAs)
  expect_equal(c(mean_ICA, sd_ICA), c(0.99124170985157, 0.00166932399608))
})
