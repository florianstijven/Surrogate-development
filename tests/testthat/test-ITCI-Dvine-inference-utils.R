# test_that("Delta method based on numerical derivatives works", {
#   # Load fitted copula model.
#   fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
#   # The computation of the ICA is way to inaccurate here to combine it with
#   # numerical differentiation. This test is added to check breaking changes.
#   estimated_variance = delta_method_log_mutinfo(
#     fitted_model = fitted_model,
#     rotation_par_unid = c(0,0,0,0),
#     copula_par_unid = 1:4,
#     copula_family2 = "clayton",
#     n_prec = 1e3,
#     composite = TRUE,
#     seed = 1,
#     eps = 1e-2
#   )
#   expect_equal(estimated_variance, matrix(0.309301300833, nrow = 1, ncol = 1))
# })

test_that("Resampling method for the ICA based on fitted model summary information works for survival-survival", {
  # Load fitted copula model.
  fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
  # The computation of the ICA is way to inaccurate here to combine it with
  # numerical differentiation. This test is added to check breaking changes.
  set.seed(1)
  suppressWarnings({
    estimated_ICAs = summary_level_bootstrap_ICA(
      fitted_model = fitted_model,
      rotation_par_unid = c(0, 0, 0, 0),
      copula_par_unid = 1:4,
      copula_family2 = "clayton",
      n_prec = 1e3,
      composite = TRUE,
      seed = 1,
      B = 10L
    )
  })

  set.seed(1)
  suppressWarnings({
    estimated_ICAs_p = summary_level_bootstrap_ICA(
      fitted_model = fitted_model,
      rotation_par_unid = c(0, 0, 0, 0),
      copula_par_unid = 1:4,
      copula_family2 = "clayton",
      n_prec = 1e3,
      composite = TRUE,
      seed = 1,
      B = 10L,
      ncores = 2
    )
  })

  mean_ICA = mean(estimated_ICAs); sd_ICA = sd(estimated_ICAs)
  mean_ICA_p = mean(estimated_ICAs_p); sd_ICA_p = sd(estimated_ICAs_p)
  expect_equal(c(mean_ICA, sd_ICA), c(0.991289115378347, 0.000645362133079), tolerance = 0.01)
  expect_equal(c(mean_ICA_p, sd_ICA_p), c(0.991289115378347, 0.000645362133079), tolerance = 0.01)
})

test_that("Resampling method for the ICA based on fitted model summary information works for ordinal/continuous settings", {
  # Load fitted copula model.
  fitted_model_contcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-ContCont.rds"))
  fitted_model_ordcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-OrdCont.rds"))
  fitted_model_ordord = readRDS(test_path("fixtures", "schizo-dvine-gaussian-OrdOrd.rds"))
  # Compute the ICA in bootstrap resamples.
  estimated_ICAs_contcont = summary_level_bootstrap_ICA(
    fitted_model = fitted_model_contcont,
    rotation_par_unid = c(0,0,0,0),
    copula_par_unid = 1:4,
    copula_family2 = "clayton",
    n_prec = 1e3,
    seed = 1,
    B = 10L
  )
  estimated_ICAs_ordcont = summary_level_bootstrap_ICA(
    fitted_model = fitted_model_ordcont,
    rotation_par_unid = c(0,0,0,0),
    copula_par_unid = 1:4,
    copula_family2 = "clayton",
    n_prec = 5e3,
    seed = 1,
    B = 10L
  )
  estimated_ICAs_ordord = summary_level_bootstrap_ICA(
    fitted_model = fitted_model_ordord,
    rotation_par_unid = c(0,0,0,0),
    copula_par_unid = 1:4,
    copula_family2 = "clayton",
    n_prec = 1e3,
    seed = 1,
    B = 10L
  )
  # ICAs correct for contcont
  expect_equal(
    mean(estimated_ICAs_contcont),
    c(0.8742195),
    tolerance = 0.05
  )
  # ICAs correct for ordcont
  expect_equal(
    mean(estimated_ICAs_ordcont),
    c(0.03083768),
    tolerance = 0.05
  )
  # ICAs correct for ordord
  expect_equal(
    mean(estimated_ICAs_ordord),
    c(0.2507742),
    tolerance = 0.1
  )
})
