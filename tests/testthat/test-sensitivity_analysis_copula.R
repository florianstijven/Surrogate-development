test_that("basic sensitivity analysis works for continuous-continuous, ordinal-continuous, and ordinal-ordinal", {
  # Load fitted copula model.
  fitted_model_contcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-ContCont.rds"))
  fitted_model_ordcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-OrdCont.rds"))
  fitted_model_ordord = readRDS(test_path("fixtures", "schizo-dvine-gaussian-OrdOrd.rds"))

  # Perform sensitivity analysis
  set.seed(1)
  results_contcont = sensitivity_analysis_copula(
    fitted_model = fitted_model_contcont,
    n_sim = 10,
    eq_cond_association = FALSE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = TRUE,
    copula_family2 = "clayton",
    n_prec = 2e3,
    ncores = 1
  )
  set.seed(1)
  results_ordcont = sensitivity_analysis_copula(
    fitted_model = fitted_model_ordcont,
    n_sim = 10,
    eq_cond_association = TRUE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = FALSE,
    copula_family2 = "frank",
    n_prec = 2e3,
    ncores = 1
  )
  set.seed(1)
  results_ordord = sensitivity_analysis_copula(
    fitted_model = fitted_model_ordord,
    n_sim = 10,
    eq_cond_association = FALSE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = TRUE,
    copula_family2 = "gumbel",
    n_prec = 2e3,
    ncores = 1
  )
  # Correct ICA
  {
  expect_equal(
    mean(results_contcont[, 1]),
    c(0.6710603),
    tolerance = 0.15
  )
  # expect_equal(
  #   as.numeric(results_contcont[1, ]),
  #   c(0.805547826084, 0.864972741243, 0.909834719459, 0.591356631839, 0.691666788417, 0.361532829383, 0.450635314659,
  #     0.923559431890, 6.105601327913, 0.649720879856, 6.596600115926, 2.467916678980, 0.563217212642, 1.070967086009,
  #     0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000
  #   ),
  #   ignore_attr = "names"
  # )

  expect_equal(
    mean(results_ordcont[, 1]),
    c(0.07629387),
    tolerance = 0.1
  )
  # expect_equal(
  #   as.numeric(results_ordcont[1, ]),
  #   c(0.141644267022, 0.435302162260,             NA,             NA,             NA,             NA,             NA,
  #     NA, 3.142725415951, 2.303975808732, 4.016510351646, 2.047511165423, 2.047511165423, 3.433422135475,
  #     0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000
  #   ),
  #   ignore_attr = "names"
  # )

  expect_equal(
    mean(results_ordord[, 1]),
    c(0.2041518),
    tolerance = 0.15
  )
  # expect_equal(
  #   as.numeric(results_contcont[1, ]),
  #   c(0.805547826084, 0.864972741243, 0.909834719459, 0.591356631839, 0.691666788417, 0.361532829383, 0.450635314659,
  #     0.923559431890, 6.105601327913, 0.649720879856, 6.596600115926, 2.467916678980, 0.563217212642, 1.070967086009,
  #     0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000
  #   ),
  #   ignore_attr = "names"
  # )
  }
})

test_that("sensitivity analysis works for continuous-continuous, ordinal-continuous, and ordinal-ordinal with
          default user-specified definition of the ICA", {
  # Load fitted copula model.
  fitted_model_contcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-ContCont.rds"))
  fitted_model_ordcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-OrdCont.rds"))
  fitted_model_ordord = readRDS(test_path("fixtures", "schizo-dvine-gaussian-OrdOrd.rds"))


  # Perform sensitivity analysis without a user-specified definition of the ICA.
  {
  set.seed(1)
  results_contcont = sensitivity_analysis_copula(
    fitted_model = fitted_model_contcont,
    n_sim = 5,
    eq_cond_association = FALSE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = TRUE,
    copula_family2 = "clayton",
    n_prec = 2e3,
    ncores = 1
  )
  set.seed(1)
  results_ordcont = sensitivity_analysis_copula(
    fitted_model = fitted_model_ordcont,
    n_sim = 5,
    eq_cond_association = TRUE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = FALSE,
    copula_family2 = "frank",
    n_prec = 2e3,
    ncores = 1
  )
  set.seed(1)
  results_ordord = sensitivity_analysis_copula(
    fitted_model = fitted_model_ordord,
    n_sim = 5,
    eq_cond_association = FALSE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = TRUE,
    copula_family2 = "gumbel",
    n_prec = 2e3,
    ncores = 1
  )
  }
  # Perform sensitivity analysis with user-specified definition of the ICA that correspond to the defaults.
  {
  set.seed(1)
  results_contcont_default = sensitivity_analysis_copula(
    fitted_model = fitted_model_contcont,
    n_sim = 5,
    eq_cond_association = FALSE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = TRUE,
    copula_family2 = "clayton",
    n_prec = 2e3,
    ncores = 1,
    ICA_estimator = constructor_ICA_estimator(
      c("continuous", "continuous"),
      function(mutinfo, H_DeltaS, H_DeltaT) {
        1 - exp(-2 * mutinfo)
      }
      )
  )
  set.seed(1)
  results_ordcont_default = sensitivity_analysis_copula(
    fitted_model = fitted_model_ordcont,
    n_sim = 5,
    eq_cond_association = TRUE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = FALSE,
    copula_family2 = "frank",
    n_prec = 2e3,
    ncores = 1,
    ICA_estimator = constructor_ICA_estimator(
      c("ordinal", "continuous"),
      function(mutinfo, H_DeltaS, H_DeltaT) {
        mutinfo / H_DeltaT
      }
    )
  )
  set.seed(1)
  results_ordord_default = sensitivity_analysis_copula(
    fitted_model = fitted_model_ordord,
    n_sim = 5,
    eq_cond_association = FALSE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = TRUE,
    copula_family2 = "gumbel",
    n_prec = 2e3,
    ncores = 1,
    ICA_estimator = constructor_ICA_estimator(
      c("ordinal", "ordinal"),
      function(mutinfo, H_DeltaS, H_DeltaT) {
        mutinfo / min(H_DeltaS, H_DeltaT)
      }
    )
  )
  }

  # The results of the sensitivity analyses above should be identitical.
  {
    expect_equal(
      results_contcont,
      results_contcont_default
    )
    expect_equal(
      results_ordcont,
      results_ordcont_default
    )
    expect_equal(
      results_ordord,
      results_ordord_default
    )
  }


})

test_that("sensitivity analysis works for continuous-continuous, ordinal-continuous, and ordinal-ordinal with
          default user-specified definition of the ICA", {
  # Load fitted copula model.
  fitted_model_contcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-ContCont.rds"))
  fitted_model_ordcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-OrdCont.rds"))
  fitted_model_ordord = readRDS(test_path("fixtures", "schizo-dvine-gaussian-OrdOrd.rds"))


  # Perform sensitivity analysis with a non-default user-specified definition of
  # the ICA.
  {
    set.seed(1)
    results_contcont = sensitivity_analysis_copula(
      fitted_model = fitted_model_contcont,
      n_sim = 5,
      eq_cond_association = FALSE,
      lower = rep(0.2, 4),
      upper = rep(0.8, 4),
      degrees = 0,
      marg_association = TRUE,
      copula_family2 = "clayton",
      n_prec = 2e3,
      ncores = 1,
      ICA_estimator = constructor_ICA_estimator(
        c("continuous", "continuous"),
        function(mutinfo, H_DeltaS, H_DeltaT) {
          sqrt(1 - exp(-2 * mutinfo))
        }
      )
    )
    set.seed(1)
    results_ordcont = sensitivity_analysis_copula(
      fitted_model = fitted_model_ordcont,
      n_sim = 5,
      eq_cond_association = TRUE,
      lower = rep(0.2, 4),
      upper = rep(0.8, 4),
      degrees = 0,
      marg_association = FALSE,
      copula_family2 = "frank",
      n_prec = 2e3,
      ncores = 1,
      ICA_estimator = constructor_ICA_estimator(
        c("ordinal", "continuous"),
        function(mutinfo, H_DeltaS, H_DeltaT) {
          (1 - exp(-2 * mutinfo)) / (1 - exp(-2 * H_DeltaT))
        }
      )
    )
    set.seed(1)
    results_ordord = sensitivity_analysis_copula(
      fitted_model = fitted_model_ordord,
      n_sim = 5,
      eq_cond_association = FALSE,
      lower = rep(0.2, 4),
      upper = rep(0.8, 4),
      degrees = 0,
      marg_association = TRUE,
      copula_family2 = "gumbel",
      n_prec = 2e3,
      ncores = 1,
      ICA_estimator = constructor_ICA_estimator(
        c("ordinal", "ordinal"),
        function(mutinfo, H_DeltaS, H_DeltaT) {
          (1 - exp(-2 * mutinfo)) / (1 - exp(-2 * min(H_DeltaS, H_DeltaT)))
        }
      )
    )
  }

  # Check the computed ICAs in the sensitivity analyses.
  {
    expect_equal(
      mean(results_contcont[, 1]),
      c(0.8484771),
      tolerance = 0.2
    )
    expect_equal(
      mean(results_ordcont[, 1]),
      c(0.161265),
      tolerance = 0.15
    )
    expect_equal(
      mean(results_ordord[, 1]),
      c(0.3747666),
      tolerance = 0.2
    )
  }


})
