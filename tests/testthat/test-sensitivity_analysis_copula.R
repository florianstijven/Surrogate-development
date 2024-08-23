test_that("basic sensitivity analysis works for continuous-continuous, ordinal-continuous, and ordinal-ordinal", {
  # Load fitted copula model.
  fitted_model_contcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-ContCont.rds"))
  fitted_model_ordcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-OrdCont.rds"))
  fitted_model_ordord = readRDS(test_path("fixtures", "schizo-dvine-gaussian-OrdOrd.rds"))

  # Perform sensitivity analysis
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
  # Correct ICA
  {
  expect_equal(
    results_contcont[, 1],
    c(
      0.805547826084,
      0.804501092570,
      0.828990437005,
      0.362778037076,
      0.849093738971
    )
  )
  expect_equal(
    as.numeric(results_contcont[1, ]),
    c(0.805547826084, 0.864972741243, 0.909834719459, 0.591356631839, 0.691666788417, 0.361532829383, 0.450635314659,
      0.923559431890, 6.105601327913, 0.649720879856, 6.596600115926, 2.467916678980, 0.563217212642, 1.070967086009,
      0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000
    ),
    ignore_attr = "names"
  )

  expect_equal(
    results_ordcont[, 1],
    c(
      0.1416477223515,
      0.1371087549909,
      0.0353045132874,
      0.0120349874214,
      0.0843796516144
    )
  )
  expect_equal(
    as.numeric(results_ordcont[1, ]),
    c(0.141647722351, 0.435308056965,             NA,             NA,             NA,             NA,             NA,
      NA, 3.142968155547, 2.303975808732, 4.017478729159, 2.047511165423, 2.047511165423, 3.433422135475,
      0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000
    ),
    ignore_attr = "names"
  )
  expect_equal(
    results_ordord[, 1],
    c(
      0.192161887674,
      0.194310196872,
      0.221638256151,
      0.158499048460,
      0.198650867242
    )
  )
  expect_equal(
    as.numeric(results_contcont[1, ]),
    c(0.805547826084, 0.864972741243, 0.909834719459, 0.591356631839, 0.691666788417, 0.361532829383, 0.450635314659,
      0.923559431890, 6.105601327913, 0.649720879856, 6.596600115926, 2.467916678980, 0.563217212642, 1.070967086009,
      0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000
    ),
    ignore_attr = "names"
  )
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
      results_contcont[, 1],
      c(
        0.897523161865,
        0.896939848914,
        0.910489119652,
        0.602310581906,
        0.921462825605
      )
    )
    expect_equal(
      results_ordcont[, 1],
      c(
        0.2475845548787,
        0.2334619837987,
        0.0546455966377,
        0.0210433689901,
        0.1345990831075
      )
    )
    expect_equal(
      results_ordord[, 1],
      c(
        0.332139509278,
        0.330488272383,
        0.354290465704,
        0.242721743253,
        0.347879741871
      )
    )
  }





})
