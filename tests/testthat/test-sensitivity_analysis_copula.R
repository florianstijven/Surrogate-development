test_that("multiplication works", {
  # Fit model
  S0 = c(2.2, 3.1, 0, -3, 0, 1, 4)
  S0 = rep(S0, 10)
  S1 = c(2, -0.5, 0, -3, 0, 2, 0)
  S1 = rep(S1, 10)

  T0 = c(1, 6, 2, 5, 3, 6, 4)
  T0 = rep(T0, 10)
  T1 = c(1, 5, 2, 6, 3, 6, 4)
  T1 = rep(T1, 10)
  K_T = 6

  data = data.frame(
    surrogate = c(S0, S1),
    true = c(T0, T1),
    treat = c(
      rep(0, 70),
      rep(1, 70)
    )
  )

  marginal = list(
    pdf_fun = function(x, para) {
      dnorm(x, mean = para[1], sd = para[2])
    },
    cdf_fun = function(x, para) {
      pnorm(x, mean = para[1], sd = para[2])
    },
    q_fun = function(p, para) {
      qnorm(p, mean = para[1], sd = para[2])
    },
    n_para = 2
  )

  fitted_model = fit_copula_OrdCont(
    data = data,
    copula_family = "frank",
    marginal_S0 = marginal,
    marginal_S1 = marginal,
    K_T = K_T,
    start_copula = 2
  )

  # Perform sensitivity analysis
  set.seed(1)
  results = sensitivity_analysis_copula(
    fitted_model = fitted_model,
    n_sim = 5,
    eq_cond_association = FALSE,
    lower = rep(0.2, 4),
    upper = rep(0.8, 4),
    degrees = 0,
    marg_association = TRUE,
    copula_family2 = "clayton",
    n_prec = 5e3,
    ncores = 1
  )
  # Correct ICA
  expect_equal(
    results[, 1],
    c(
      0.0825988143900,
      0.0821616252617,
      0.1136904411633,
      0.0488619958643,
      0.1155981371652
    )
  )
  # Correct Spearman's rho
  expect_equal(
    results[, 2],
    c(
      -0.547959472635,
      -0.540196293513,
      -0.608194702531,
      -0.395349607610,
      -0.614067609887
    )
  )
})
