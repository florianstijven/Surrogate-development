test_that("the ordinal-ordinal loglikelihood works", {
  para = c(-2, 0:3, 0:4, 2)
  X = c(1, 2, 2, 5, 3, 6, 4)
  Y = c(2, 1, 2, 5, 6, 3, 4)
  K = 6
  expect_equal(
    ordinal_ordinal_loglik(para, X, Y, "clayton", K, K),
    -47.6053539307
  )
  expect_equal(
    ordinal_ordinal_loglik(para, X, Y, "frank", K, K),
    -42.5758336041
  )
  # Correlation cannot be 2.
  expect_equal(
    ordinal_ordinal_loglik(c(para[-1], 0.5), X, Y, "gaussian", K, K),
    -28.1181684056
  )
  expect_equal(
    ordinal_ordinal_loglik(para, X, Y, "gumbel", K, K),
    -52.5339288404
  )
})

test_that("fit_copula_submodel_OrdOrd() works for the twostep estimator", {
  X = c(1, 2, 2, 5, 3, 6, 4)
  Y = c(2, 1, 2, 5, 3, 5, 4)
  K_X = 6
  K_Y = 5
  fitted_submodel = fit_copula_submodel_OrdOrd(
    X = X,
    Y = Y,
    copula_family = "clayton",
    start_copula = 2,
    K_X = K_X,
    K_Y = K_Y,
    twostep = TRUE
  )
  expect_equal(
    logLik(fitted_submodel$ml_fit),
    -20.1225445332,
    ignore_attr = "df"
  )
  expect_equal(
    fitted_submodel$marginal_X$pmf(1:5),
    c(0.142857142857, 0.285714285714, 0.142857142857, 0.142857142857, 0.142857142857),
    ignore_attr = "names"
  )
  expect_equal(
    fitted_submodel$marginal_Y$inv_cdf((1:5) / 5),
    c(2, 2, 4, 5, 5)
  )
})

test_that("fit_copula_submodel_OrdOrd() works for the full estimator", {
  X = c(1, 2, 2, 5, 3, 6, 4)
  X = rep(X, 10)
  Y = c(2, 1, 2, 5, 3, 5, 4)
  Y = rep(Y, 10)
  K_X = 6
  K_Y = 5

  fitted_submodel = fit_copula_submodel_OrdOrd(
    X = X,
    Y = Y,
    copula_family = "clayton",
    K_X = K_X,
    K_Y = K_Y,
    start_copula = 2
  )
  expect_equal(
    logLik(fitted_submodel$ml_fit),
    -193.005983363,
    ignore_attr = "df"
  )
  expect_equal(
    fitted_submodel$marginal_X$pmf(1:5),
    c(0.163139303754, 0.137925386782, 0.165417272812, 0.200230588170, 0.185809273741),
    ignore_attr = "names"
  )
  expect_equal(
    fitted_submodel$marginal_Y$inv_cdf((1:5) / 5),
    c(2, 3, 4, 5, 5)
  )
})


test_that("fit_copula_OrdOrd() works for the full estimator", {
  S0 = c(1, 2, 2, 5, 3, 6, 4)
  S0 = rep(S0, 10)
  S1 = c(1, 4, 2, 5, 2, 6, 4)
  S1 = rep(S1, 10)
  K_S = 6

  T0 = c(2, 1, 2, 5, 3, 5, 4)
  T0 = rep(T0, 10)
  T1 = c(1, 3, 2, 4, 2, 5, 4)
  T1 = rep(T1, 10)
  K_T = 5

  data = data.frame(
    surrogate = c(S0, S1),
    true = c(T0, T1),
    treat = c(
      rep(0, 70),
      rep(1, 70)
    )
  )

  fitted_model = fit_copula_OrdOrd(
    data = data,
    copula_family = "clayton",
    K_S = K_S,
    K_T = K_T,
    start_copula = 0.5
  )
  expect_equal(
    fitted_model$fit_0$ml_fit$maximum,
    -193.005983363
  )
})
