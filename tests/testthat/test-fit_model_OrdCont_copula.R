test_that("ordinal to cutpoints conversion works", {
  cutpoints = -1:3
  x = c(1, 2, 2, 5, 3, 6, 4)
  expect_equal(
    ordinal_to_cutpoints(x, cutpoints, TRUE),
    c(-Inf, -1, -1, 2, 0, 3, 1)
  )
  expect_equal(
    ordinal_to_cutpoints(x, cutpoints, FALSE),
    c(-1, 0, 0, 3, 1, +Inf, 2)
  )
})

test_that("the ordinal-continuous loglikelihood works", {
  para = c(-1:3, 0, 2, 2)
  X = c(1, 2, 2, 5, 3, 6, 4)
  Y = c(2.2, 3.1, 0, -3, 0, 1, 4)
  copula_family = "clayton"
  marginal_Y = list(
    pdf_fun = function(x, para) {
      dnorm(x, mean = para[1], sd = para[2])
    },
    cdf_fun = function(x, para) {
      pnorm(x, mean = para[1], sd = para[2])
    },
    q_fun = NA,
    n_para = 2
  )
  K = 6
  expect_equal(
    ordinal_continuous_loglik(para, X, Y, copula_family, marginal_Y, K),
    -40.644492067
  )
})

test_that("the estimation of marginal distrbutions works", {
  Y = c(2.2, 3.1, 0, -3, 0, 1, 4)
  marginal_Y = list(
    pdf_fun = function(x, para) {
      dnorm(x, mean = para[1], sd = para[2])
    },
    cdf_fun = function(x, para) {
      pnorm(x, mean = para[1], sd = para[2])
    },
    q_fun = NA,
    n_para = 2
  )
  expect_equal(
    estimate_marginal(Y, marginal_Y, c(1, 1)),
    c(mean(Y), sd(Y)*sqrt(((length(Y) - 1) / length(Y)))),
    tolerance = 1e-4
  )
})

test_that("marginal_ord_constructor() works", {
  # WIP, OUTPUT IS STILL WRONG
  cut_points = 0:3
  marginal_distribution = marginal_ord_constructor(cut_points)
  # pmf is correct.
  expect_equal(
    marginal_distribution[[1]](1:5),
    c(0.500000000, 0.341344746, 0.135905122, 0.021400234, 0.001349898)
  )
  # cdf is correct.
  expect_equal(
    marginal_distribution[[2]](1:5),
    c(0.50000000, 0.84134475, 0.97724987, 0.99865010, 1.00000000)
  )
  # inverse cdf is correct.
  expect_equal(
    marginal_distribution[[3]]((1:5) / 5),
    c(1, 1, 2, 2, 5)
  )
})

test_that("fit_copula_submodel_OrdCont() works for the twostep estimator", {
  X = c(1, 2, 2, 5, 3, 6, 4)
  Y = c(2.2, 3.1, 0, -3, 0, 1, 4)
  K = 6
  marginal_Y = list(
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
  fitted_submodel = fit_copula_submodel_OrdCont(
    X = X,
    Y = Y,
    copula_family = "clayton",
    marginal_Y = marginal_Y,
    start_Y = c(mean(Y), sd(Y)),
    start_copula = 2,
    K = K,
    twostep = TRUE
  )
  expect_equal(
  logLik(fitted_submodel$ml_fit),
  -27.5774225982,
  ignore_attr = "df",
  tolerance = 1e-5
  )
  expect_equal(
    fitted_submodel$marginal_X$pmf(1:5),
    c(0.142857142857, 0.285714285714, 0.142857142857, 0.142857142857, 0.142857142857),
    ignore_attr = "names",
    tolerance = 1e-5
  )
  expect_equal(
    fitted_submodel$marginal_Y$inv_cdf((1:5) / 5),
    c(-0.779990392239, 0.494138641766,  1.591575646335, 2.865704680340, Inf),
    tolerance = 1e-5
  )
})

test_that("fit_copula_submodel_OrdCont() works for the full estimator", {
  X = c(1, 3, 2, 2, 3, 4, 4, 3, 3, 2)
  X = rep(X, 10)
  Y = c(2.2, 3.1, 0, -3, 0, 1, 4, 2.3, 0.5, 3.1)
  Y = rep(Y, 10)
  K = 4
  marginal_Y = list(
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
  fitted_submodel = fit_copula_submodel_OrdCont(
    X = X,
    Y = Y,
    copula_family = "clayton",
    marginal_Y = marginal_Y,
    start_Y = c(mean(Y), sd(Y)),
    start_copula = 3,
    K = K,
  )
  expect_equal(
    logLik(fitted_submodel$ml_fit),
    c(-335.9523),
    ignore_attr = "df",
    tolerance = 0.01
  )
  expect_equal(
    fitted_submodel$marginal_X$pmf(1:4),
    c(
      0.09757605,
      0.29212648,
      0.40887020,
      0.20142728
      ),
      ignore_attr = "names",
      tolerance = 1e-5
    )
    expect_equal(
      fitted_submodel$marginal_Y$inv_cdf((1:5) / 5),
      c(-0.3360741, 0.8170597, 1.8102807, 2.9634145, Inf),
      tolerance = 1e-5
    )
})

test_that("fit_copula_OrdCont() works for the full estimator", {
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
  expect_equal(
    fitted_model$fit_0$ml_fit$maximum,
    -275.637127024
  )
})

test_that("GoF functions work", {
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
  # Conditional mean function
  expect_equal(
    conditional_mean_copula_OrdCont(fitted_model$fit_0, grid = 1:5),
    c(3.86302464243, 3.92597904964, 3.97729279147, 4.01121752671, 4.02942336283),
    tolerance = 1
  )
}
)
