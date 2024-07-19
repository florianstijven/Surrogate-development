test_that("the ordinal-continuous loglikelihood works", {
  para = c(1:2, 1:2, 2)
  X = c(1, 2, 2, 5, 3, 6, 4)
  Y = c(2.2, 3.1, 0, -3, 0, 1, 4)
  copula_family = "clayton"
  marginal_X = list(
    pdf_fun = function(x, para) {
      dweibull(x, shape = para[1], scale = para[2])
    },
    cdf_fun = function(x, para) {
      pweibull(x, shape = para[1], scale = para[2])
    },
    NA,
    n_para = 2
  )
  marginal_Y = list(
    pdf_fun = function(x, para) {
      dnorm(x, mean = para[1], sd = para[2])
    },
    cdf_fun = function(x, para) {
      pnorm(x, mean = para[1], sd = para[2])
    },
    NA,
    n_para = 2
  )
  K = 6
  expect_equal(
    continuous_continuous_loglik(para, X, Y, "clayton", marginal_X, marginal_Y),
    -38.082448704
  )
  expect_equal(
    continuous_continuous_loglik(c(para[-1], 0.5), X, Y, "gaussian", marginal_X, marginal_Y),
    -139.6712877
  )
  expect_equal(
    continuous_continuous_loglik(para, X, Y, "frank", marginal_X, marginal_Y),
    -32.57395345
  )
  expect_equal(
    continuous_continuous_loglik(para, X, Y, "gumbel", marginal_X, marginal_Y),
    -36.81896141
  )
})

test_that("fit_copula_submodel_ContCont() works for the twostep estimator", {
  X = c(1, 2, 2, 5, 3, 6, 4)
  Y = c(2.2, 3.1, 0, -3, 0, 1, 4)
  marginal_X = list(
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
  marginal_Y = marginal_X
  fitted_submodel = fit_copula_submodel_ContCont(
    X = X,
    Y = Y,
    copula_family = "clayton",
    marginal_X = marginal_X,
    marginal_Y = marginal_Y,
    start_X = c(mean(X), sd(X)),
    start_Y = c(mean(Y), sd(Y)),
    start_copula = 2,
    twostep = TRUE
  )
  expect_equal(
    logLik(fitted_submodel$ml_fit),
    -25.3351941499,
    ignore_attr = "df"
  )
  expect_equal(
    fitted_submodel$marginal_X$pdf(1:5),
    c(0.0934309799828, 0.1777904010307, 0.2359673786127, 0.2184348251465, 0.1410321230219),
    ignore_attr = "names"
  )
  expect_equal(
    fitted_submodel$marginal_Y$inv_cdf((1:5) / 5),
    c(-0.779990392239, 0.494138641766,  1.591575646335, 2.865704680340, Inf)
  )
})

test_that("fit_copula_submodel_OrdCont() works for the full estimator", {
  X = c(1, 6, 2, 5, 3, 6, 4)
  X = rep(X, 10)
  Y = c(2.2, 3.1, 0, -3, 0, 1, 4)
  Y = rep(Y, 10)

  marginal_X = list(
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
  marginal_Y = marginal_X

  fitted_submodel = fit_copula_submodel_ContCont(
    X = X,
    Y = Y,
    copula_family = "clayton",
    marginal_X = marginal_X,
    marginal_Y = marginal_Y,
    start_X = c(mean(X), sd(X)),
    start_Y = c(mean(Y), sd(Y)),
    start_copula = 2
  )
  expect_equal(
    logLik(fitted_submodel$ml_fit),
    -275.337419571,
    ignore_attr = "df"
  )
  expect_equal(
    fitted_submodel$marginal_X$pdf(1:5),
    c(0.0653537833828, 0.1218625519486, 0.1768479630656, 0.1997377321968, 0.1755700474436),
    ignore_attr = "names"
  )
  expect_equal(
    fitted_submodel$marginal_Y$inv_cdf((1:5) / 5),
    c(-0.963697255915, 0.178130044707,  1.161612485722, 2.303439786344, Inf)
  )
})
