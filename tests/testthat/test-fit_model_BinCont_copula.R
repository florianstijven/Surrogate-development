test_that("binary_continuous_loglik() works with clayton copula and logistic margins", {
  copula = "clayton"
  marginal = "logistic"
  x = c(-0.6264538, 0.1836433,-0.8356286, 1.595280, 0.3295078)
  y = c(0, 1, 0, 1, 0)
  loglik = binary_continuous_loglik(
    para = c(0, 0, 2, 1.2),
    X = x,
    Y = y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  expect_equal(loglik, -5.75446910)
})

test_that("binary_continuous_loglik() works with clayton copula and logistic margins", {
  copula = "gumbel"
  marginal = "logistic"
  x = c(-0.6264538, 0.1836433,-0.8356286, 1.595280, 0.3295078)
  y = c(0, 1, 0, 1, 0)
  loglik = binary_continuous_loglik(
    para = c(0, 0, 2, 1.2),
    X = x,
    Y = y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  expect_equal(loglik, -5.5347204)
})

test_that("binary_continuous_loglik() works with clayton copula and logistic margins", {
  copula = "gaussian"
  marginal = "logistic"
  x = c(-0.6264538, 0.1836433,-0.8356286, 1.595280, 0.3295078)
  y = c(0, 1, 0, 1, 0)
  loglik = binary_continuous_loglik(
    para = c(0, 0, 2, 0.75),
    X = x,
    Y = y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  expect_equal(loglik, -5.4788387)
})

test_that("binary_continuous_loglik() works with clayton copula and logistic margins", {
  copula = "frank"
  marginal = "logistic"
  x = c(-0.6264538, 0.1836433,-0.8356286, 1.595280, 0.3295078)
  y = c(0, 1, 0, 1, 0)
  loglik = binary_continuous_loglik(
    para = c(0, 0, 2, 0.75),
    X = x,
    Y = y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  expect_equal(loglik, -5.4788387)
})

test_that("fit_copula_submodel_BinCont() works with clayton copula and normal margins", {
  copula_family = "gaussian"
  marginal_surrogate = "normal"
  data("Ovarian")
  X = -1 * Ovarian$Pfs[1:200]
  Y = Ovarian$SurvInd[1:200]
  BinCont_starting_values(
    X = X,
    Y = Y,
    copula_family = copula_family,
    marginal_surrogate = marginal_surrogate
  )
  binary_continuous_loglik(
    para = c(0, 0, 2, 0.7),
    X = X,
    Y = Y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  ml_fit = fit_copula_submodel_BinCont(
    X = x,
    Y = y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  summary(ml_fit)
  ml_fit$hessian
  expect_equal(loglik, -5.75446910)
}
)

test_that("Simulate data", {
  data = copula::rCopula(copula = copula::claytonCopula(param = 3), n = 500)
  data = qnorm(data)
  X = data[, 1]
  Y = ifelse(data[, 2] < 0, 0, 1)
  copula_family = "clayton"
  marginal_surrogate = "normal"
  start = BinCont_starting_values(
    X = X,
    Y = Y,
    copula_family = copula_family,
    marginal_surrogate = marginal_surrogate
  )
  binary_continuous_loglik(
    para = start,
    X = X,
    Y = Y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  binary_continuous_loglik(
    para = c(0, 0, 1, 3),
    X = X,
    Y = Y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  ml_fit = fit_copula_submodel_BinCont(
    X = X,
    Y = Y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  summary(ml_fit)
  ml_fit$hessian
  expect_equal(loglik, -5.75446910)
}
)
