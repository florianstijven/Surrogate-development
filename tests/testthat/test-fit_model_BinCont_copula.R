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
