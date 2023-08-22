test_that("Clayton copula likelihood works with right-censoring", {
  u = c(0.2, 0.5, 0.3, 0.25, 0.98)
  v = c(0.54, 0.25, 0.01, 0.99, 0.5)
  d1 = c(0, 1, 0, 1, 0)
  d2 = c(0, 0, 1, 1, 0)
  theta = 0.001
  log_lik = clayton_loglik_copula_scale(theta, u, v, d1, d2)
  expect_equal(log_lik, -6.2503256)
})

test_that("Frank copula likelihood works with right-censoring", {
  u = c(0.2, 0.5, 0.3, 0.25, 0.98)
  v = c(0.54, 0.25, 0.01, 0.99, 0.5)
  d1 = c(0, 1, 0, 1, 0)
  d2 = c(0, 0, 1, 1, 0)
  theta = 0.001
  log_lik = frank_loglik_copula_scale(theta, u, v, d1, d2)
  expect_equal(log_lik, -6.2492926)
})

test_that("Gumbel copula likelihood works with right-censoring", {
  u = c(0.2, 0.5, 0.3, 0.25, 0.98)
  v = c(0.54, 0.25, 0.01, 0.99, 0.5)
  d1 = c(0, 1, 0, 1, 0)
  d2 = c(0, 0, 1, 1, 0)
  theta = 1
  log_lik = gumbel_loglik_copula_scale(theta, u, v, d1, d2)
  expect_equal(log_lik, -6.2491995)
})

test_that("Gaussian copula likelihood works with right-censoring", {
  u = c(0.2, 0.5, 0.3, 0.25, 0.98)
  v = c(0.54, 0.25, 0.01, 0.99, 0.5)
  d1 = c(0, 1, 0, 1, 0)
  d2 = c(0, 0, 1, 1, 0)
  theta = 0
  log_lik = gaussian_loglik_copula_scale(theta, u, v, d1, d2)
  expect_equal(log_lik, -6.2491995)
})

test_that("Clayton copula likelihood works with left-censoring of second variable", {
  u = c(0.2, 0.5, 0.3, 0.25, 0.98)
  v = c(0.54, 0.25, 0.01, 0.99, 0.5)
  d1 = c(0, 1, 0, 1, 0)
  d2 = c(0, 0, 1, -1, 0)
  theta = 0.001
  log_lik = clayton_loglik_copula_scale(theta, u, v, d1, d2)
  expect_equal(log_lik, -6.2599892)
})

test_that("Frank copula likelihood works with left-censoring of second variable", {
  u = c(0.2, 0.5, 0.3, 0.25, 0.98)
  v = c(0.54, 0.25, 0.01, 0.99, 0.5)
  d1 = c(0, 1, 0, 1, 0)
  d2 = c(0, 0, 1, -1, 0)
  theta = 0.001
  log_lik = frank_loglik_copula_scale(theta, u, v, d1, d2)
  expect_equal(log_lik, -6.2590954)
})

test_that("Gumbel copula likelihood works with left-censoring of second variable", {
  u = c(0.2, 0.5, 0.3, 0.25, 0.98)
  v = c(0.54, 0.25, 0.01, 0.99, 0.5)
  d1 = c(0, 1, 0, 1, 0)
  d2 = c(0, 0, 1, -1, 0)
  theta = 1
  log_lik = gumbel_loglik_copula_scale(theta, u, v, d1, d2)
  expect_equal(log_lik, -6.2592499)
})

test_that("Gaussian copula likelihood works with left-censoring of second variable", {
  u = c(0.2, 0.5, 0.3, 0.25, 0.98)
  v = c(0.54, 0.25, 0.01, 0.99, 0.5)
  d1 = c(0, 1, 0, 1, 0)
  d2 = c(0, 0, 1, -1, 0)
  theta = 0
  log_lik = gaussian_loglik_copula_scale(theta, u, v, d1, d2)
  expect_equal(log_lik, -6.2592499)
})

