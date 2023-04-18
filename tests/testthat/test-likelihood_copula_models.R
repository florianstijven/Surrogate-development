test_that(
  "loglikelihood function for bivariate copula model works for clayton copula and normal margins",
  {
    copula = "clayton"
    marginal = "normal"
    set.seed(1)
    x = rnorm(5)
    y = rnorm(5)
    d1 = c(0, 1, 0, 1, 0)
    d2 = c(0, 0, 1, 1, 0)
    theta = 0.5
    loglik = log_likelihood_copula_model(
      theta = theta,
      X = x,
      Y = y,
      d1 = d1,
      d2 = d2,
      copula = copula,
      cdf_X = cdf_fun(1, 1, "normal"),
      cdf_Y = cdf_fun(1, 1, "normal"),
      pdf_X = pdf_fun(1, 1, "normal"),
      pdf_Y = pdf_fun(1, 1, "normal")
    )
    expect_equal(loglik, -16.551776, tolerance = 1e-6)
  }
)

test_that(
  "loglikelihood function for bivariate copula model works for gaussian copula and normal margins",
  {
    copula = "gaussian"
    marginal = "normal"
    set.seed(1)
    x = rnorm(5)
    y = rnorm(5)
    d1 = c(0, 1, 0, 1, 0)
    d2 = c(0, 0, 1, 1, 0)
    theta = 0.5
    loglik = log_likelihood_copula_model(
      theta = theta,
      X = x,
      Y = y,
      d1 = d1,
      d2 = d2,
      copula = copula,
      cdf_X = cdf_fun(1, 1, "normal"),
      cdf_Y = cdf_fun(1, 1, "normal"),
      pdf_X = pdf_fun(1, 1, "normal"),
      pdf_Y = pdf_fun(1, 1, "normal")
    )
    expect_equal(loglik, -16.10984, tolerance = 1e-6)
  }
)

test_that(
  "loglikelihood function for bivariate copula model works for gaussian copula and logistic margins",
  {
    copula = "gaussian"
    marginal = "logistic"
    set.seed(1)
    x = rnorm(5)
    y = rnorm(5)
    d1 = c(0, 1, 0, 1, 0)
    d2 = c(0, 0, 1, 1, 0)
    theta = 0.5
    loglik = log_likelihood_copula_model(
      theta = theta,
      X = x,
      Y = y,
      d1 = d1,
      d2 = d2,
      copula = copula,
      cdf_X = cdf_fun(1, 1, marginal),
      cdf_Y = cdf_fun(1, 1, marginal),
      pdf_X = pdf_fun(1, 1, marginal),
      pdf_Y = pdf_fun(1, 1, marginal)
    )
    expect_equal(loglik, -13.321792)
  }
)
