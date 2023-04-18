test_that(
  "loglikelihood function for bivariate copula model works for clayton copula and normal margins",
  {
    copula = "clayton"
    marginal = "normal"
    x = c(-0.6264538, 0.1836433, -0.8356286, 1.595280, 0.3295078)
    y = c(-0.8204684, 0.4874291, 0.7383247, 0.5757814, -0.3053884)
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
    x = c(-0.6264538, 0.1836433, -0.8356286, 1.595280, 0.3295078)
    y = c(-0.8204684, 0.4874291, 0.7383247, 0.5757814, -0.3053884)
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
    x = c(-0.6264538, 0.1836433, -0.8356286, 1.595280, 0.3295078)
    y = c(-0.8204684, 0.4874291, 0.7383247, 0.5757814, -0.3053884)
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
    expect_equal(loglik, -13.3217916)
  }
)

test_that(
  "loglikelihood function for bivariate copula model works for gaussian copula and logistic margins, with left-censoring",
  {
    copula = "gaussian"
    marginal = "logistic"
    x = c(-0.6264538, 0.1836433, -0.8356286, 1.595280, 0.3295078)
    y = c(-0.8204684, 0.4874291, 0.7383247, 0.5757814, -0.3053884)
    d1 = c(0, 1, 0, 1, 0)
    d2 = c(0, 0, 1, -1, 0)
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
    expect_equal(loglik, -11.122886, tolerance = 1e-5)
  }
)
