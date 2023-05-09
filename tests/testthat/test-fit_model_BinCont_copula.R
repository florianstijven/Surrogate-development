test_that("binary_continuous_loglik() works with clayton copula and logistic margins", {
  copula = "clayton"
  marginal = "logistic"
  X = c(-0.6264538, 0.1836433,-0.8356286, 1.595280, 0.3295078)
  Y = c(0, 1, 0, 1, 0)
  loglik = binary_continuous_loglik(
    para = c(0, 0, 2, 1.2),
    X = X,
    Y = Y,
    copula_family = copula,
    marginal_surrogate = marginal
  )
  expect_equal(loglik, -13.8160346)
})

# test_that("binary_continuous_loglik() works with clayton copula and logistic margins", {
#   copula = "gumbel"
#   marginal = "logistic"
#   x = c(-0.6264538, 0.1836433,-0.8356286, 1.595280, 0.3295078)
#   y = c(0, 1, 0, 1, 0)
#   loglik = binary_continuous_loglik(
#     para = c(0, 0, 2, 1),
#     X = x,
#     Y = y,
#     copula_family = copula,
#     marginal_surrogate = marginal
#   )
#   expect_equal(loglik, -5.5347204)
# })

test_that("binary_continuous_loglik() works with gaussian copula and logistic margins", {
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
  expect_equal(loglik, -13.416364)
})

test_that("binary_continuous_loglik() works with frank copula and logistic margins", {
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
  expect_equal(loglik, -13.9664643)
})

test_that("twostep_BinCont() works with clayton copula and gamma margins", {
  copula_family = "clayton"
  marginal_surrogate = "gamma"
  data("Schizo_BinCont")
  na = is.na(Schizo_BinCont$CGI_Bin) | is.na(Schizo_BinCont$PANSS)
  X = 120 - abs(Schizo_BinCont$PANSS[!na])
  Y = Schizo_BinCont$CGI_Bin[!na]
  BinCont_starting_values(
    X = X,
    Y = Y,
    copula_family = copula_family,
    marginal_surrogate = marginal_surrogate
  )
  twostep_fit = twostep_BinCont(
    X = X,
    Y = Y,
    copula_family = copula_family,
    marginal_surrogate = marginal_surrogate
  )
  expect_equal(coef(summary(twostep_fit))[4], 0.9675053)
}
)

test_that("twostep_BinCont() works with gaussian copula and weibull margins", {
  copula_family = "gaussian"
  marginal_surrogate = "weibull"
  data("Schizo_BinCont")
  na = is.na(Schizo_BinCont$CGI_Bin) | is.na(Schizo_BinCont$PANSS)
  X = 120 - abs(Schizo_BinCont$PANSS[!na])
  Y = Schizo_BinCont$CGI_Bin[!na]
  BinCont_starting_values(
    X = X,
    Y = Y,
    copula_family = copula_family,
    marginal_surrogate = marginal_surrogate
  )
  twostep_fit = twostep_BinCont(
    X = X,
    Y = Y,
    copula_family = copula_family,
    marginal_surrogate = marginal_surrogate
  )
  expect_equal(coef(summary(twostep_fit))[4], 0.58564315)
}
)
