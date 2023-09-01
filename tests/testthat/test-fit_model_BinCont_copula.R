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
  expect_equal(coef(summary(twostep_fit$ml_fit))[4], 0.96837258)
}
)

test_that("twostep_BinCont() works with gaussian copula and weibull margins", {
  copula_family = "gaussian"
  marginal_surrogate = "lognormal"
  data("Schizo_BinCont")
  na = is.na(Schizo_BinCont$CGI_Bin) | is.na(Schizo_BinCont$PANSS)
  X = 6 - log(abs(Schizo_BinCont$PANSS[!na]) + 2)
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
  expect_equal(coef(summary(twostep_fit$ml_fit))[4], 0.63754401)
}
)

test_that("fit_copula_submodel_BinCont() works with clayton copula and lognormal margins", {
  copula_family = "clayton"
  marginal_surrogate = "weibull"
  data("Schizo_BinCont")
  na = is.na(Schizo_BinCont$CGI_Bin) | is.na(Schizo_BinCont$PANSS)
  X = 6 - log(abs(Schizo_BinCont$PANSS[!na]) + 2)
  Y = Schizo_BinCont$CGI_Bin[!na]
  fit = fit_copula_submodel_BinCont(
    X = X,
    Y = Y,
    copula_family = copula_family,
    marginal_surrogate = marginal_surrogate,
    method = "BFGS"
  )
  # Reference vector
  check_values = c(-792.949260, 3.328630, 3.235286)
  # Values from running the functions
  output_values = c(fit$ml_fit$maximum,
                    fit$marginal_S_dist$estimate)
  expect_equal(output_values, check_values, ignore_attr = TRUE)
}
)

test_that("fit_model_binCont_copula() works with clayton copula and lognormal margins and twostep estimator", {
  copula_family = "clayton"
  marginal_surrogate = "lognormal"
  data("Schizo_BinCont")
  na = is.na(Schizo_BinCont$CGI_Bin) | is.na(Schizo_BinCont$PANSS)
  X = 6 - log(abs(Schizo_BinCont$PANSS[!na]) + 2)
  Y = Schizo_BinCont$CGI_Bin[!na]
  Treat = Schizo_BinCont$Treat[!na]
  Treat = ifelse(Treat == 1, 1, 0)
  data = data.frame(X,
                    Y,
                    Treat)
  full_model = fit_copula_model_BinCont(data,
                                        copula_family,
                                        marginal_surrogate,
                                        twostep = TRUE,
                                        method = "BFGS")
  coef(full_model$submodel0$ml_fit)
  coefs_check = c(0.01958429, 1.08101496, 0.31669380, 1.86676923)
  expect_equal(coef(full_model$submodel0$ml_fit), coefs_check, ignore_attr = TRUE)
}
)

test_that("fit_model_binCont_copula() works with clayton copula and lognormal margins and full ML estimator", {
  copula_family = "clayton"
  marginal_surrogate = "lognormal"
  data("Schizo_BinCont")
  na = is.na(Schizo_BinCont$CGI_Bin) | is.na(Schizo_BinCont$PANSS)
  X = 6 - log(abs(Schizo_BinCont$PANSS[!na]) + 2)
  Y = Schizo_BinCont$CGI_Bin[!na]
  Treat = Schizo_BinCont$Treat[!na]
  Treat = ifelse(Treat == 1, 1, 0)
  data = data.frame(X,
                    Y,
                    Treat)
  full_model = fit_copula_model_BinCont(data,
                                        copula_family,
                                        marginal_surrogate,
                                        twostep = FALSE,
                                        method = "BFGS")
  coef(full_model$submodel0$ml_fit)
  coefs_check = c(0.00773068, 1.08557400, 0.32093667, 1.88178947)
  expect_equal(coef(full_model$submodel0$ml_fit), coefs_check, ignore_attr = TRUE)
}
)

test_that("fit_model_binCont_copula() works with clayton copula and normal margins and full ML estimator", {
  copula_family = "clayton"
  marginal_surrogate = "normal"
  data("Schizo_BinCont")
  na = is.na(Schizo_BinCont$CGI_Bin) | is.na(Schizo_BinCont$PANSS)
  X = Schizo_BinCont$PANSS[!na]
  Y = Schizo_BinCont$CGI_Bin[!na]
  Treat = Schizo_BinCont$Treat[!na]
  Treat = ifelse(Treat == 1, 1, 0)
  data = data.frame(X,
                    Y,
                    Treat)
  full_model = fit_copula_model_BinCont(data,
                                        copula_family,
                                        marginal_surrogate,
                                        twostep = FALSE,
                                        method = "BFGS")
  coef(full_model$submodel0$ml_fit)
  coefs_check = c(-0.04323186 ,-16.63365493, 28.63091060, 3.14363815)
  expect_equal(
    coef(full_model$submodel0$ml_fit),
    coefs_check,
    ignore_attr = TRUE,
    tolerance = 10e-3
  )
}
)
