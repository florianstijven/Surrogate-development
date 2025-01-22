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
    -38.082448704,
    tolerance = 0.01
  )
  expect_equal(
    continuous_continuous_loglik(c(para[-1], 0.5), X, Y, "gaussian", marginal_X, marginal_Y),
    -139.6712877,
    tolerance = 0.01
  )
  expect_equal(
    continuous_continuous_loglik(para, X, Y, "frank", marginal_X, marginal_Y),
    -32.57395345,
    tolerance = 0.01
  )
  expect_equal(
    continuous_continuous_loglik(para, X, Y, "gumbel", marginal_X, marginal_Y),
    -36.81896141,
    tolerance = 0.01
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
    -25.3366761945,
    ignore_attr = "df",
    tolerance = 0.01
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
    -275.28200947,
    ignore_attr = "df",
    tolerance = 0.1
  )
  expect_equal(
    fitted_submodel$marginal_X$pdf(1:5),
    c(0.0653530354649, 0.1218597549174, 0.1768435023308, 0.1997340854162, 0.1755696866447),
    ignore_attr = "names",
    tolerance = 0.01
  )
  expect_equal(
    fitted_submodel$marginal_Y$inv_cdf((1:5) / 5),
    c(-0.963708006673, 0.178103097350,1.161571587856, 2.303382691879, Inf),
    tolerance = 0.01
  )
})

test_that("fit_copula_ContCont() works", {
  S0 = c(1, 6, 2, 5, 3, 6, 4)
  S0 = rep(S0, 10)
  S1 = c(1, 2, 2, 5, 3, 2, 4)
  S1 = rep(S1, 10)

  T0 = c(2.2, 3.1, 0, -3, 0, 1, 4)
  T0 = rep(T0, 10)
  T1 = c(0, 3.1, 0.5, -3, 2, 1, 1)
  T1 = rep(T1, 10)

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

  fitted_model = fit_copula_ContCont(
    data = data,
    copula_family = "frank",
    marginal_S0 = marginal,
    marginal_S1 = marginal,
    marginal_T0 = marginal,
    marginal_T1 = marginal,
    start_copula = 2
  )
  expect_equal(
    fitted_model$fit_0$ml_fit$maximum,
    -294.088304226
  )
}
)

test_that("fit_copula_ContCont() works with the non-central ", {
  data("Schizo")
  Schizo = Schizo[1:100, ]
  na = is.na(Schizo$BPRS) | is.na(Schizo$PANSS)
  X = Schizo$BPRS[!na]
  Y = Schizo$PANSS[!na]
  Treat = Schizo$Treat[!na]
  Treat = ifelse(Treat == 1, 1, 0)
  data = data.frame(X,
                    Y,
                    Treat)


  marginal_S = list(
    pdf_fun = function(x, para) {
      sn::dsn(x, xi = para[1], omega = para[2], tau = para[3])
    },
    cdf_fun = function(x, para) {
      sn::psn(x, xi = para[1], omega = para[2], tau = para[3])
    },
    q_fun = function(p, para) {
      sn::qsn(x, xi = para[1], omega = para[2], tau = para[3])
    },
    n_para = 3,
    start = c(mean(data$X), sd(data$X), 0)
  )

  marginal_T = list(
    pdf_fun = function(x, para) {
      sn::dsn(x, xi = para[1], omega = para[2], tau = para[3])
    },
    cdf_fun = function(x, para) {
      sn::psn(x, xi = para[1], omega = para[2], tau = para[3])
    },
    q_fun = function(p, para) {
      sn::qsn(x, xi = para[1], omega = para[2], tau = para[3])
    },
    n_para = 3,
    start = c(mean(data$Y), sd(data$Y), 0)
  )

  fitted_model = fit_copula_ContCont(
    data = data,
    copula_family = "clayton",
    marginal_S0 = marginal_S,
    marginal_S1 = marginal_S,
    marginal_T0 = marginal_T,
    marginal_T1 = marginal_T,
    start_copula = 5,
    method = "BHHH"
  )
  # Maximized loglikelihoods match
  expect_equal(
    c(
      fitted_model$fit_0$ml_fit$maximum,
      fitted_model$fit_1$ml_fit$maximum
    ),
    c(-296.933821938, -530.076269924),
    tolerance = 0.01
  )
  # Estimated coefficients match
  expect_equal(
    c(
      coef(fitted_model$fit_0$ml_fit),
      coef(fitted_model$fit_1$ml_fit)
    ),
    c(-15.111858211,
      21.802586887,
      0.313607345,
      -19.593466595,
      37.031933928,
      0.007206932,
      4.692048418,
      -16.9450853,
      23.8469725,
      0.1226650,
      -25.5717622,
      42.4446427,
      -0.2563797,
      8.7131708
    ),
    ignore_attr = "names",
    tolerance = 0.1
  )
}
)

# test_that("fit_copula_ContCont() works with the skewed normal ", {
#   data("Schizo")
#   na = is.na(Schizo$BPRS) | is.na(Schizo$PANSS)
#   X = Schizo$BPRS[!na]
#   Y = Schizo$PANSS[!na]
#   Treat = Schizo$Treat[!na]
#   Treat = ifelse(Treat == 1, 1, 0)
#   data = data.frame(X,
#                     Y,
#                     Treat)
#
#
#   marginal_S0 = list(
#     pdf_fun = function(x, para) {
#       sn::dsn(x, xi = para[1], omega = exp(para[2]), alpha = para[3])
#     },
#     cdf_fun = function(x, para) {
#       sn::psn(x, xi = para[1], omega = exp(para[2]), alpha = para[3])
#     },
#     q_fun = function(p, para) {
#       sn::qsn(x, xi = para[1], omega = exp(para[2]), alpha = para[3])
#     },
#     n_para = 3,
#     start = c(13.67444, 3.24893, -3.47)
#   )
#   marginal_S1 = marginal_S0
#   marginal_S1[[5]] = c(mean(data$X[data$Treat == 1]), log(sd(data$X[data$Treat == 1])), 0)
#
#   marginal_T0 = list(
#     pdf_fun = function(x, para) {
#       sn::dsn(x, xi = para[1], omega = exp(para[2]), alpha = para[3])
#     },
#     cdf_fun = function(x, para) {
#       sn::psn(x, xi = para[1], omega = exp(para[2]), alpha = para[3])
#     },
#     q_fun = function(p, para) {
#       sn::qsn(x, xi = para[1], omega = exp(para[2]), alpha = para[3])
#     },
#     n_para = 3,
#     start = c(23.14603, 3.78565, -3.16)
#   )
#   marginal_T1 = marginal_T0
#   marginal_T1[[5]] = c(mean(data$Y[data$Treat == 1]), log(sd(data$Y[data$Treat == 1])), 0)
#
#   fitted_model = fit_copula_ContCont(
#     data = data,
#     copula_family = "clayton",
#     marginal_S0 = marginal_S0,
#     marginal_S1 = marginal_S1,
#     marginal_T0 = marginal_T0,
#     marginal_T1 = marginal_T1,
#     start_copula = 6.5,
#     method = "BFGS",
#     finalHessian = "BHHH"
#   )
#   print.vine_copula_fit(fitted_model)
#   expect_equal(
#     fitted_model$fit_0$ml_fit$maximum,
#     -294.088304226
#   )
# }
# )

test_that("GoF functions work", {
  S0 = c(1, 6, 2, 5, 3, 6, 4)
  S0 = rep(S0, 10)
  S1 = c(1, 2, 2, 5, 3, 2, 4)
  S1 = rep(S1, 10)

  T0 = c(2.2, 3.1, 0, -3, 0, 1, 4)
  T0 = rep(T0, 10)
  T1 = c(0, 3.1, 0.5, -3, 2, 1, 1)
  T1 = rep(T1, 10)

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

  fitted_model = fit_copula_ContCont(
    data = data,
    copula_family = "frank",
    marginal_S0 = marginal,
    marginal_S1 = marginal,
    marginal_T0 = marginal,
    marginal_T1 = marginal,
    start_copula = 2
  )
  # Conditional mean function
  expect_equal(
    conditional_mean_copula_ContCont(fitted_model$fit_0, grid = 1:5),
    c(0.924384620663, 0.952002129594, 1.000133909041, 1.062347011343, 1.121985370807),
    tolerance = 0.01
  )
}
)
