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
    n_para = 2
  )
  K = 6
  expect_equal(
    ordinal_continuous_loglik(para, X, Y, copula_family, marginal_Y, K),
    -40.644492067
  )
})
