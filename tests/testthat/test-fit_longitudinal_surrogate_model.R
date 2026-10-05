library(testthat)

test_that("fit_longitudinal_surrogate_model() works on a small example", {
  set.seed(1)

  raw_dat <- expand.grid(
    Id = 1:12,
    Time = c(1, 2),
    Endpoint = c("BPRS", "PANSS")
  )

  subject_treat <- rep(c(0, 1), each = 6)
  raw_dat$Treat <- subject_treat[raw_dat$Id]

  raw_dat$Baseline <- ifelse(raw_dat$Endpoint == "BPRS", 40, 80) +
    rnorm(nrow(raw_dat), 0, 2)
  raw_dat$Response <- raw_dat$Baseline -
    2 * raw_dat$Treat -
    0.5 * raw_dat$Time +
    ifelse(raw_dat$Endpoint == "PANSS", 5, 0) +
    rnorm(nrow(raw_dat), 0, 1)

  raw_dat$Response_log <- log(raw_dat$Response)
  raw_dat$Baseline_log <- log(raw_dat$Baseline)
  raw_dat$Time_cont <- raw_dat$Time

  dat <- prepare_longitudinal_surrogate_data(
    data = raw_dat,
    id = "Id",
    treatment = "Treat",
    time = "Time_cont",
    endpoint = "Endpoint",
    outcome = "Response_log",
    baseline = "Baseline_log",
    surrogate_label = "BPRS",
    true_label = "PANSS",
    treatment_levels = c(0, 1)
  )

  mean_formula <- Y ~ endpoint + visitF + treatF + Ybase

  spec_galecki <- make_longitudinal_re_none(
    name = "Galecki's marginal model"
  )

  fit <- fit_longitudinal_surrogate_model(
    data = dat,
    formula = mean_formula,
    model_spec = spec_galecki,
    residual_type = "independence",
    control = list(eval.max = 100, iter.max = 100)
  )

  expect_s3_class(fit, "surrogate_longitudinal_fit")
  expect_equal(fit$model_spec$name, "Galecki's marginal model")
  expect_equal(fit$residual_type, "independence")

  pars <- extract_longitudinal_parameters(fit)
  expect_true(is.list(pars))
  expect_true(all(c("V0", "V1") %in% names(pars)))

  stats <- longitudinal_fit_stats(fit)
  expect_s3_class(stats, "data.frame")
  expect_equal(nrow(stats), 1)
  expect_true(all(c("logLik", "AIC", "convergence") %in% names(stats)))

  diagnostics <- longitudinal_matrix_diagnostics(fit)
  expect_s3_class(diagnostics, "data.frame")
  expect_true(all(c("matrix", "determinant", "min_eig", "status") %in% names(diagnostics)))
})
