library(testthat)

test_that("Schizo_PANSS_BPRS has the expected structure", {
  data(Schizo_PANSS_BPRS)

  expect_s3_class(Schizo_PANSS_BPRS, "data.frame")
  expect_equal(ncol(Schizo_PANSS_BPRS), 9)
  expect_equal(nrow(Schizo_PANSS_BPRS), 3962)

  expect_named(
    Schizo_PANSS_BPRS,
    c("Id", "Treat", "Time", "Endpoint", "Response",
      "Response_log", "Time_cont", "Baseline", "Baseline_log")
  )

  expect_equal(length(unique(Schizo_PANSS_BPRS$Id)), 450)
  expect_setequal(as.character(unique(Schizo_PANSS_BPRS$Endpoint)), c("BPRS", "PANSS"))
  expect_true(all(Schizo_PANSS_BPRS$Treat %in% c(0, 1)))
})

test_that("R2Lambda_from_covariance() and summary_longitudinal_ica() work", {
  Sigma_TT <- matrix(c(1, 0.3,
                       0.3, 1),
                     nrow = 2, byrow = TRUE)

  Sigma_SS <- matrix(c(1, 0.2,
                       0.2, 1),
                     nrow = 2, byrow = TRUE)

  Sigma_TS <- matrix(c(0.4, 0.1,
                       0.1, 0.5),
                     nrow = 2, byrow = TRUE)

  r2 <- R2Lambda_from_covariance(
    Sigma_TT = Sigma_TT,
    Sigma_SS = Sigma_SS,
    Sigma_TS = Sigma_TS
  )

  expect_type(r2, "double")
  expect_length(r2, 1)
  expect_true(r2 >= 0)
  expect_true(r2 <= 1)

  out <- list(R2_Lambda = c(0.70, 0.80, 0.90, 0.95))

  summ <- summary_longitudinal_ica(out, threshold = 0.90)

  expect_s3_class(summ, "data.frame")
  expect_equal(nrow(summ), 1)
  expect_equal(summ$minimum, 0.70)
  expect_equal(summ$maximum, 0.95)
  expect_equal(summ$threshold, 0.90)
  expect_equal(summ$fraction_below_threshold, 0.5)
})
