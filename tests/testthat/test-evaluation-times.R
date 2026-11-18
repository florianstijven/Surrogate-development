mock_fit_for_times <- function() {
  structure(
    list(
      observed_times = list(T = c(1, 2, 4), S = c(1, 2, 4)),
      observed_time_ranges = list(T = c(1, 4), S = c(1, 4))
    ),
    class = "surrogate_longitudinal_fit"
  )
}

test_that("common and endpoint-specific evaluation times are accepted", {
  fit <- mock_fit_for_times()
  common <- validate_longitudinal_evaluation_times(fit, c(1, 2, 4))
  expect_true(common$common)
  expect_equal(common$T, c(1, 2, 4))

  different <- validate_longitudinal_evaluation_times(
    fit,
    list(T = c(2, 4), S = c(1, 2, 4))
  )
  expect_false(different$common)
  expect_equal(different$T, c(2, 4))
  expect_equal(different$S, c(1, 2, 4))
})

test_that("common_observed is explicit and safe", {
  fit <- mock_fit_for_times()
  out <- validate_longitudinal_evaluation_times(fit, "common_observed")
  expect_equal(out$T, c(1, 2, 4))

  fit$observed_times$S <- c(1, 3, 4)
  expect_error(
    validate_longitudinal_evaluation_times(fit, "common_observed"),
    "schedules differ"
  )
})

test_that("extrapolation warns but continues by default", {
  fit <- mock_fit_for_times()
  expect_warning(
    out <- validate_longitudinal_evaluation_times(fit, c(1, 2, 5)),
    "extrapolated"
  )
  expect_equal(out$T, c(1, 2, 5))
  expect_silent(
    validate_longitudinal_evaluation_times(
      fit, c(1, 2, 5), extrapolation = "allow"
    )
  )
  expect_error(
    validate_longitudinal_evaluation_times(
      fit, c(1, 2, 5), extrapolation = "error"
    )
  )
})

test_that("duplicated and unsorted evaluation times are rejected", {
  fit <- mock_fit_for_times()
  expect_error(validate_longitudinal_evaluation_times(fit, c(1, 1, 2)))
  expect_error(validate_longitudinal_evaluation_times(fit, c(2, 1, 4)))
})

test_that("baseline time is excluded from every ICA evaluation grid", {
  fit <- mock_fit_for_times()
  expect_error(
    validate_longitudinal_evaluation_times(fit, c(0, 1, 2)),
    "strictly positive"
  )
  expect_error(
    validate_longitudinal_evaluation_times(
      fit, list(T = c(1, 2), S = c(0, 1, 2))
    ),
    "strictly positive"
  )
})
