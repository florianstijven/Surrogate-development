test_that("baseline is represented only as a time-zero response", {
  dat <- data.frame(
    subject = rep(1:2, each = 6),
    arm = rep(c("control", "active"), each = 6),
    week = rep(c(0, 0, 1, 1, 2, 2), 2),
    endpoint_name = rep(c("BPRS", "PANSS"), 6),
    response = seq_len(12)
  )

  out <- prepare_longitudinal_surrogate_data(
    data = dat,
    id = "subject",
    treatment = "arm",
    time = "week",
    endpoint = "endpoint_name",
    outcome = "response",
    surrogate_label = "BPRS",
    true_label = "PANSS",
    treatment_levels = c("control", "active")
  )

  expect_equal(sum(out$is_baseline), 4)
  expect_equal(sum(out$is_postbaseline), 8)
  expect_equal(out$ctime, out$time_cont)
  expect_false("Ybase" %in% names(out))
})

test_that("the single baseline convention is enforced", {
  dat <- data.frame(
    subject = c(1, 1, 1, 1),
    arm = 0,
    week = c(0, 0, 1, 1),
    endpoint_name = c("BPRS", "PANSS", "BPRS", "PANSS"),
    response = 1:4,
    old_baseline_column = 10:13
  )

  expect_error(
    prepare_longitudinal_surrogate_data(
      data = dat,
      id = "subject",
      treatment = "arm",
      time = "week",
      endpoint = "endpoint_name",
      outcome = "response",
      surrogate_label = "BPRS",
      true_label = "PANSS",
      baseline = "old_baseline_column"
    ),
    "unused argument"
  )
  dat$week[3] <- -1
  expect_error(
    prepare_longitudinal_surrogate_data(
      data = dat,
      id = "subject",
      treatment = "arm",
      time = "week",
      endpoint = "endpoint_name",
      outcome = "response",
      surrogate_label = "BPRS",
      true_label = "PANSS"
    ),
    "nonnegative"
  )
})

test_that("duplicated subject-endpoint-time rows are rejected", {
  dat <- data.frame(
    subject = c(1, 1, 1),
    arm = 0,
    week = c(0, 0, 1),
    endpoint_name = c("BPRS", "BPRS", "PANSS"),
    response = 1:3
  )
  expect_error(
    prepare_longitudinal_surrogate_data(
      data = dat,
      id = "subject",
      treatment = "arm",
      time = "week",
      endpoint = "endpoint_name",
      outcome = "response",
      surrogate_label = "BPRS",
      true_label = "PANSS"
    ),
    "at most one observation"
  )
})
