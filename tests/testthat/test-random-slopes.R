test_that("endpoint RI/RS design uses the fitted time centering", {
  object <- structure(
    list(
      model_spec = longitudinal_model_spec("endpoint_ri_rs"),
      endpoint_order = c("S", "T"),
      time_center = 4
    ),
    class = "surrogate_longitudinal_fit"
  )
  evaluation <- list(T = c(2, 6), S = c(1, 8))
  design <- Surrogate:::.surr_evaluation_design(object, evaluation, 0)$Z

  expect_equal(dim(design), c(4, 4))
  expect_equal(
    colnames(design),
    c("ri_S", "rs_S", "ri_T", "rs_T")
  )
  expect_equal(design[1:2, "rs_T"], c(-2, 2))
  expect_equal(design[3:4, "rs_S"], c(-3, 4))
})

test_that("invalid RI/RS cross-world completions are removed", {
  D0 <- diag(2)
  D1 <- diag(2)
  valid <- matrix(0, 2, 2)
  invalid <- matrix(2, 2, 2)
  out <- Surrogate:::.surr_validate_random_completions(
    D0, D1, list(valid, invalid)
  )
  expect_length(out, 1)
  expect_equal(out[[1]], valid)
})
