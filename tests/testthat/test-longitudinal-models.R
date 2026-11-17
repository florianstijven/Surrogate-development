test_that("standard longitudinal model specifications are constructed", {
  expect_identical(longitudinal_model_spec("galecki")$re_mode, "none")
  expect_identical(
    longitudinal_model_spec("endpoint_ri")$re_struct$T_0,
    "ri"
  )
  expect_identical(
    longitudinal_model_spec("endpoint_ri_rs")$re_struct$S_1,
    "rislope"
  )
  expect_identical(
    longitudinal_model_spec("shared_ri")$shared_re[["0"]],
    "ri"
  )
  expect_identical(
    longitudinal_model_spec("shared_ri_rs")$shared_re[["1"]],
    "rislope"
  )
  expect_error(longitudinal_model_spec("ri_nested"))
})
