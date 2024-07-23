test_that("GoF plots for OrdCont work", {
  # Load fitted copula model.
  fitted_OrdCont = readRDS(test_path("fixtures", "schizo-dvine-clayton-OrdCont.rds"))

  # Compare the plotted functions with a snapshot. Note that this
  # comparison is fragile, and unimportant changes could lead to a
  # failure of this test.
  vdiffr::expect_doppelganger(
    title = "Marginal GoF T0 - OrdCont - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_OrdCont$fit_0$marginal_X, fitted_OrdCont$fit_0$data$X, fitted_OrdCont$fit_0$names_XY[1], fitted_OrdCont$endpoint_types[1], 0)
    }
  )
  vdiffr::expect_doppelganger(
    title = "Marginal GoF T1 - OrdCont - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_OrdCont$fit_1$marginal_X, fitted_OrdCont$fit_1$data$X, fitted_OrdCont$fit_1$names_XY[1], fitted_OrdCont$endpoint_types[1], 1)
    }
  )

  vdiffr::expect_doppelganger(
    title = "Marginal GoF S0 - OrdCont - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_OrdCont$fit_0$marginal_Y, fitted_OrdCont$fit_0$data$Y, fitted_OrdCont$fit_0$names_XY[2], fitted_OrdCont$endpoint_types[2], 0)
    }
  )
  vdiffr::expect_doppelganger(
    title = "Marginal GoF S1 - OrdCont - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_OrdCont$fit_1$marginal_Y, fitted_OrdCont$fit_1$data$Y, fitted_OrdCont$fit_1$names_XY[2], fitted_OrdCont$endpoint_types[2], 1)
    }
  )

  vdiffr::expect_doppelganger(
    title = "Association GoF 0 - OrdCont - Schizo",
    fig = function() {
      association_gof_copula(fitted_OrdCont$fit_0, 0, c("ordinal", "continuous"))
    }
  )
  vdiffr::expect_doppelganger(
    title = "Association GoF 1 - OrdCont - Schizo",
    fig = function() {
      association_gof_copula(fitted_OrdCont$fit_1, 1, c("ordinal", "continuous"))
    }
  )
})

test_that("GoF plots for ContCont work", {
  # Load fitted copula model.
  fitted_ContCont = readRDS(test_path("fixtures", "schizo-dvine-clayton-ContCont.rds"))

  # Compare the plotted functions with a snapshot. Note that this
  # comparison is fragile, and unimportant changes could lead to a
  # failure of this test.
  vdiffr::expect_doppelganger(
    title = "Marginal GoF T0 - ContCont - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_ContCont$fit_0$marginal_X, fitted_ContCont$fit_0$data$X, fitted_ContCont$fit_0$names_XY[1], fitted_ContCont$endpoint_types[1], 0)
    }
  )
  vdiffr::expect_doppelganger(
    title = "Marginal GoF T1 - ContCont - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_ContCont$fit_1$marginal_X, fitted_ContCont$fit_1$data$X, fitted_ContCont$fit_1$names_XY[1], fitted_ContCont$endpoint_types[1], 1)
    }
  )

  vdiffr::expect_doppelganger(
    title = "Marginal GoF S0 - ContCont - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_ContCont$fit_0$marginal_Y, fitted_ContCont$fit_0$data$Y, fitted_ContCont$fit_0$names_XY[2], fitted_ContCont$endpoint_types[2], 0)
    }
  )
  vdiffr::expect_doppelganger(
    title = "Marginal GoF S1 - ContCont - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_ContCont$fit_1$marginal_Y, fitted_ContCont$fit_1$data$Y, fitted_ContCont$fit_1$names_XY[2], fitted_ContCont$endpoint_types[2], 1)
    }
  )

  vdiffr::expect_doppelganger(
    title = "Association GoF 0 - ContCont - Schizo",
    fig = function() {
      association_gof_copula(fitted_ContCont$fit_0, 0, c("continuous", "continuous"))
    }
  )
  vdiffr::expect_doppelganger(
    title = "Association GoF 1 - ContCont - Schizo",
    fig = function() {
      association_gof_copula(fitted_ContCont$fit_1, 1, c("continuous", "continuous"))
    }
  )
})

test_that("GoF plots for OrdOrd work", {
  # Load fitted copula model.
  fitted_OrdOrd = readRDS(test_path("fixtures", "schizo-dvine-gaussian-OrdOrd.rds"))

  # Compare the plotted functions with a snapshot. Note that this
  # comparison is fragile, and unimportant changes could lead to a
  # failure of this test.
  vdiffr::expect_doppelganger(
    title = "Marginal GoF T0 - OrdOrd - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_OrdOrd$fit_0$marginal_X, fitted_OrdOrd$fit_0$data$X, fitted_OrdOrd$fit_0$names_XY[1], fitted_OrdOrd$endpoint_types[1], 0)
    }
  )
  vdiffr::expect_doppelganger(
    title = "Marginal GoF T1 - OrdOrd - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_OrdOrd$fit_1$marginal_X, fitted_OrdOrd$fit_1$data$X, fitted_OrdOrd$fit_1$names_XY[1], fitted_OrdOrd$endpoint_types[1], 1)
    }
  )

  vdiffr::expect_doppelganger(
    title = "Marginal GoF S0 - OrdOrd - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_OrdOrd$fit_0$marginal_Y, fitted_OrdOrd$fit_0$data$Y, fitted_OrdOrd$fit_0$names_XY[2], fitted_OrdOrd$endpoint_types[2], 0)
    }
  )
  vdiffr::expect_doppelganger(
    title = "Marginal GoF S1 - OrdOrd - Schizo",
    fig = function() {
      marginal_gof_copula(fitted_OrdOrd$fit_1$marginal_Y, fitted_OrdOrd$fit_1$data$Y, fitted_OrdOrd$fit_1$names_XY[2], fitted_OrdOrd$endpoint_types[2], 1)
    }
  )

  vdiffr::expect_doppelganger(
    title = "Association GoF 0 - OrdOrd - Schizo",
    fig = function() {
      association_gof_copula(fitted_OrdOrd$fit_0, 0, c("ordinal", "ordinal"))
    }
  )
  vdiffr::expect_doppelganger(
    title = "Association GoF 1 - OrdOrd - Schizo",
    fig = function() {
      association_gof_copula(fitted_OrdOrd$fit_1, 1, c("ordinal", "ordinal"))
    }
  )
})
