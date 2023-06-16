test_that("ica_SurvSurv_sens() works", {
  data("Ovarian")
  # For simplicity, data is not recoded to semi-competing risks format, but the
  # data are left in the composite event format.
  data = data.frame(
    Ovarian$Pfs,
    Ovarian$Surv,
    Ovarian$Treat,
    Ovarian$PfsInd,
    Ovarian$SurvInd
  )
  ovarian_fitted =
      fit_model_SurvSurv(data = data,
                         copula_family = "clayton",
                         nknots = 1)
  # Illustration with small number of replications and low precision
  set.seed(1)
  sens_results = ica_SurvSurv_sens(
    ovarian_fitted,
    n_sim = 5,
    n_prec = 500,
    copula_family2 = "clayton",
    minfo_prec = 2e3
  )
  output_vector = c(sens_results$kendall[1],
                    sens_results$minfo[2],
                    sens_results$c23[3])
  check_vector = c(0.93848497, 2.21133646, 1.66162051)
  expect_equal(output_vector, check_vector)
})
