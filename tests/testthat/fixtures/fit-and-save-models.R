# D-vine models for Ovarian data ------------------------------------------
# Load the Ovarian data.
data("Ovarian")
# For simplicity, data is not recoded to semi-competing risks format, but is
# left in the composite event format.
data = data.frame(Ovarian$Pfs,
                  Ovarian$Surv,
                  Ovarian$Treat,
                  Ovarian$PfsInd,
                  Ovarian$SurvInd)
# Fit and save the D-vine copula model with Clayton copula
fitted_model = fit_model_SurvSurv(data = data,
                                  copula_family = "clayton",
                                  n_knots = 1)
saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-clayton.rds")

# Fit and save the D-vine copula model with Frank copula
fitted_model = fit_model_SurvSurv(data = data,
                                  copula_family = "frank",
                                  n_knots = 1)
saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-frank.rds")

# Fit and save the D-vine copula model with Gaussian copula
fitted_model = fit_model_SurvSurv(data = data,
                                  copula_family = "gaussian",
                                  n_knots = 1)
saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-gaussian.rds")


# D-vine model for Schizo data --------------------------------------------
copula_family = "clayton"
marginal_surrogate = "normal"
data("Schizo_BinCont")
na = is.na(Schizo_BinCont$CGI_Bin) | is.na(Schizo_BinCont$PANSS)
X = Schizo_BinCont$PANSS[!na]
Y = Schizo_BinCont$CGI_Bin[!na]
Treat = Schizo_BinCont$Treat[!na]
Treat = ifelse(Treat == 1, 1, 0)
data = data.frame(X,
                  Y,
                  Treat)
fitted_model = fit_copula_model_BinCont(data, copula_family, marginal_surrogate, twostep = FALSE)
saveRDS(fitted_model, file = "tests/testthat/fixtures/schizo-dvine-clayton.rds")


# D-vine copula model for Ovarian data with semi-competing risks ----------

# Recode the Ovarian data in the semi-competing risks format.
data_scr = data.frame(
  ttp = Ovarian$Pfs,
  os = Ovarian$Surv,
  treat = Ovarian$Treat,
  ttp_ind = ifelse(
    Ovarian$Pfs == Ovarian$Surv &
      Ovarian$SurvInd == 1,
    0,
    Ovarian$PfsInd
  ),
  os_ind = Ovarian$SurvInd
)
# Save data in semi-competing risks format.
saveRDS(data_scr, file = "tests/testthat/fixtures/ovarian-data-scr.rds")
# Fit and save the D-vine copula model with Clayton copula
fitted_model = fit_model_SurvSurv(data = data_scr,
                                  copula_family = "clayton",
                                  n_knots = 1)
saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-clayton-scr.rds")

fitted_model = fit_model_SurvSurv(data = data_scr,
                                  copula_family = "frank",
                                  n_knots = 1)
saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-frank-scr.rds")

fitted_model = fit_model_SurvSurv(data = data_scr,
                                  copula_family = "gaussian",
                                  n_knots = 1)
saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-gaussian-scr.rds")

fitted_model = fit_model_SurvSurv(data = data_scr,
                                  copula_family = "gumbel",
                                  n_knots = 1)
saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-gumbel-scr.rds")

