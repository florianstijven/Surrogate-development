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

# # Fit and save the D-vine copula model with Frank copula
# fitted_model = fit_model_SurvSurv(data = data,
#                                   copula_family = "frank",
#                                   n_knots = 1)
# saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-frank.rds")

# # Fit and save the D-vine copula model with Gaussian copula
# fitted_model = fit_model_SurvSurv(data = data,
#                                   copula_family = "gaussian",
#                                   n_knots = 1)
# saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-gaussian.rds")

# Fit and save the D-vine copula model with Clayton and Gumbel copula and
# variable number of knots
fitted_model = fit_model_SurvSurv(data = data,
                                  copula_family = c("clayton", "gumbel"),
                                  n_knots = c(0, 1, 2, 1))
saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-variable.rds")


# D-vine model for Schizo data --------------------------------------------

copula_family = "clayton"
marginal_surrogate = "normal"

## Continuous-Continuous Setting
data("Schizo")
na = is.na(Schizo$BPRS) | is.na(Schizo$PANSS)
X = Schizo$BPRS[!na]
Y = Schizo$PANSS[!na]
Treat = Schizo$Treat[!na]
Treat = ifelse(Treat == 1, 1, 0)
data = data.frame(X,
                  Y,
                  Treat)

marginal_X = list(
  pdf_fun = function(x, para) {
    dnorm(x, mean = para[1], sd = para[2])
  },
  cdf_fun = function(x, para) {
    pnorm(x, mean = para[1], sd = para[2])
  },
  NA,
  n_para = 2,
  start = c(-25, 28)
)
fitted_model = fit_copula_ContCont(
  data = data,
  copula_family = copula_family,
  marginal_S0 = marginal_X,
  marginal_S1 = marginal_X,
  marginal_T0 = marginal_X,
  marginal_T1 = marginal_X,
  start_copula = 3
)
saveRDS(fitted_model, file = "tests/testthat/fixtures/schizo-dvine-clayton-ContCont.rds")

## Binary-Continuous setting
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


data$Y = data$Y + 1
fitted_model = fit_copula_OrdCont(
  data = data,
  copula_family = copula_family,
  marginal_S0 = marginal_X,
  marginal_S1 = marginal_X,
  K_T = 2,
  start_copula = 3
)
saveRDS(fitted_model, file = "tests/testthat/fixtures/schizo-dvine-clayton-OrdCont.rds")

## Binary-Binary setting
data("Schizo_Bin")
na = is.na(Schizo_Bin$CGI_Bin) | is.na(Schizo_Bin$PANSS_Bin)
X = Schizo_Bin$CGI_Bin[!na]
Y = Schizo_Bin$PANSS_Bin[!na]
Treat = Schizo_Bin$Treat[!na]
Treat = ifelse(Treat == 1, 1, 0)
data = data.frame(X,
                  Y,
                  Treat)

data$Y = data$Y + 1
data$X = data$X + 1
fitted_model = fit_copula_OrdOrd(
  data = data,
  copula_family = "gaussian",
  K_S = 2,
  K_T = 2,
  start_copula = -0.5
)
saveRDS(fitted_model, file = "tests/testthat/fixtures/schizo-dvine-gaussian-OrdOrd.rds")

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

# fitted_model = fit_model_SurvSurv(data = data_scr,
#                                   copula_family = "gumbel",
#                                   n_knots = 1)
# saveRDS(fitted_model, file = "tests/testthat/fixtures/ovarian-dvine-gumbel-scr.rds")

