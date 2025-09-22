# test_that("Survival-continuous model works for prostate data", {
#   skip_on_cran()
#   data("prostate")
#
#   fitted_model = MetaAnalyticSurvCont(
#     data = prostate,
#     true = SURVTIME,
#     trueind = SURVIND,
#     surrog = PSA,
#     trt = TREAT,
#     center = TRIAL,
#     trial = TRIAL,
#     patientid = PATID,
#     copula = "Hougaard",
#     adjustment = "weighted"
#   )
#   rsquared_theta = c(
#     round(fitted_model$Trial.R2[1, 1], digits = 4),
#     round(fitted_model$Indiv.Surrogacy[1, 1], digits = 4)
#   )
#   expect_equal(rsquared_theta, c(0.0066, 0.2764), tolerance =
#                  1e-2)
# })
