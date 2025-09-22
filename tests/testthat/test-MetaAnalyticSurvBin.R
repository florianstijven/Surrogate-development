# test_that("Survival-binary model works for Colorectal data", {
#   skip_on_cran()
#   data("colorectal")
#
#   fitted_model = MetaAnalyticSurvBin(
#     data = colorectal,
#     true = surv,
#     trueind = SURVIND,
#     surrog = responder,
#     trt = TREAT,
#     center = CENTER,
#     trial = TRIAL,
#     patientid = patientid,
#     adjustment = "unadjusted"
#   )
#   rsquared_theta = c(
#     round(fitted_model$Trial.R2[1, 1], digits = 4),
#     round(fitted_model$Indiv.Surrogacy[1, 1], digits = 4)
#   )
#   expect_equal(rsquared_theta, c(0.4419, 4.9109), tolerance =
#                  1e-2)
# })
