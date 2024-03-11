test_that("Survival-ordinal model works for Colorectal4 data",
          {
            data("colorectal4")

            fitted_model = survcat(data = colorectal4, true = truend, trueind = trueind,
                                          surrog = surrogend, trt = treatn, center = center,
                                          trial = trialend, patientid = patid)
            rsquared_theta = c(round(fitted_model$Trial.R2[1,1], digits = 4), round(fitted_model$Indiv.GlobalOdds[1,1], digits = 4))
            expect_equal(rsquared_theta, c(0.0971,  0.1471), tolerance=1e-3)
          })
