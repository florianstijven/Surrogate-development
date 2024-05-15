test_that("Survival- Survival model works for ovarian data",
          {
            skip_on_cran()
            data("Ovarian")

            fitted_model = MetaAnalyticSurvSurv(data=Ovarian,true=Surv,trueind=SurvInd,surrog=Pfs,surrogind=PfsInd,
                                                trt=Treat,center=Center,trial=Center,patientid=Patient,
                                                copula="Plackett",adjustment="unadjusted")
            rsquared_theta = c(round(fitted_model$Trial.R2[1,1], digits = 4), round(fitted_model$Indiv.Surrogacy[1,1], digits = 4))
            expect_equal(rsquared_theta, c(0.8721,  0.9671), tolerance=1e-2)
          })
