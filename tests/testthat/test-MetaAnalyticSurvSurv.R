test_that("Survival- Survival model works for ovarian data",
          {
            skip_on_cran()
            data("Ovarian")

            Ovarian2 <- Ovarian %>%
              group_by(Center) %>%
              filter(sum(Treat == 1) >= 3 & sum(Treat == 0) >= 3) %>%
              ungroup()

            fitted_model =  MetaAnalyticSurvSurv(data=Ovarian2,true=Surv,trueind=SurvInd,surrog=Pfs,surrogind=PfsInd,
                                                 trt=Treat,center=Center,trial=Center,patientid=Patient,
                                                 copula="Plackett",adjustment="weighted")
            rsquared_theta = c(round(fitted_model$Trial.R2[1,1], digits = 4), round(fitted_model$Indiv.Surrogacy[1,1], digits = 4))
            expect_equal(rsquared_theta, c(0.8721,  0.9671), tolerance=1e-2)
          })
