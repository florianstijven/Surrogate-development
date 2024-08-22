# Load fitted copula model.
fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
# Perform sensitivity analysis with a small number of replications
set.seed(1)
sens_results = sensitivity_analysis_SurvSurv_copula(
  fitted_model,
  composite = TRUE,
  lower = c(-1, 0, 0, -1),
  upper = c(1, 0, 0, 1),
  eq_cond_association = FALSE,
  n_sim = 20,
  n_prec = 1e4
)
# Save results to file.
saveRDS(sens_results, file = "tests/testthat/fixtures/sens-results-ovarian-clayton.rds")


# Load fitted copula model.
fitted_model_comp = readRDS(test_path("fixtures", "ovarian-dvine-gaussian-scr.rds"))
# Perform sensitivity analysis with a small number of replications
set.seed(1)
sens_results_comp = sensitivity_analysis_SurvSurv_copula(
  fitted_model_comp,
  composite = TRUE,
  lower = c(-1, 0, 0, -1),
  upper = c(1, 0, 0, 1),
  eq_cond_association = FALSE,
  n_sim = 20,
  n_prec = 1e4
)
# Save results to file.
saveRDS(sens_results_comp, file = "tests/testthat/fixtures/sens-results-ovarian-gaussian-comp.rds")

# Load fitted copula model.
fitted_model_comp = readRDS(test_path("fixtures", "ovarian-dvine-gaussian-scr.rds"))
# Perform sensitivity analysis with a small number of replications
set.seed(1)
sens_results_comp = sensitivity_analysis_SurvSurv_copula(
  fitted_model_comp,
  composite = TRUE,
  lower = c(-1, 0, 0, -1),
  upper = c(1, 0, 0, 1),
  eq_cond_association = FALSE,
  n_sim = 20,
  n_prec = 1e4,
  restr_time = 1,
  mutinfo_estimator = function(x, y) -0.5 *  log(1 - stats::cor(x, y, method = "spearman"))
)
# Save results to file.
saveRDS(sens_results_comp, file = "tests/testthat/fixtures/sens-results-ovarian-gaussian-comp-sprho-restr.rds")

# Load fitted copula models.
fitted_model_contcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-ContCont.rds"))
fitted_model_ordcont = readRDS(test_path("fixtures", "schizo-dvine-clayton-OrdCont.rds"))
fitted_model_ordord = readRDS(test_path("fixtures", "schizo-dvine-gaussian-OrdOrd.rds"))
# Perform sensitivity analyses with a small number of replications.
set.seed(1)
sens_results_contcont = sensitivity_analysis_copula(
  fitted_model = fitted_model_contcont,
  n_sim = 20,
  eq_cond_association = TRUE,
  lower = c(0, 0, 0, 0),
  upper = c(0.3, 0, 0, 0.95),
  degrees = 0,
  marg_association = TRUE,
  copula_family2 = "clayton",
  n_prec = 1e4,
  ncores = 1
)
set.seed(1)
sens_results_ordcont = sensitivity_analysis_copula(
  fitted_model = fitted_model_ordcont,
  n_sim = 20,
  eq_cond_association = TRUE,
  lower = c(0, 0, 0, 0),
  upper = c(0.3, 0, 0, 0.95),,
  degrees = 0,
  marg_association = TRUE,
  copula_family2 = "clayton",
  n_prec = 1e4,
  ncores = 1
)
set.seed(1)
sens_results_ordord = sensitivity_analysis_copula(
  fitted_model = fitted_model_ordord,
  n_sim = 20,
  eq_cond_association = TRUE,
  lower = c(0, 0, 0, 0),
  upper = c(0.3, 0, 0, 0.95),
  degrees = 0,
  marg_association = TRUE,
  copula_family2 = "clayton",
  n_prec = 1e4,
  ncores = 1
)
# Save results to file.
saveRDS(sens_results_contcont, file = "tests/testthat/fixtures/sens-results-schizo-clayton-ContCont.rds")
saveRDS(sens_results_ordcont, file = "tests/testthat/fixtures/sens-results-schizo-clayton-OrdCont.rds")
saveRDS(sens_results_ordord, file = "tests/testthat/fixtures/sens-results-schizo-gaussian-OrdOrd.rds")
