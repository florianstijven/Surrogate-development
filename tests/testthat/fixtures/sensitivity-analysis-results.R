# Load fitted copula model.
fitted_model = readRDS(test_path("fixtures", "ovarian-dvine-clayton.rds"))
# Perform sensitivity analysis with a small number of replications
set.seed(1)
sens_results = sensitivity_analysis_SurvSurv_copula(
  fitted_model,
  composite = TRUE,
  cond_ind = TRUE,
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
  cond_ind = TRUE,
  n_sim = 20,
  n_prec = 1e4
)
# Save results to file.
saveRDS(sens_results_comp, file = "tests/testthat/fixtures/sens-results-ovarian-gaussian-comp.rds")
