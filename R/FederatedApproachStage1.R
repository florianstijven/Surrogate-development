
#' Fits the first stage model in the two-stage federated data analysis approach.
#'
#'The function 'FederatedApproachStage1()' fits the first stage model of the two-stage federated data analysis approach to assess surrogacy.
#'
#' @details
#'
#' # Model
#'
#' The two-stage federated data analysis approach developed by XXX can be used to assess surrogacy in the meta-analytic multiple-trial setting
#' (Continuous-continuous case), but without the need of sharing data. Instead, each organization conducts separate analyses on their data set
#' using a so-called "first stage" model. The results of these analyses are then aggregated at a central analysis hub,
#' where the aggregated results are analyzed using a "second stage" model and the necessary metrics (\eqn{R^2_{trial}} and \eqn{R^2_{indiv}})
#' for the validation of the surrogate endpoint are obtained. This function fits the first stage model, where a linear model is fitted,
#' allowing estimation of the fixed effects.
#'
#' # Data Format
#'
#' The data frame must contain the following columns:
#'
#' * a column with the true endpoint
#' * a column with the surrogate endpoint
#' * a column with the treatment indicator: 0 or 1
#' * a column with the trial indicator
#' * a column with the patient indicator
#'
#' @references Florez, A. J., Molenberghs G, Verbeke G, Alonso, A. (2019). A closed-form estimator for metaanalysis and surrogate markers evaluation. Journal of Biopharmaceutical Statistics, 29(2) 318-332.
#'
#' @param Dataset A data frame with the correct columns (See Data Format).
#' @param Surr Surrogate endpoint.
#' @param True True endpoint.
#' @param Treat Treatment indicator.
#' @param Trial.ID Trial indicator.
#' @param Min.Treat.Size The minimum number of patients in each group (control or experimental) that a
#' trial should contain to be included in the analysis. If the number of patients in a
#' group of a trial is smaller than the value specified by Min.Treat.Size, the data
#' of the trial are excluded from the analysis. Default 2.
#' @param Alpha The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}} and \eqn{R^2_{indiv}}. Default 0.05.
#'
#' @return Returns an object of class "FederatedApproachStage1()" that can be used to evaluate surrogacy in the second stage model and contains the following elements:
#'
#' * Results.Stage.1: a data frame that contains the estimated fixed effects and the elements of \eqn{\Sigma_i}.
#' * R.i: the variance-covariance matrix of the estimated fixed effects.
#'
#' @export
#'
#' @author Dries De Witte
#'
#' @examples
#' \dontrun{
#' #As an example, the federated data analysis approach can be applied to the Schizo data set
#' data(Schizo)
#' Schizo <-  Schizo[order(Schizo$InvestId, Schizo$Id),]
#' #Create separate datasets for each investigator
#' Schizo_datasets <- list()
#'
#' for (invest_id in 1:198) {
#' Schizo_datasets[[invest_id]] <- Schizo[Schizo$InvestId == invest_id, ]
#' assign(paste0("Schizo", invest_id), Schizo_datasets[[invest_id]])
#' }
#' #Fit the first stage model for each dataset separately
#' results_stage1 <- list()
#' invest_ids <- list()
#' i <- 1
#' for (invest_id in 1:198) {
#'   dataset <- Schizo_datasets[[invest_id]]
#'
#'   skip_to_next <- FALSE
#'   tryCatch(FederatedApproachStage1(dataset, Surr=CGI, True=PANSS, Treat=Treat, Trial.ID = InvestId,
#'                                    Min.Treat.Size = 5, Alpha = 0.05),
#'                                    error = function(e) { skip_to_next <<- TRUE})
#'   #if the trial does not have the minimum required number, skip to the next
#'   if(skip_to_next) { next }
#'
#'   results_stage1[[invest_id]] <- FederatedApproachStage1(dataset, Surr=CGI, True=PANSS, Treat=Treat,
#'                                                          Trial.ID = InvestId, Min.Treat.Size = 5,
#'                                                          Alpha = 0.05)
#'   assign(paste0("stage1_invest", invest_id), results_stage1[[invest_id]])
#'   invest_ids[[i]] <- invest_id #keep a list of ids with datasets with required number of patients
#'   i <- i+1
#' }
#'
#' invest_ids <- unlist(invest_ids)
#' invest_ids
#'
#' #Combine the results of the first stage models
#' for (invest_id in invest_ids) {
#'   dataset <- results_stage1[[invest_id]]$Results.Stage.1
#'   if (invest_id == invest_ids[1]) {
#'     all_results_stage1<- dataset
#'  } else {
#'     all_results_stage1 <- rbind(all_results_stage1,dataset)
#'   }
#' }
#'
#' all_results_stage1 #that combines the results of the first stage models
#'
#' R.list <- list()
#' i <- 1
#' for (invest_id in invest_ids) {
#'   R <- results_stage1[[invest_id]]$R.i
#'   R.list[[i]] <- as.matrix(R[1:4,1:4])
#'   i <- i+1
#' }
#'
#' R.list #list that combines all the variance-covariance matrices of the fixed effects
#'
#' fit <- FederatedApproachStage2(Dataset = all_results_stage1, Intercept.S = Intercept.S,
#'                                alpha = alpha, Intercept.T = Intercept.T, beta = beta,
#'                                sigma.SS = sigma.SS, sigma.ST = sigma.ST,
#'                                sigma.TT = sigma.TT, Obs.per.trial = n,
#'                                Trial.ID = Trial.ID, R.list = R.list)
#' summary(fit)
#' }
FederatedApproachStage1 <- function(Dataset, Surr, True, Treat, Trial.ID, Min.Treat.Size = 2,
                                    Alpha = 0.05){
  Dataset = Dataset[, c(paste(substitute(Trial.ID)), paste(substitute(Treat)),
                        paste(substitute(Surr)), paste(substitute(True)))]
  Dataset = Dataset[!apply(Dataset, 1, anyNA), ]
  Trial.ID <- Dataset[, 1]
  Treat <- Dataset[, 2]
  Surr <- Dataset[, 3]
  True <- Dataset[, 4]
  Cluster.ID = Trial.ID[order(Trial.ID)]
  n = c(table(Cluster.ID))
  N = length(n)

  Y_i = matrix(cbind(Surr, True), length(cbind(Surr, True))/2, 2)

  X_i = matrix(cbind(1, Treat[order(Trial.ID)]), length(cbind(1, Treat[order(Trial.ID)]))/2, 2)

  Y = Y_i
  Z = X_i
  n.trial = min(table(Z[, 2]))
  n = nrow(Y)
  if (n.trial < Min.Treat.Size) {
    stop("This trial does not have the minimum required number of patients per treatment arm.")
  }
  BetaH <- tryCatch(solve(crossprod(Z), crossprod(Z, Y)),
                    error = function(e) {
                      NULL
                    })
  if (is.null(BetaH)) {
    stop("This trial does not have the minimum required number of patients per treatment arm.")
  }
  e <- Y - Z %*% BetaH
  SigmaH <- crossprod(e)/(n - 2)
  Output_beta_sigma <- c(c(BetaH), ks::vech(SigmaH))
  Output_beta_sigma <- t(as.data.frame(Output_beta_sigma))
  colnames(Output_beta_sigma) <- c("Intercept.S", "alpha", "Intercept.T", "beta", "sigma.SS", "sigma.ST", "sigma.TT")
  rownames(Output_beta_sigma) <- NULL
  R_i = kronecker(SigmaH, solve(crossprod(X_i)))
  Trial.ID <- Trial.ID[1]
  Results.Stage.1 <- data.frame(Trial.ID, n, Output_beta_sigma, stringsAsFactors = TRUE)
  rownames(Results.Stage.1) <- NULL
  output_list <- list(Results.Stage.1 = Results.Stage.1, R.i = R_i,
                      Call = match.call())
  class(output_list) <- "FederatedApproachStage1"

  return(output_list)

}
