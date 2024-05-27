
#' Fits the second stage model in the two-stage federated data analysis approach.
#'
#'The function 'FederatedApproachStage2()' fits the second stage model of the two-stage federated data analysis approach to assess surrogacy.
#'
#' @details
#'
#' # Model
#'
#' The two-stage federated data analysis approach developed by XXX can be used to assess surrogacy in the meta-analytic multiple-trial setting
#' (Continuous-continuous case), but without the need of sharing data. Instead, each organization conducts separate analyses on their data set
#' using the so-called "first stage" model. The results of these analyses are then aggregated at a central analysis hub,
#' where the aggregated results are analyzed using a "second stage" model and the necessary metrics (\eqn{R^2_{trial}} and \eqn{R^2_{indiv}})
#' for the validation of the surrogate endpoint are obtained. This function fits the second stage model, where a method-of-moments estimator is
#' used to obtain the variance-covariance matrix \eqn{D} from which the \eqn{R^2_{trial}} can be derived. The \eqn{R^2_{indiv}} is obtained with
#' a weighted average of the elements in \eqn{\Sigma_i}.
#'
#' # Data Format
#'
#' A data frame that combines the results of the first stage models and contains:
#'
#' * a column with the trial indicator
#' * a column with the number of subjects in the trial
#' * a column with the estimated intercepts for the surrogate
#' * a column with the estimated treatment effects for the surrogate
#' * a column with the estimated intercepts for the true endpoint
#' * a column with the estimated treatment effects for the true endpoint
#' * a column with the variances of the error term for the surrogate endpoint
#' * a column with the covariances between the error terms of the surrogate and true endpoint
#' * a column with the variances of the error term for the true endpoint
#'
#' A list that combines all the variance-covariance matrices of the fixed effects obtained using the first stage model
#'
#' @references Florez, A. J., Molenberghs G, Verbeke G, Alonso, A. (2019). A closed-form estimator for metaanalysis and surrogate markers evaluation. Journal of Biopharmaceutical Statistics, 29(2) 318-332.
#'
#' @param Dataset A data frame with the correct columns (See Data Format).
#' @param Intercept.S Estimated intercepts for the surrogate endpoint.
#' @param alpha Estimated treatment effects for the surrogate endpoint.
#' @param Intercept.T Estimated intercepts for the true endpoint.
#' @param beta Estimated treatment effects for the true endpoint.
#' @param sigma.SS Estimated variance of the error terms for the surrogate endpoint.
#' @param sigma.ST Estimated covariance between the error terms of the surrogate and true endpoint.
#' @param sigma.TT Estimated variance of the error terms for the true endpoint.
#' @param Obs.per.trial Number of subjects in the trial.
#' @param Trial.ID Trial indicator.
#' @param R.list List of the variance-covariance matrices of the fixed effects.
#' @param Alpha The \eqn{\alpha}-level that is used to determine the confidence intervals around \eqn{R^2_{trial}} and \eqn{R^2_{indiv}}. Default 0.05.
#'
#' @return Returns an object of class "FederatedApproachStage2()" that can be used to evaluate surrogacy.
#'
#' * Indiv.R2: a data frame that contains the \eqn{R^2_{indiv}} and 95% confidence interval to evaluate surrogacy at the trial level.
#' * Trial.R2: a data frame that contains the \eqn{R^2_{trial}} and 95% confidence interval to evaluate surrogacy at the trial level.
#' * Fixed.Effects: a data frame that contains the average of the estimated fixed effects.
#' * D: estimated \eqn{D} matrix.
#' * Obs.Per.Trial: number of observations in each trial.
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
FederatedApproachStage2 <- function(Dataset, Intercept.S,
                                    alpha, Intercept.T, beta,
                                    sigma.SS, sigma.ST, sigma.TT,
                                    Obs.per.trial, Trial.ID, R.list, Alpha=0.05){
  R2TrialFun <- function(D) {
    A <- matrix(c(D[1, 4], D[2, 4]), 2, 1)
    B <- matrix(c(D[1, 1], D[1, 2], D[1, 2], D[2, 2]), 2,
                2)
    C <- D[4, 4]
    R2.trial <- crossprod(A, solve(B, A))/C
    return(R2.trial)
  }
  R2IndFun <- function(Sigma) {
    R2.ind <- Sigma[1, 2]^2/(Sigma[1, 1] * Sigma[2, 2])
    return(R2.ind)
  }
  pdDajustment = function(D) {
    EigenD <- eigen(D)
    E <- diag(EigenD$values)
    diag(E)[diag(E) <= 0] <- 1e-04
    L <- EigenD$vector
    DH <- L %*% tcrossprod(E, L)
    return(DH)
  }

  Dataset = Dataset[, c(paste(substitute(Trial.ID)), paste(substitute(Obs.per.trial)), paste(substitute(Intercept.S)),
                        paste(substitute(alpha)), paste(substitute(Intercept.T)), paste(substitute(beta)),
                        paste(substitute(sigma.SS)), paste(substitute(sigma.ST)), paste(substitute(sigma.TT)))]


  n <- Dataset[,2]
  Trial.ID <- Dataset[,1]
  N = length(n)

  BetaH_i = Dataset[,3:6]
  SigmaH_i = Dataset[,-c(1:6)]
  a_i = n/sum(n)
  a_i.2 = (n - 2)/sum(n - 2)
  BetaH = apply(t(BetaH_i), 1, weighted.mean, w = a_i)
  SigmaH = ks::invvech(apply(t(SigmaH_i), 1, weighted.mean, w = a_i.2))
  N = length(a_i)
  b_i <- t(BetaH_i) - tcrossprod(BetaH, matrix(1, N, 1))
  Sb = tcrossprod(b_i)
  num = 1 - 2 * a_i + sum(a_i^2)
  denom = sum(num)
  DH = (Sb - Reduce("+", mapply(function(X, x) {
    X * x
  }, R.list, num, SIMPLIFY = F)))/denom
  DH.pd = min(eigen(DH, only.values = T)$values) > 0
  if (!DH.pd) {
    DH = pdDajustment(DH)
    warning(paste("The estimate of D is non-positive definite. Adjustment for non-positive definiteness was required"))
  }
  Var.BetaH_i = mapply("+", R.list, MoreArgs = list(DH = DH),
                       SIMPLIFY = F)
  VarBetaH_i.inv = mapply(solve, Var.BetaH_i, SIMPLIFY = FALSE)
  A_i.denom = solve(Reduce("+", VarBetaH_i.inv))
  A_i = mapply("%*%", list(A_i.denom), VarBetaH_i.inv, SIMPLIFY = FALSE)
  BetaH = apply(mapply(function(x) {
    A_i[[x]] %*% t(BetaH_i)[, x]
  }, x = 1:N), 1, sum)
  VarBetaH = Reduce("+", mapply(tcrossprod, mapply("%*%",
                                                   A_i, Var.BetaH_i, SIMPLIFY = FALSE), A_i, SIMPLIFY = FALSE))
  R2trial = R2TrialFun(DH)
  R2ind = R2IndFun(SigmaH)
  R2ind.sd <- sqrt((4 * R2ind * ((1 - R2ind)^2))/(sum(n) - 3))
  R2ind.lb <- max(0, R2ind + qnorm(Alpha/2) * R2ind.sd)
  R2ind.ub <- min(1, R2ind + qnorm(1 - Alpha/2) * R2ind.sd)
  Indiv.R2 <- data.frame(cbind(R2ind, R2ind.sd, R2ind.lb,
                               R2ind.ub), stringsAsFactors = TRUE)
  colnames(Indiv.R2) <- c("R2 Indiv", "Standard Error", "CI lower limit",
                          "CI upper limit")
  R2trial.sd <- sqrt((4 * R2trial * (1 - R2trial)^2)/(N -
                                                        3))
  R2trial.lb <- max(0, R2trial + qnorm(Alpha/2) * (R2trial.sd))
  R2trial.ub <- min(1, R2trial + qnorm(1 - Alpha/2) * (R2trial.sd))
  Trial.R2 <- data.frame(cbind(R2trial, R2trial.sd, R2trial.lb,
                               R2trial.ub), stringsAsFactors = TRUE)
  colnames(Trial.R2) <- c("R2 Trial", "Standard Error", "CI lower limit",
                          "CI upper limit")
  rownames(Trial.R2) = rownames(Indiv.R2) <- c(" ")
  Obs.Per.Trial <- as.data.frame(t(rbind(Trial.ID, n)))
  Fixed.Effects = cbind(BetaH, sqrt(diag(VarBetaH)))[c(3,
                                                       1, 4, 2), ]
  rownames(Fixed.Effects) = c("mu_S", "mu_T", "alpha", "beta")
  colnames(Fixed.Effects) = c("estimate", "standard error")
  Fixed.Effects = as.data.frame(Fixed.Effects)
  DH = DH[c(3, 1, 4, 2), c(3, 1, 4, 2)]
  colnames(DH) = rownames(DH) = c("mu_S_i", "mu_T_i", "a_i",
                                  "b_i")
  SigmaH = SigmaH
  colnames(SigmaH) = rownames(SigmaH) = c("S_ij", "T_ij")

  Output = list(Obs.Per.Trial = Obs.Per.Trial,
                Fixed.Effects = Fixed.Effects, Trial.R2 = Trial.R2,
                Indiv.R2 = Indiv.R2, D = DH,
                Call = match.call())
  class(Output) <- "FederatedApproachStage2"

  return(Output)

}

#' Provides a summary of the surrogacy measures for an object fitted with the 'FederatedApproachStage2()' function.
#'
#' @method summary FederatedApproachStage2
#'
#' @param object An object of class 'FederatedApproachStage2' fitted with the 'FederatedApproachStage2()' function.
#' @param ... ...
#'
#' @return The surrogacy measures with their 95% confidence intervals.
#' @export
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

summary.FederatedApproachStage2 <- function(object,...){
  cat("Surrogacy measures with 95% confidence interval \n\n")
  cat("Individual level surrogacy: ", "\n\n")
  cat("R Square", ": ", sprintf("%.4f", object$Indiv.R2[1,1]), "[", sprintf("%.4f", object$Indiv.R2[1,3]),";", sprintf("%.4f", object$Indiv.R2[1,4]) , "]", "\n\n")
  cat("Trial level surrogacy: ", "\n\n")
  cat("R Square: ", sprintf("%.4f", object$Trial.R2[1,1]),"[", sprintf("%.4f", object$Trial.R2[1,3]),";", sprintf("%.4f", object$Trial.R2[1,4]) , "]", "\n\n")
}
