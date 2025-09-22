#' Assess surrogacy using a Rényi divergence based family of metrics in the causal-inference single-trial setting in normal case
#'
#' The function [ICA_alpha_ContCont()] is a set of metrics to evaluate surrogacy. \eqn{ICA_{\alpha}} have the similar
#' mathematical properties with [ICA.ContCont()].
#'
#' @param alpha (numeric) is order `alpha in [0, infinity]`
#' @param T0S0 A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the control treatment condition
#' @param T1S1 A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the control treatment condition
#' @param T0T0 A scalar that specifies the variance of the true endpoint in the control treatment condition
#' @param T1T1 A scalar that specifies the variance of the true endpoint in the control treatment condition
#' @param S0S0 A scalar that specifies the variance of the true endpoint in the control treatment condition
#' @param S1S1 A scalar that specifies the variance of the true endpoint in the control treatment condition
#' @param T0T1 A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1
#' @param T0S1 A scalar or vector that contains the correlation(s) between the counterfactuals T0 and S1
#' @param T1S0 A scalar or vector that contains the correlation(s) between the counterfactuals T1 and S0
#' @param S0S1 A scalar or vector that contains the correlation(s) between the counterfactuals S0 and S1




#' @return
#'
#' * Total.Num.Matrices: An object of class numeric that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.
#' * Pos.Def: A data.frame that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the \eqn{\rho_{\Delta}} values.
#' * rho: A scalar or vector that contains the individual causal association \eqn{\rho_{\Delta}}
#' * ICA: A scalar or vector that contains the individual causal association \eqn{\rho_{\Delta}^2=ICA}
#' * ICA_alpha: A scalar or vector that contains the individual causal association \eqn{ICA_{\alpha}}
#' * Sigmas: A data.frame that contains the \eqn{\sigma_{\Delta T}} and \eqn{\sigma_{\Delta S}}
#'
#'
#' @details
#' # New family of surrogacy metrics:
#' ICA was extended as \eqn{ICA_{\alpha}} by proposing a new family of surrogacy metrics derived from Rényi divergence.
#' Based on this family of divergences, one can propose a more general set of surrogacy metrics defined as:
#'
#' \deqn{ ICA_{\alpha} = 1 - e^{-2 D_{\alpha}[f(\Delta T,\Delta S), f(\Delta T)f(\Delta S)]} }
#'
#' when \eqn{\alpha=1} then
#'
#' \eqn{ICA_{1}=1-e^{-2D_1[f(\Delta T,\Delta S),f(\Delta T)f(\Delta S)]}=1-e^{-2I(\Delta T, \Delta S)}=\rho_{\Delta}^2=ICA}
#'
#' @references
#' Deliorman, G., Alonso, A., & Pardo, M. C. (2025).
#' A Rényi-Divergence-Based Family of Metrics for the Evaluation of Surrogate Endpoints
#' in a Causal Inference Framework. *Statistics in Biopharmaceutical Research*.
#' \url{https://doi.org/10.1080/19466315.2025.2484009}
#'
#' @examples
#'
#'ICA_alpha_ContCont(alpha=0.5, T0S0 = 0.9597334, T1S1=0.9644139, T0T0=544.3285, T1T1=550.6597, S0S0=180.6831, S1S1=180.9433, T0T1=-0.9,
#' T0S1=-0.9 , T1S0=-0.9 , S0S1=-0.9)
#'
#'More than one alpha value can be calculated simultaneously:
#'
#'ICA_alpha_ContCont(alpha=c(1.25,1.5,2), T0S0 = 0.9597334, T1S1=0.9644139, T0T0=544.3285, T1T1=550.6597, S0S0=180.6831, S1S1=180.9433, T0T1=-0.9,
#'T0S1=-0.9 , T1S0=-0.9 , S0S1=-0.9)
#'

#' @export
ICA_alpha_ContCont <- function(alpha=numeric(), T0S0, T1S1, T0T0=1, T1T1=1, S0S0=1, S1S1=1,
                               T0T1=seq(-1, 1, by=.1), T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1),
                               S0S1=seq(-1, 1, by=.1)) {

  T0S0_hier <- T0S0[1]
  T1S1_hier <- T1S1[1]
  alpha <- as.numeric(alpha)
  length_alpha<- length(alpha)
  Results <- na.exclude(matrix(NA, 1, (10+length_alpha)))
  colnames_result <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "sigma.delta.T", "sigma.delta.S", "rho", "ICA")
  ica_alpha_cols <- paste0("ICA_", as.character(alpha))
  colnames(Results) <- c(colnames_result, ica_alpha_cols)

  combins <- expand.grid(T0T1, T0S0_hier, T0S1, T1S0, T1S1_hier, S0S1)
  lengte <- dim(combins)[1]

  if (length(T0S0)>1){
    if (length(T0S0)<lengte){stop("The specified vector for T0S0 should be larger than ", lengte) }
    T0S0_vector <- T0S0[1:lengte] # sample
    combins[,2] <- T0S0_vector
  }
  if (length(T1S1)>1){
    if (length(T1S1)<lengte){stop("The specified vector for T1S1 should be larger than ", lengte) }
    T1S1_vector <- T1S1[1:lengte] # sample
    combins[,5] <- T1S1_vector
  }


  for (i in 1: nrow(combins)) {
    T0T1 <- combins[i, 1]
    T0S0 <- combins[i, 2]
    T0S1 <- combins[i, 3]
    T1S0 <- combins[i, 4]
    T1S1 <- combins[i, 5]
    S0S1 <- combins[i, 6]
    Sigma_c <- diag(4)
    Sigma_c[2,1] <- Sigma_c[1,2] <- T0T1 * (sqrt(T0T0)*sqrt(T1T1))
    Sigma_c[3,1] <- Sigma_c[1,3] <- T0S0 * (sqrt(T0T0)*sqrt(S0S0))
    Sigma_c[4,1] <- Sigma_c[1,4] <- T0S1 * (sqrt(T0T0)*sqrt(S1S1))
    Sigma_c[3,2] <- Sigma_c[2,3] <- T1S0 * (sqrt(T1T1)*sqrt(S0S0))
    Sigma_c[4,2] <- Sigma_c[2,4] <- T1S1 * (sqrt(T1T1)*sqrt(S1S1))
    Sigma_c[4,3] <- Sigma_c[3,4] <- S0S1 * (sqrt(S0S0)*sqrt(S1S1))
    Sigma_c[1,1] <- T0T0
    Sigma_c[2,2] <- T1T1
    Sigma_c[3,3] <- S0S0
    Sigma_c[4,4] <- S1S1
    Cor_c <- cov2cor(Sigma_c)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)

    if (Min.Eigen.Cor > 0) {
      rho <- ((sqrt(S0S0*T0T0)*Cor_c[3,1])+(sqrt(S1S1*T1T1)*Cor_c[4,2])-(sqrt(S0S0*T1T1)*Cor_c[3,2])-(sqrt(S1S1*T0T0)*Cor_c[4,1]))/(sqrt((T0T0+T1T1-(2*sqrt(T0T0*T1T1)*Cor_c[2,1]))*(S0S0+S1S1-(2*sqrt(S0S0*S1S1)*Cor_c[4,3]))))
      ICA<- rho^2
      if ((is.finite(ICA))==TRUE){
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * Cor_c[2,1])
        sigma.delta.S <- S0S0 + S1S1 - (2 * sqrt(S0S0*S1S1) * Cor_c[3,4])

        ICA_alpha_fun<- function(a) {1-(1-rho^2)*(1-(1-a)^2*rho^2)^(-1/(1-a)) }
        ICA_alphas <- data.frame()  # Initialize an empty data frame
        for (j in alpha) {
          if (j > 1 + (1 / abs(rho))) {
            stop("one of the alpha values does not meet the specific condition")
            rownames(Results) <- NULL # Warning message
            next  # Skip the current iteration
          }
          else if (j == 1) {
            ICA_alpha <- ICA  # Assign ICA directly if alpha is 1
          } else {
            ICA_alpha <- ICA_alpha_fun(j)  # Otherwise, call the function
          }
          ICA_alphas <- rbind(ICA_alphas, data.frame(ICA_alpha = ICA_alpha))
        }

        ICA_alphas<- ICA_alphas[1:nrow(ICA_alphas),]
        ICA_alphas<- as.numeric(ICA_alphas)
        results.part_1 <- as.vector(cbind(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1,  sigma.delta.T, sigma.delta.S, rho, ICA))
        results.part_2<- as.vector(ICA_alphas)
        results.part<- c(results.part_1,results.part_2)

        if (!exists("Results")) {
          Results <- data.frame(matrix(ncol = length(results.part), nrow = 0))
        }

        Results <- rbind(Results, results.part)
      }

      # Results <- rbind(Results, results.part)
      # rownames(Results) <- NULL}
    }
  }
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  Total.Num.Matrices <- nrow(combins)

  fit <-
    list(Total.Num.Matrices=Total.Num.Matrices, Pos.Def=Results[,1:6], rho=Results$rho, Sigmas=Results[,7:8], ICA=Results$ICA, ICA_alpha=Results[,11:length(results.part)], Call=match.call())

  class(fit) <- "ICA_alpha_ContCont"
  fit

}


