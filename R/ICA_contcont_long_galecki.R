#' Assess surrogacy in the information-theoretic causal-inference framework (Individual Causal Association, ICA)
#' for a continuous surrogate and true endpoint measured repeatedly over time in a single-trial setting
#'
#'
#' This function quantifies surrogacy under Galecki's model. It is a special
#' case without subject-specific random effects. See Details below.
#'
#' @param p (numeric) number of time points
#' @param T0S0 A scalar that specifies the variance of the between the surrogate and the true endpoint in the control treatment condition
#' @param T1S1 A scalar that specifies the variance of the between the surrogate and the true endpoint in the experimental treatment condition
#' @param T0T0 A scalar that specifies the variance of the true endpoint in the control treatment condition
#' @param T1T1 A scalar that specifies the variance of the true endpoint in the experimental treatment condition
#' @param S0S0 A scalar that specifies the variance of the surrogate endpoint in the control treatment condition
#' @param S1S1 A scalar that specifies the variance of the surroagte endpoint in the experimental treatment condition
#' @param T0T1 A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1
#' @param T0S1 A scalar or vector that contains the correlation(s) between the counterfactuals T0 and S1
#' @param T1S0 A scalar or vector that contains the correlation(s) between the counterfactuals T1 and S0
#' @param S0S1 A scalar or vector that contains the correlation(s) between the counterfactuals S0 and S1
#'
#' @return
#'
#' * Total.Num.Matrices: An object of class numeric that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.
#' * Pos.Def: A data.frame that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the \eqn{\rho_{\Delta}} values.
#' * rho_u: A scalar or vector that contains the individual causal association \eqn{\rho_u}
#' * rho_u2: A scalar or vector that contains the individual causal association \eqn{\rho_u^2}
#' * R2_Lambda: A scalar or vector that contains the individual causal association \eqn{R_{\Lambda}^2}
#' * Sigmas: A data.frame that contains the \eqn{\sigma_{\Delta T}} and \eqn{\sigma_{\Delta S}}
#'
##' @details
#' # Galecki's model:
#' The evolution of \eqn{y_j} over time is defined using a marginal model without subject-specific random effects.
#' \deqn{ T_{0j} = x_{T0j} \gamma_{T0} + \varepsilon_{T0j}}
#' \deqn{T_{1j} = x_{T1j} \gamma_{T1} + \varepsilon_{T1j}}
#' \deqn{S_{0j} = x_{S0j} \gamma_{S0} + \varepsilon_{S0j}}
#' \deqn{S_{1j} = x_{S1j} \gamma_{S1} + \varepsilon_{S1j}}
#'
#'It can be written as \eqn{y=x \gamma + \varepsilon}. It is assummed that \eqn{\varepsilon \sim N(0, \Sigma)},
#' or, equivalently, \eqn{y \sim N(x \gamma, \Sigma_y)}. Essentially,
#' Galecki (1994) assumed that the covariance matrix is
#' separable, i.e. it can be written as \eqn{\Sigma_y = \Sigma = R \otimes V} , where \eqn{V} is assumed to be unstructured and
#' invariant across time, and \eqn{R} reflects a general correlation matrix.
#'
#'
#' Under Galecki’s model,
#' the vector of individual causal treatment effects \eqn{\Delta} is:
#' \deqn{ \Delta = \mu_{\Delta}+ \varepsilon_{\Delta} }
#' and the marginal covariance matrix of the individual causal treatment effects is
#' \eqn{\Sigma_{\Delta} = \Sigma_u \otimes R}.
#'
#' Under Galecki’s model, the ICA is given by:
#' \deqn{ ICA_{\Lambda}^2 = 1 - (1-\rho_{u}^2)^p }
#'
#' In this model \eqn{\rho_u=\rho_{\Delta}=corr(\Delta T, \Delta S)}, which is the metric used in the cross-sectional setting.
#'
#' \deqn{
#'\rho_u =
#'  \frac{
#'    \sqrt{\sigma_{T0T0} \, \sigma_{S0S0}} \, \rho_{T0S0} +
#'      \sqrt{\sigma_{T1T1} \, \sigma_{S1S1}} \, \rho_{T1S1} -
#'      \sqrt{\sigma_{T1T1} \, \sigma_{S0S0}} \, \rho_{T1S0} -
#'      \sqrt{\sigma_{T0T0} \, \sigma_{S1S1}} \, \rho_{T0S1}
#'  }{
#'    \sqrt{
#'      (\sigma_{T0T0} + \sigma_{T1T1} - 2 \sqrt{\sigma_{T0T0} \, \sigma_{T1T1}} \, \rho_{T0T1})
#'      (\sigma_{S0S0} + \sigma_{S1S1} - 2 \sqrt{\sigma_{S0S0} \, \sigma_{S1S1}} \, \rho_{S0S1})
#'    }
#'  }
#'}
#'
#' @references
#' Deliorman, G., Pardo, M. C., Van der Elst W., Molenberghs G., & Alonso, A. (submitted in 2026).
#' Assessing Surrogate Endpoints in Longitudinal Studies: An Information-Theoretic
#' Causal Inference Approach.
#'
#'
#' @examples
#' # Example output for Galecki's model (symbolic 4x4 V matrix)
#' # V = | σT0T0  σT0T1  σT0S0  σT0S1 | = | 544.32  σT0T1  266.56  σT0S1 |
#' #     | σT0T1  σT1T1  σT1S0  σT1S1 |   | σT0T1  550.65  σT1S0  205.17 |
#' #     | σT0S0  σT1S0  σS0S0  σS0S1 |   | 266.56  σT1S0  180.68  σS0S1 |
#' #     | σT0S1  σT1S1  σS0S1  σS1S1 |   | σT0S1  205.17  σS0S1  180.94 |
#'
#' ICA_contcont_long_galecki(p=4,  T0S0 = 266.56, T1S1=205.17,
#'                   T0T0=544.3285, T1T1=550.6597,
#'                   S0S0=180.6831, S1S1=180.9433,
#'                   T0T1=seq(-1, 1, by=.2),  T0S1=seq(-1, 1, by=.2),
#'                   T1S0=seq(-1, 1, by=.2) , S0S1=seq(-1, 1, by=.2))
#'
#'
#'
#' @export
ICA_contcont_long_galecki <- function(p=numeric(), T0S0, T1S1, T0T0, T1T1, S0S0, S1S1,
                                      T0T1=seq(-1, 1, by=.1), T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1),
                                      S0S1=seq(-1, 1, by=.1)) {


  Sigma_c <- matrix(NA,nrow=4, ncol=4)

  # Fill covariances
  Sigma_c[4,2] <- Sigma_c[2,4] <- T1S1
  Sigma_c[3,1] <- Sigma_c[1,3] <- T0S0

  Sigma_c[1,1] <-  T0T0
  Sigma_c[2,2] <-  T1T1
  Sigma_c[3,3] <-  S0S0
  Sigma_c[4,4] <-  S1S1

  Cor_c <- cov2cor(Sigma_c)

  T0S0 <- Cor_c[1,3]
  T1S1 <- Cor_c[2,4]
  T0T0 <- 1
  T1T1 <- 1
  S0S0 <- 1
  S1S1 <- 1

  T0S0_hier <- T0S0[1]
  T1S1_hier <- T1S1[1]
  p <- as.numeric(p)
  Results <- na.exclude(matrix(NA, 1, 11))
  colnames_result <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "sigma.delta.T", "sigma.delta.S", "rho_u", "rho_u2", "R2_Lambda")
  colnames(Results) <- c(colnames_result)

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
      rho_u <- ((sqrt(S0S0*T0T0)*Cor_c[3,1])+(sqrt(S1S1*T1T1)*Cor_c[4,2])-(sqrt(S0S0*T1T1)*Cor_c[3,2])-(sqrt(S1S1*T0T0)*Cor_c[4,1]))/(sqrt((T0T0+T1T1-(2*sqrt(T0T0*T1T1)*Cor_c[2,1]))*(S0S0+S1S1-(2*sqrt(S0S0*S1S1)*Cor_c[4,3]))))
      rho_u2<- rho_u^2
      if ((is.finite(rho_u2))==TRUE){
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * Cor_c[2,1])
        sigma.delta.S <- S0S0 + S1S1 - (2 * sqrt(S0S0*S1S1) * Cor_c[3,4])

        R2_Lambda<-1 - (1-rho_u2)^p
        R2_Lambda<- as.numeric(R2_Lambda)
        results.part_1 <- as.vector(cbind(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1,  sigma.delta.T, sigma.delta.S, rho_u, rho_u2, R2_Lambda))

        if (!exists("Results")) {
          Results <- data.frame(matrix(ncol = length(results.part_1), nrow = 0))
        }

        Results <- rbind(Results, results.part_1)
      }

      # Results <- rbind(Results, results.part)
      # rownames(Results) <- NULL}
    }
  }
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  Total.Num.Matrices <- nrow(combins)

  fit <-
    list(Total.Num.Matrices=Total.Num.Matrices, Pos.Def=Results[,1:6], rho_u=Results$rho_u, Sigmas=Results[,7:8], rho_u2=Results$rho_u2, R2_Lambda=Results$R2_Lambda, Call=match.call())

  class(fit) <- "ICA_contcont_long_galecki"
  fit

}

