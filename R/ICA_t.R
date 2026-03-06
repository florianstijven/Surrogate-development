#' ICA under the t-causal model
#'
#' This function evaluates surrogacy in the single-trial causal-inference framework under the t-causal model.
#'
#' @param df (numeric) is degree of freedom \eqn{\nu}. The maximum value for df is 342. When df exceeds this threshold, the model behavior aligns with the Individual Causal Association (ICA) under the normal causal model.
#' @param T0S0 A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the control treatment condition (\eqn{\rho_{T0S0}})
#' @param T1S1 A scalar or vector that specifies the correlation(s) between the surrogate and the true endpoint in the experimental treatment condition (\eqn{\rho_{T1S1}})
#' @param T0T0 A scalar that specifies the variance of the true endpoint in the control treatment condition
#' @param T1T1 A scalar that specifies the variance of the true endpoint in the experimental treatment condition
#' @param S0S0 A scalar that specifies the variance of the surrogate endpoint in the control treatment condition
#' @param S1S1 A scalar that specifies the variance of the surrogate endpoint in the experimental treatment condition
#' @param T0T1 A scalar or vector that contains the correlation(s) between the counterfactuals T0 and T1
#' @param T0S1 A scalar or vector that contains the correlation(s) between the counterfactuals T0 and S1
#' @param T1S0 A scalar or vector that contains the correlation(s) between the counterfactuals T1 and S0
#' @param S0S1 A scalar or vector that contains the correlation(s) between the counterfactuals S0 and S1



#' @return
#'
#' * Total.Num.Matrices: An object of class numeric that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.
#' * Pos.Def: A data.frame that contains the positive definite matrices that can be formed based on the user-specified correlations. These matrices are used to compute the vector of the \eqn{\rho_{\Delta}} values.
#' * rho: A scalar or vector that contains the individual causal association \eqn{\rho_{\Delta}}
#' * ICA: A scalar or vector that contains the individual causal association under the normal causal model \eqn{\rho_{\Delta}^2=ICA_N}
#' * ICA_t: A scalar or vector that contains the individual causal association under the t-causal model \eqn{ICA_t}
#' * Sigmas: A data.frame that contains the \eqn{\sigma_{\Delta T}} and \eqn{\sigma_{\Delta S}}
#'
#'
#' @details
#' The multivariate t-distribution has parameters \eqn{\boldsymbol{\Sigma}, \boldsymbol{\mu}, \nu} and is denoted as \eqn{\mathbf{x} \sim t_p(\boldsymbol{\mu}, \boldsymbol{\Sigma}, \nu)}.
#'
#' Mutual information between t-distributed random vectors \eqn{\mathbf{x}_1} and \eqn{\mathbf{x}_2} is given by \cite{Arellano-Valle et al. (2013)}:
#'
#' \eqn{I_t(\mathbf{x}_1,\mathbf{x}_2) = I_N(\mathbf{x}_1,\mathbf{x}_2) + \zeta(\nu)}
#'
#' The individual causal association under the t-causal model can be computed as:
#'
#' \eqn{ICA_t = 1 - e^{-2 I_t(\Delta T, \Delta S)}=1 - e^{-2 I_N(\Delta T, \Delta S)-2\zeta(\nu)}}
#'
#' \eqn{ICA_t = 1 - (1-ICA_N) e^{-2 \zeta(\nu)}}
#'
#'
#' Note that when \eqn{\nu} approaches infinity, \eqn{\zeta(\nu)} converges to zero and \eqn{ICA_t=ICA_N}.
#'
#' @references
#' Arellano-Valle RB, Conreras-Reyes JE, and Genton MG. (2013).
#' Shannon Entropy and Mutual Information for Multivariate Skew-Elliptical Distributions.
#' *Scandinavian Journal of Statistics*.
#' \url{https://doi.org/10.1111/sjos.12718}
#'
#'
#' Deliorman, G., Stijven, F., Van der Elst W., Pardo, M. C., & Alonso, A. (submitted in 2025, under review).
#' The Effects of Causal Model Misspecification on Surrogate Evaluation.
#'
#' @export
ICA_t <- function(df,
                  T0S0,
                  T1S1,
                  T0T0 = 1,
                  T1T1 = 1,
                  S0S0 = 1,
                  S1S1 = 1,
                  T0T1 = seq(-1, 1, by = .1),
                  T0S1 = seq(-1, 1, by = .1),
                  T1S0 = seq(-1, 1, by = .1),
                  S0S1 = seq(-1, 1, by = .1)) {
  # The pracma package should be installed. If note, raise an error.
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("The 'pracma' package is required but not installed. Please install it.")
  }

  T0S0_hier <- T0S0[1]
  T1S1_hier <- T1S1[1]

  Results <- na.exclude(matrix(NA, 1, 11))
  colnames(Results) <- c(
    "T0T1",
    "T0S0",
    "T0S1",
    "T1S0",
    "T1S1",
    "S0S1",
    "rho",
    "ICA",
    "ICA_t",
    "sigma.delta.T",
    "sigma.delta.S"
  )
  combins <- expand.grid(T0T1, T0S0_hier, T0S1, T1S0, T1S1_hier, S0S1)
  lengte <- dim(combins)[1]
  zeta_fun <- function(nu)
    2 * log((gamma(nu / 2) * gamma(1 / 2)) / (sqrt(pi) * gamma((1 + nu) / 2)) *
              sqrt(nu / 2)) - (nu / 2) * pracma::psi(nu / 2) + (nu + 1) * pracma::psi((1 +
                                                                                         nu) / 2) - ((nu + 2) / 2) * pracma::psi((nu + 2) / 2)

  if (length(T0S0) > 1) {
    if (length(T0S0) < lengte) {
      stop("The specified vector for T0S0 should be larger than ", lengte)
    }
    T0S0_vector <- T0S0[1:lengte] # sample
    combins[, 2] <- T0S0_vector
  }
  if (length(T1S1) > 1) {
    if (length(T1S1) < lengte) {
      stop("The specified vector for T1S1 should be larger than ", lengte)
    }
    T1S1_vector <- T1S1[1:lengte] # sample
    combins[, 5] <- T1S1_vector
  }


  for (i in 1:nrow(combins)) {
    T0T1 <- combins[i, 1]
    T0S0 <- combins[i, 2]
    T0S1 <- combins[i, 3]
    T1S0 <- combins[i, 4]
    T1S1 <- combins[i, 5]
    S0S1 <- combins[i, 6]
    Sigma_c <- diag(4)
    Sigma_c[2, 1] <- Sigma_c[1, 2] <- T0T1 * (sqrt(T0T0) * sqrt(T1T1))
    Sigma_c[3, 1] <- Sigma_c[1, 3] <- T0S0 * (sqrt(T0T0) * sqrt(S0S0))
    Sigma_c[4, 1] <- Sigma_c[1, 4] <- T0S1 * (sqrt(T0T0) * sqrt(S1S1))
    Sigma_c[3, 2] <- Sigma_c[2, 3] <- T1S0 * (sqrt(T1T1) * sqrt(S0S0))
    Sigma_c[4, 2] <- Sigma_c[2, 4] <- T1S1 * (sqrt(T1T1) * sqrt(S1S1))
    Sigma_c[4, 3] <- Sigma_c[3, 4] <- S0S1 * (sqrt(S0S0) * sqrt(S1S1))
    Sigma_c[1, 1] <- T0T0
    Sigma_c[2, 2] <- T1T1
    Sigma_c[3, 3] <- S0S0
    Sigma_c[4, 4] <- S1S1
    Cor_c <- cov2cor(Sigma_c)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)

    if (Min.Eigen.Cor > 0) {
      rho <- ((sqrt(S0S0 * T0T0) * Cor_c[3, 1]) + (sqrt(S1S1 * T1T1) * Cor_c[4, 2]) -
                (sqrt(S0S0 * T1T1) * Cor_c[3, 2]) - (sqrt(S1S1 * T0T0) * Cor_c[4, 1])) /
        (sqrt((T0T0 + T1T1 - (
          2 * sqrt(T0T0 * T1T1) * Cor_c[2, 1]
        )) * (S0S0 + S1S1 - (
          2 * sqrt(S0S0 * S1S1) * Cor_c[4, 3]
        ))))
      ICA <- rho^2
      if ((is.finite(ICA)) == TRUE) {
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0 * T1T1) * Cor_c[2, 1])
        sigma.delta.S <- S0S0 + S1S1 - (2 * sqrt(S0S0 * S1S1) * Cor_c[3, 4])
        ICA_t <- 1 - (1 - rho^2) * exp(-2 * zeta_fun(df))
        results.part <- as.vector(
          cbind(
            T0T1,
            T0S0,
            T0S1,
            T1S0,
            T1S1,
            S0S1,
            rho,
            ICA,
            ICA_t,
            sigma.delta.T,
            sigma.delta.S
          )
        )
        Results <- rbind(Results, results.part)
        if ((df > 343) == TRUE)  {
          stop(
            "Maximum df is 342. For larger df, use normal ICA"
          )
        }
        rownames(Results) <- NULL
      }
    }
  }
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  Total.Num.Matrices <- nrow(combins)

  fit <-
    list(
      Total.Num.Matrices = Total.Num.Matrices,
      Pos.Def = Results[, 1:6],
      rho = Results$rho,
      ICA = Results$ICA,
      ICA_t = Results$ICA_t,
      Sigmas = Results[, 10:11],
      Call = match.call()
    )

  class(fit) <- "ICA_t"
  fit

}
