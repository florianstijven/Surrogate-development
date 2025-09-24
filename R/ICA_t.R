#' The function [ICA_t()] is to evaluate surrogacy in the single-trial causal-inference framework.
#'
#' @param df (numeric) is degree of freedom \eqn{\nu}. The maximum value for \eqn{df} is 342. When \eqn{df} exceeds this threshold, the model behavior aligns with the Individual Causal Association (ICA) under the normal causal model.
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
#' * ICA_t: A scalar or vector that contains the individual causal association \eqn{ICA_{t}}
#' * Sigmas: A data.frame that contains the \eqn{\sigma_{\Delta T}} and \eqn{\sigma_{\Delta S}}
#'


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
          print(
            " The maximum value of df is 342. When it is greater than 342, it is the same as ICA in the normal causal model. "
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


#' @details
#' # t-causal model \eqn{ICA_{t}}
#'The definition of a p-dimensional multivariate t-distribution is based on the fact that if \eqn{\y \sim \Norm(\0, \bSigma)} and \eqn{u \sim \chi_{\nu}^2} with \eqn{\y \perp u$ then $\x= \bmu+\y/\sqrt{u+\nu}}
#'follows a multivariate t-distribution with density function:

#' \eqn{  f(\x|\bmu, \bSigma, \nu)=\frac{\Gamma[(\nu+p)/2]}{\Gamma(\nu/2)(\nu \pi)^{p/2} |\bSigma|^{1/2} } \Big [ 1+ \frac{1}{\nu}(\x-\bmu)^T \bSigma^{-1}(\x-\bmu) \Big ]^{-(\nu+p)/2}}
#'
#'The multivariate t-distribution has parameters \eqn{\bSigma, \bmu, \nu} and it is denoted as \eqn{\x \sim t_p(\bmu, \bSigma, \nu)}.

#' \cite{Arellano 2013} provided an expression for the mutual information between \eqn{\x_1$ and $\x_2}:
#' \eqn{I_t(\x_1,\x_2)=&I_N(\x_1,\x_2)+\zeta(\nu)\quad\mbox{with}\\[5mm]
#' \zeta(\nu)=&\log\left[\dfrac{\Gamma(\nu/2)\Gamma[(\nu+p_1+p_2)/2]}{\Gamma[(\nu+p_1)/2]\Gamma[(\nu+p_2)/2]}\right]+\dfrac{\nu+p_2}{2}\psi\left(\dfrac{\nu+p_2}{2}\right)+
#'  \dfrac{\nu+p_1}{2}\psi\left(\dfrac{\nu+p_1}{2}\right)-\nonumber\\[8mm]
#' &\dfrac{\nu+p_1+p_2}{2}\psi\left(\dfrac{\nu+p_1+p_2}{2}\right)-\dfrac{\nu}{2}\psi\left(\dfrac{\nu}{2}\right)\nonumber
#' \end{align}}

#'where \eqn{\psi(x)=d/dx[\Gamma(x)]} is the so-called digamma function and
#' \eqn{I_N(\x_1,\x_2)=-\dfrac{1}{2}\log\left(\dfrac{\left|\bSigma\right|}{\left|\bSigma_{11}\right|\left|\bSigma_{22}\right|}\right)}

#' \eqn{ICA_t} can be defined as using the Squared Information Correlation Coefficient (SICC):
#'\eqn{ R_{Ht}^2=1-e^{-2I_t(\Delta T, \Delta S)}=1-e^{-2I_N(\Delta T, \Delta S)-2\zeta(\nu)}=1-(1-\rho_{\Delta}^2)e^{-2\zeta(\nu)}=ICA_t}
#' where
#' \eqn{\zeta(\nu)=2\log\left[\dfrac{\Gamma(\nu/2)}{\Gamma[(1+\nu)/2]}\sqrt{\dfrac{\nu}{2}}\,\right]+\left(1+\nu\right)\psi\left(\dfrac{1+\nu}{2}\right)-\left(1+\nu\right)\psi\left(\dfrac{\nu}{2}\right)-\dfrac{2+\nu}{\nu}}
#'
#' Note that when \eqn{\nu} approaches infinity, \eqn{\zeta(\nu)} converges to zero and \eqn{ICA_N=\rho_{\Delta}^2=ICA_t}.
#'



