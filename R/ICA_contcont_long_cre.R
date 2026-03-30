#' Assess surrogacy in the information-theoretic causal-inference framework (Individual Causal Association, ICA)
#' for a continuous surrogate and true endpoint measured repeatedly over time in a single-trial setting
#'
#' This function quantifies surrogacy under a linear model with
#' common random-effects with separable residual covariance structure. See Details below.
#'
#' @param p (numeric) number of time points
#' @param time_points Numeric vector of measurement times. Its length must be  `p`.
#' @param R_structure Residual correlation structure for `R`. Supported
#'   values are `"AR1"` and `"CS"`.
#' @param R_rho (numeric) correlation parameter for `R`. It must be
#'   a single value in `(-1, 1)`.
#' @param VT0S0 A scalar that specifies the covariance between the surrogate and
#'   the true endpoint in the control treatment condition for the `V` matrix.
#' @param VT1S1 A scalar that specifies the covariance between the surrogate and
#'   the true endpoint in the experimental treatment condition for the `V` matrix.
#' @param VT0T0 A scalar that specifies the variance of the true endpoint in the
#'   control treatment condition for the `V` matrix.
#' @param VT1T1 A scalar that specifies the variance of the true endpoint in the
#'   experimental treatment condition for the `V` matrix.
#' @param VS0S0 A scalar that specifies the variance of the surrogate endpoint in
#'   the control treatment condition for the `V` matrix.
#' @param VS1S1 A scalar that specifies the variance of the surrogate endpoint in
#'   the experimental treatment condition for the `V` matrix.
#' @param VT0T1 A scalar or vector that contains the correlation(s) between the
#'   counterfactuals `T0` and `T1` for the `V` matrix.
#' @param VT0S1 A scalar or vector that contains the correlation(s) between the
#'   counterfactuals `T0` and `S1` for the `V` matrix.
#' @param VT1S0 A scalar or vector that contains the correlation(s) between the
#'   counterfactuals `T1` and `S0` for the `V` matrix.
#' @param VS0S1 A scalar or vector that contains the correlation(s) between the
#'   counterfactuals `S0` and `S1` for the `V` matrix.
#' @param D An `8 x 8` covariance-matrix for the random effects.
#'   Observed entries from the SAS output fit should be supplied on the covariance
#'   scale and unidentified entries should be coded as `NA`. The row/column
#'   order must be `(bT0, aT0, bT1, aT1, bS0, aS0, bS1, aS1)`, where `b`
#'   denotes the random intercept and `a` denotes the random slope.
#' @param N Sample size used by [ICA.ContCont.MultS.MPC()] when generating
#'   positive definite `D`-matrices.
#' @param M Number of Monte Carlo draws used by [ICA.ContCont.MultS.MPC()].
#' @param prob Optional probability vector passed to
#'   [ICA.ContCont.MultS.MPC()].
#' @param Seed Random seed passed to [ICA.ContCont.MultS.MPC()].
#' @param Show.Progress Logical; should [ICA.ContCont.MultS.MPC()] display
#'   progress information.
#'
#' @return
#'
#' * Num.Pos.Def.V: An object of class numeric that contains the number of
#'   positive definite `V` matrices.
#' * Num.Pos.Def.D: An object of class numeric that contains the number of
#'   positive definite `D` matrices.
#' * Num.Pos.Def.Pairs: An object of class numeric that contains the number of
#'   valid `V-D` pairs used in the final calculations.
#' * R2_Lambda: A scalar or vector that contains the individual causal association \eqn{R_{\Lambda}^2}.
#'
#' @details
#' # Common random-effects with separable residual covariance structure:
#' The evolution of \eqn{y_j} over time is defined using a linear mixed-effects model.
#' It is considered the important special case in which all four potential-outcome models share the
#' same random-effects structure, for instance a random intercept and slope.
#'
#' \deqn{T_{0j} = x_{T0j} \gamma_{T0} + z_{T0j}b_{T0} + \varepsilon_{T0j}}
#' \deqn{T_{1j} = x_{T1j} \gamma_{T1} + z_{T1j}b_{T1} + \varepsilon_{T1j}}
#' \deqn{S_{0j} = x_{S0j} \gamma_{S0} + z_{S0j}b_{S0} + \varepsilon_{S0j}}
#' \deqn{S_{1j} = x_{S1j} \gamma_{S1} + z_{S1j}b_{S1} + \varepsilon_{S1j}}
#'
#' In this setting, the row vectors of covariates associated with the random
#' effects satisfy \eqn{\mathbf{z}_{kzj} = \mathbf{z}_j} for
#' \eqn{k = T, S} and \eqn{z = 0, 1}, and we write
#' \eqn{\dim(\mathbf{z}_j) = q}. Moreover,
#' \eqn{\mathbf{Z}_{kz} = \mathbf{Z}} for \eqn{k = T, S}.
#'
#' Within this framework, the vector of individual causal treatment effects
#' \eqn{\mathbf{\Delta}} can be expressed as
#' \deqn{
#' \mathbf{\Delta} =
#' \mathbf{\mu}_{\Delta} + (\mathbf{I}_2 \otimes \mathbf{Z}) \,
#' \boldsymbol{\delta b} + \boldsymbol{\varepsilon}_{\Delta},
#' }
#' where \eqn{\boldsymbol{\delta b}} denotes the vector of random-effect
#' differences
#' \deqn{
#' \boldsymbol{\delta b} =
#' (\boldsymbol{\delta b}_T^{\top}, \boldsymbol{\delta b}_S^{\top})^{\top},
#' \qquad
#' \boldsymbol{\delta b}_k = \mathbf{b}_{k1} - \mathbf{b}_{k0},
#' \quad k = T, S.
#' }
#'
#' It follows that
#' \eqn{\boldsymbol{\delta b} \sim N(\mathbf{0}, \mathbf{D}_L)},
#' \eqn{\boldsymbol{\varepsilon}_{\Delta} \sim N(\mathbf{0},
#' \mathbf{\Sigma}_{\varepsilon_{\Delta}})}, and
#' \eqn{\mathbf{\Sigma}_{\varepsilon_{\Delta}} =
#' \mathbf{\Sigma}_u \otimes \mathbf{R}}. Therefore, the marginal covariance
#' matrix is given by
#' \deqn{
#' \mathbf{\Sigma}_{\Delta} =
#' \mathbf{\Sigma}_u \otimes \mathbf{R} +
#' (\mathbf{I}_2 \otimes \mathbf{Z}) \, \mathbf{D}_L \,
#' (\mathbf{I}_2 \otimes \mathbf{Z})^{\top},
#' }
#' where
#' \deqn{
#' \mathbf{D}_L =
#' \left(
#' \begin{array}{cc}
#' \mathbf{D}_{TT} & \mathbf{D}_{TS} \\
#' \mathbf{D}_{TS}^{\top} & \mathbf{D}_{SS}
#' \end{array}
#' \right),
#' }
#' with \eqn{\mathbf{D}_{TT} = \mathrm{Var}(\boldsymbol{\delta b}_T)},
#' \eqn{\mathbf{D}_{TS} = \mathrm{Cov}(\boldsymbol{\delta b}_T,
#' \boldsymbol{\delta b}_S)}, and
#' \eqn{\mathbf{D}_{SS} = \mathrm{Var}(\boldsymbol{\delta b}_S)}.
#'
#' Equivalently,
#' \deqn{
#' \mathbf{\Sigma}_{\Delta} =
#' \left(
#' \begin{array}{cc}
#' \mathbf{\Sigma}_{\Delta T \Delta T} &
#' \mathbf{\Sigma}_{\Delta T \Delta S} \\
#' \mathbf{\Sigma}_{\Delta T \Delta S}^{\top} &
#' \mathbf{\Sigma}_{\Delta S \Delta S}
#' \end{array}
#' \right),
#' }
#' where
#' \deqn{
#' \mathbf{\Sigma}_{\Delta T \Delta T} =
#' \mathbf{Z} \mathbf{D}_{TT} \mathbf{Z}^{\top} + u_{11}\mathbf{R},
#' }
#' \deqn{
#' \mathbf{\Sigma}_{\Delta T \Delta S} =
#' \mathbf{Z} \mathbf{D}_{TS} \mathbf{Z}^{\top} + u_{12}\mathbf{R},
#' }
#' \deqn{
#' \mathbf{\Sigma}_{\Delta S \Delta S} =
#' \mathbf{Z} \mathbf{D}_{SS} \mathbf{Z}^{\top} + u_{22}\mathbf{R}.
#' }
#'
#' The closed-form expression for \eqn{\Lambda} is
#' \deqn{
#' \Lambda =
#' (1 - \rho_u^2)^p
#' \frac{
#' \left| \mathbf{I}_{2q} +
#' \mathbf{D}_L (\mathbf{\Sigma}_u^{-1} \otimes \mathbf{E}) \right|
#' }{
#' \left| \mathbf{I}_q + u_{11}^{-1}\mathbf{D}_{TT}\mathbf{E} \right|
#' \left| \mathbf{I}_q + u_{22}^{-1}\mathbf{D}_{SS}\mathbf{E} \right|
#' },
#' }
#' where
#' \deqn{
#' \rho_u^2 =
#' 1 - \frac{|\mathbf{\Sigma}_u|}{u_{11}u_{22}} =
#' \frac{u_{12}^2}{u_{11}u_{22}},
#' \qquad
#' \mathbf{E} = \mathbf{Z}^{\top}\mathbf{R}^{-1}\mathbf{Z}
#' \in \mathbb{R}^{q \times q}.
#' }
#' The ICA metric is given by
#' \eqn{R_{\Lambda}^2 = 1 - \Lambda.}
#'
#' @references
#' Deliorman, G., Pardo, M. C., Van der Elst, W., Molenberghs, G., & Alonso, A.
#' (submitted in 2026). Assessing Surrogate Endpoints in Longitudinal Studies:
#' An Information-Theoretic Causal Inference Approach.
#'
#' @examples
#' # Example structure of V matrix for the random-effects model (4x4 V matrix)
#' # V = | VT0T0  VT0T1  VT0S0  VT0S1 |
#' #     | VT0T1  VT1T1  VT1S0  VT1S1 |
#' #     | VT0S0  VT1S0  VS0S0  VS0S1 |
#' #     | VT0S1  VT1S1  VS0S1  VS1S1 |
#'
#' # Example structure of D matrix for the random-effects model (8x8 D matrix)
#' # D = | d(bT0,bT0)  d(bT0,aT0)  ...  d(bT0,bS1)  d(bT0,aS1) |
#' #     | d(aT0,bT0)  d(aT0,aT0)  ...  d(aT0,bS1)  d(aT0,aS1) |
#' #     | ...                                           ...   |
#' #     | d(aS1,bT0)  d(aS1,aT0)  ...  d(aS1,bS1)  d(aS1,aS1) |
#'
#' # The SAS outputs are inserted as follows:
#'
#' D <- matrix(NA_real_, nrow = 8, ncol = 8)
#' un11 <- D[1, 1] <- 0.04162
#' un22 <- D[3, 3] <- 0.04140
#' un31 <- D[5, 1] <- D[1, 5] <- 0.02306
#' un33 <- D[5, 5] <- 0.03991
#' un42 <- D[7, 3] <- D[3, 7] <- 0.02896
#' un44 <- D[7, 7] <- 0.04811
#' un51 <- D[1, 2] <- D[2, 1] <- 0.001000
#' un53 <- D[2, 5] <- D[5, 2] <- 0.001224
#' un55 <- D[2, 2] <- 0.001086
#' un62 <- D[3, 4] <- D[4, 3] <- 0.000439
#' un64 <- D[4, 7] <- D[7, 4] <- 0.000025
#' un66 <- D[4, 4] <- 0.001452
#' un71 <- D[1, 6] <- D[6, 1] <- 0.001114
#' un73 <- D[5, 6] <- D[6, 5] <- 0.000374
#' un75 <- D[2, 6] <- D[6, 2] <- 0.000684
#' un77 <- D[6, 6] <- 0.000919
#' un82 <- D[3, 8] <- D[8, 3] <- 0.000918
#' un84 <- D[7, 8] <- D[8, 7] <- -0.00053
#' un86 <- D[4, 8] <- D[8, 4] <- 0.001472
#' un88 <- D[8, 8] <- 0.002100
#' \dontrun{
#' model <- ICA_contcont_long_cre(
#'   p = 6,
#'   time_points = c(0, 1, 2, 4, 6, 8),
#'   R_structure = "CS",
#'   R_rho = 0.000236,
#'   VT0S0 = 0.009203,
#'   VT1S1 = 0.008759,
#'   VT0T0 = 0.01365,
#'   VT1T1 = 0.01379,
#'   VS0S0 = 0.01819,
#'   VS1S1 = 0.01854,
#'   VT0T1 = seq(-1, 1, by = 0.2),
#'   VT0S1 = seq(-1, 1, by = 0.2),
#'   VT1S0 = seq(-1, 1, by = 0.2),
#'   VS0S1 = seq(-1, 1, by = 0.2),
#'   D = D,
#'   N = 568,
#'   M = 100,
#'   Seed = 123
#' )}
#'
#'
#'
#' @export
ICA_contcont_long_cre <- function(p,
                                  time_points,
                                  R_structure = c("AR1", "CS"),
                                  R_rho = numeric(),
                                  VT0S0,
                                  VT1S1,
                                  VT0T0,
                                  VT1T1,
                                  VS0S0,
                                  VS1S1,
                                  VT0T1 = seq(-1, 1, by = 0.1),
                                  VT0S1 = seq(-1, 1, by = 0.1),
                                  VT1S0 = seq(-1, 1, by = 0.1),
                                  VS0S1 = seq(-1, 1, by = 0.1),
                                  D,
                                  N,
                                  M = 10000,
                                  prob = NULL,
                                  Seed = 123,
                                  Show.Progress = FALSE) {
  build_covariance_matrix <- function(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1,
                                      T0T0, T1T1, S0S0, S1S1) {
    Sigma_c <- diag(c(T0T0, T1T1, S0S0, S1S1))

    Sigma_c[2, 1] <- Sigma_c[1, 2] <- T0T1 * sqrt(T0T0 * T1T1)
    Sigma_c[3, 1] <- Sigma_c[1, 3] <- T0S0
    Sigma_c[4, 1] <- Sigma_c[1, 4] <- T0S1 * sqrt(T0T0 * S1S1)
    Sigma_c[3, 2] <- Sigma_c[2, 3] <- T1S0 * sqrt(T1T1 * S0S0)
    Sigma_c[4, 2] <- Sigma_c[2, 4] <- T1S1
    Sigma_c[4, 3] <- Sigma_c[3, 4] <- S0S1 * sqrt(S0S0 * S1S1)

    Sigma_c
  }

  is_positive_definite <- function(mat) {
    chol_result <- try(chol(mat), silent = TRUE)
    !inherits(chol_result, "try-error")
  }

  build_residual_correlation <- function(p, structure, rho) {
    if (structure == "AR1") {
      return(fastmatrix::corAR1(rho, p))
    }

    fastmatrix::corCS(rho, p = p)
  }

  p <- as.numeric(p)
  R_structure <- match.arg(R_structure)
  R_rho <- as.numeric(R_rho)
  time_points <- as.numeric(time_points)
  p <- as.integer(p)

  if (length(time_points) != p) {
    stop("`time_points` must have length `p`.")
  }

  cor_matrix <- build_residual_correlation(p = p, structure = R_structure, rho = R_rho)
  v_grid <- expand.grid(
    T0T1 = VT0T1,
    T0S0 = VT0S0,
    T0S1 = VT0S1,
    T1S0 = VT1S0,
    T1S1 = VT1S1,
    S0S1 = VS0S1,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  complete_V <- list()
  if (nrow(v_grid) > 0) {
    for (i in seq_len(nrow(v_grid))) {
      V_i <- build_covariance_matrix(
        T0T1 = v_grid$T0T1[i],
        T0S0 = v_grid$T0S0[i],
        T0S1 = v_grid$T0S1[i],
        T1S0 = v_grid$T1S0[i],
        T1S1 = v_grid$T1S1[i],
        S0S1 = v_grid$S0S1[i],
        T0T0 = VT0T0,
        T1T1 = VT1T1,
        S0S0 = VS0S0,
        S1S1 = VS1S1
      )

      if (is_positive_definite(V_i)) {
        complete_V[[length(complete_V) + 1]] <- V_i
      }
    }
  }

  D_corr <- cov2cor(D)
  lower_template <- D_corr[lower.tri(D_corr)]
  missing_lower <- which(is.na(lower_template))
  D_sd <- sqrt(diag(D))
  complete_D <- list()

  if (length(missing_lower) == 0) {
    D_matrix <- D
    D_matrix[upper.tri(D_matrix)] <- t(D_matrix)[upper.tri(D_matrix)]
    if (is_positive_definite(D_matrix)) {
      complete_D[[1]] <- D_matrix
    }
  } else {
    ica_mpc <- suppressWarnings(
      ICA.ContCont.MultS.MPC(
        M = M,
        N = N,
        Sigma = D,
        prob = prob,
        Seed = Seed,
        Save.Corr = TRUE,
        Show.Progress = Show.Progress
      )
    )

    valid_rows <- which(complete.cases(ica_mpc$Lower.Dig.Corrs.All[, missing_lower, drop = FALSE]))

    for (i in valid_rows) {
      D_corr_i <- D_corr
      lower_values <- D_corr_i[lower.tri(D_corr_i)]
      lower_values[missing_lower] <- as.numeric(ica_mpc$Lower.Dig.Corrs.All[i, missing_lower])
      D_corr_i[lower.tri(D_corr_i)] <- lower_values
      D_corr_i[upper.tri(D_corr_i)] <- t(D_corr_i)[upper.tri(D_corr_i)]
      diag(D_corr_i) <- 1

      D_matrix <- D_corr_i * outer(D_sd, D_sd)
      D_matrix[upper.tri(D_matrix)] <- t(D_matrix)[upper.tri(D_matrix)]

      if (is_positive_definite(D_matrix)) {
        complete_D[[length(complete_D) + 1]] <- D_matrix
      }
    }
  }

  num_pos_def_v <- length(complete_V)
  num_pos_def_d <- length(complete_D)
  r2_lambda <- numeric(0)
  Z0 <- cbind(1, time_points)
  Q <- matrix(c(-1, 0, 1, 0,
                0, -1, 0, 1), nrow = 2, ncol = 4)
  Z_full <- kronecker(diag(4), Z0)
  A <- kronecker(Q, diag(p)) %*% Z_full

  if (num_pos_def_v > 0 && num_pos_def_d > 0) {
    for (i in seq_len(num_pos_def_v)) {
      sigma_u <- Q %*% complete_V[[i]] %*% t(Q)

      for (j in seq_len(num_pos_def_d)) {
        sigma_delta <- A %*% complete_D[[j]] %*% t(A) + kronecker(sigma_u, cor_matrix)
        sigma_delta <- round(sigma_delta, digits = 6)
        sigma_tt <- sigma_delta[1:p, 1:p, drop = FALSE]
        sigma_ss <- sigma_delta[(p + 1):(2 * p), (p + 1):(2 * p), drop = FALSE]

        if (is_positive_definite(sigma_delta) && is_positive_definite(sigma_u)) {
          pair_value <- 1 - det(sigma_delta) / (det(sigma_tt) * det(sigma_ss))
          r2_lambda <- c(r2_lambda, pair_value)
        }
      }
    }
  }

  fit <- list(
    Num.Pos.Def.V = num_pos_def_v,
    Num.Pos.Def.D = num_pos_def_d,
    Num.Pos.Def.Pairs = length(r2_lambda),
    R2_Lambda = r2_lambda,
    Call = match.call()
  )

  class(fit) <- "ICA_contcont_long_cre"
  fit
}
