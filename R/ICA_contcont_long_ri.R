#' Assess surrogacy in the information-theoretic causal-inference framework (Individual Causal Association, ICA)
#' for a continuous surrogate and true endpoint measured repeatedly over time in a single-trial setting
#'
#' This function quantifies surrogacy under a linear model with subject-specific
#' random intercepts and separable residual covariance structure. See Details below.
#'
#' @param p (numeric) number of time points
#' @param AR1 (logical) should be `TRUE`. Other correlation structures for `R`
#'   can be added later.
#' @param AR1_rho (numeric) AR(1) correlation parameter for \eqn{R}. It must be
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
#' @param DT0S0 A scalar that specifies the covariance between the surrogate and
#'   the true endpoint in the control treatment condition for the `D` matrix.
#' @param DT1S1 A scalar that specifies the covariance between the surrogate and
#'   the true endpoint in the experimental treatment condition for the `D` matrix.
#' @param DT0T0 A scalar that specifies the variance of the true endpoint in the
#'   control treatment condition for the `D` matrix.
#' @param DT1T1 A scalar that specifies the variance of the true endpoint in the
#'   experimental treatment condition for the `D` matrix.
#' @param DS0S0 A scalar that specifies the variance of the surrogate endpoint in
#'   the control treatment condition for the `D` matrix.
#' @param DS1S1 A scalar that specifies the variance of the surrogate endpoint in
#'   the experimental treatment condition for the `D` matrix.
#' @param DT0T1 A scalar or vector that contains the correlation(s) between the
#'   counterfactuals `T0` and `T1` for the `D` matrix.
#' @param DT0S1 A scalar or vector that contains the correlation(s) between the
#'   counterfactuals `T0` and `S1` for the `D` matrix.
#' @param DT1S0 A scalar or vector that contains the correlation(s) between the
#'   counterfactuals `T1` and `S0` for the `D` matrix.
#' @param DS0S1 A scalar or vector that contains the correlation(s) between the
#'   counterfactuals `S0` and `S1` for the `D` matrix.
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
#' # Random intercepts with separable residual covariance structure:
#' The evolution of \eqn{y_j} over time is defined using a linear model with a
#' subject-specific random intercept for each potential outcome.
#' \deqn{T_{0j} = x_{T0j} \gamma_{T0} + b_{T0} + \varepsilon_{T0j}}
#' \deqn{T_{1j} = x_{T1j} \gamma_{T1} + b_{T1} + \varepsilon_{T1j}}
#' \deqn{S_{0j} = x_{S0j} \gamma_{S0} + b_{S0} + \varepsilon_{S0j}}
#' \deqn{S_{1j} = x_{S1j} \gamma_{S1} + b_{S1} + \varepsilon_{S1j}}
#'
#' It can be written as \eqn{y = x \gamma + b + \varepsilon}, where
#' \eqn{b \sim N(0, D)} with \eqn{D} an unstructured \eqn{4 \times 4}
#' covariance matrix and \eqn{\varepsilon \sim N(0, \Sigma)}. The marginal
#' model follows that \eqn{y \sim N(0, \Sigma_y)} where
#' \eqn{\Sigma_y = zDz^\top + \Sigma} with \eqn{\Sigma = R \otimes V}. Here
#' `V` is assumed to be unstructured and invariant across time.
#'
#' The marginal covariance matrix of the individual causal treatment effects is
#' \deqn{\Sigma_{\Delta} = \Sigma_u \otimes R + \Sigma_D \otimes 1_p 1_p^\top}
#' with \eqn{\Sigma_u = Q V Q^\top} and \eqn{\Sigma_D = Q D Q^\top}, where
#' \eqn{Q} maps \eqn{(T_0, T_1, S_0, S_1)} to \eqn{(\Delta T, \Delta S)}.
#'
#' Under the random-intercept model, the ICA metric is given by:
#' \deqn{R_{\Lambda}^{2} = 1 - (1 - \rho_u^2)^{p - 1}(1 - \rho_{ud}^2)}
#' where
#' \deqn{
#' \rho_u^2 = \frac{u_{12}^2}{u_{11} u_{22}},
#' \qquad
#' \rho_{ud}^2 = \frac{(u_{12} + \alpha d_{12})^2}
#' {(u_{11} + \alpha d_{11})(u_{22} + \alpha d_{22})}
#' }
#' and \eqn{\alpha = 1_p^\top R^{-1} 1_p}.
#'
#' If \eqn{R} is assumed to follow a first-order autoregressive structure, then
#' \deqn{\alpha = \frac{p - \rho p + 2\rho}{1 + \rho}}
#' where `rho = AR1_rho`.
#'
#'
#' @references
#' Deliorman, G., Pardo, M. C., Van der Elst, W., Molenberghs, G., & Alonso, A.
#' (submitted in 2026). Assessing Surrogate Endpoints in Longitudinal Studies:
#' An Information-Theoretic Causal Inference Approach.
#'
#' @examples
#' # Example output for the random-intercept model (symbolic 4x4 V matrix)
#' # V = | VT0T0  VT0T1  VT0S0  VT0S1 |
#' #     | VT0T1  VT1T1  VT1S0  VT1S1 |
#' #     | VT0S0  VT1S0  VS0S0  VS0S1 |
#' #     | VT0S1  VT1S1  VS0S1  VS1S1 |
#'
#' # Example output for the random-intercept model (symbolic 4x4 D matrix)
#' # D = | DT0T0  DT0T1  DT0S0  DT0S1 |
#' #     | DT0T1  DT1T1  DT1S0  DT1S1 |
#' #     | DT0S0  DT1S0  DS0S0  DS0S1 |
#' #     | DT0S1  DT1S1  DS0S1  DS1S1 |
#'
#' ICA_contcont_long_ri(
#'   p = 6,
#'   AR1 = TRUE,
#'   AR1_rho = 0.689,
#'   VT0S0 = 0.02377,
#'   VT1S1 = 0.03410,
#'   VT0T0 = 0.03685,
#'   VT1T1 = 0.04968,
#'   VS0S0 = 0.04581,
#'   VS1S1 = 0.06830,
#'   VT0T1 = seq(-1, 1, by=.2),
#'   VT0S1 = seq(-1, 1, by=.2),
#'   VT1S0 = seq(-1, 1, by=.2),
#'   VS0S1 = seq(-1, 1, by=.2),
#'   DT0S0 = 0.02131,
#'   DT1S1 = 0.02289,
#'   DT0T0 = 0.03301,
#'   DT1T1 = 0.02588,
#'   DS0S0 = 0.02801,
#'   DS1S1 = 0.02428,
#'   DT0T1 = seq(-1, 1, by=.2),
#'   DT0S1 = seq(-1, 1, by=.2),
#'   DT1S0 = seq(-1, 1, by=.2),
#'   DS0S1 = seq(-1, 1, by=.2)
#' )
#'
#'
#'
#' @export
ICA_contcont_long_ri <- function(p = numeric(),
                                 AR1 = TRUE,
                                 AR1_rho = numeric(),
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
                                 DT0S0,
                                 DT1S1,
                                 DT0T0,
                                 DT1T1,
                                 DS0S0,
                                 DS1S1,
                                 DT0T1 = seq(-1, 1, by = 0.1),
                                 DT0S1 = seq(-1, 1, by = 0.1),
                                 DT1S0 = seq(-1, 1, by = 0.1),
                                 DS0S1 = seq(-1, 1, by = 0.1)) {

  validate_scalar <- function(x, name, positive = FALSE) {
    if (length(x) != 1 || !is.finite(x)) {
      stop("`", name, "` must be a single finite numeric value.")
    }

    if (positive && x <= 0) {
      stop("`", name, "` must be strictly positive.")
    }
  }

  validate_correlation_vector <- function(x, name) {
    if (length(x) < 1 || any(!is.finite(x))) {
      stop("`", name, "` must contain at least one finite value.")
    }

    if (any(x < -1 | x > 1)) {
      stop("`", name, "` must contain values between -1 and 1.")
    }
  }

  validate_covariance_scalar <- function(cov_xy, var_x, var_y, name) {
    upper_bound <- sqrt(var_x * var_y)

    if (!is.finite(cov_xy) || abs(cov_xy) > upper_bound) {
      stop("`", name, "` must satisfy |cov| <= sqrt(var_x * var_y).")
    }
  }

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

  summarize_delta_components <- function(valid_candidates, T0T0, T1T1, S0S0, S1S1,
                                         T0S0_cov, T1S1_cov) {
    if (nrow(valid_candidates) == 0) {
      return(data.frame(
        c11 = numeric(0),
        c12 = numeric(0),
        c22 = numeric(0),
        rho2 = numeric(0)
      ))
    }

    cov_t0t1 <- valid_candidates$T0T1 * sqrt(T0T0 * T1T1)
    cov_t0s1 <- valid_candidates$T0S1 * sqrt(T0T0 * S1S1)
    cov_t1s0 <- valid_candidates$T1S0 * sqrt(T1T1 * S0S0)
    cov_s0s1 <- valid_candidates$S0S1 * sqrt(S0S0 * S1S1)

    c11 <- T0T0 + T1T1 - (2 * cov_t0t1)
    c12 <- T0S0_cov + T1S1_cov - cov_t0s1 - cov_t1s0
    c22 <- S0S0 + S1S1 - (2 * cov_s0s1)
    det_2x2 <- (c11 * c22) - (c12^2)

    keep <- is.finite(c11) & is.finite(c12) & is.finite(c22) &
      c11 > 0 & c22 > 0 & det_2x2 > 0

    if (!any(keep)) {
      return(data.frame(
        c11 = numeric(0),
        c12 = numeric(0),
        c22 = numeric(0),
        rho2 = numeric(0)
      ))
    }

    data.frame(
      c11 = c11[keep],
      c12 = c12[keep],
      c22 = c22[keep],
      rho2 = (c12[keep]^2) / (c11[keep] * c22[keep])
    )
  }

  enumerate_valid_matrices <- function(T0S0_cov, T1S1_cov, T0T0, T1T1, S0S0, S1S1,
                                       T0T1, T0S1, T1S0, S0S1) {
    combins <- expand.grid(
      T0T1 = T0T1,
      T0S0 = T0S0_cov,
      T0S1 = T0S1,
      T1S0 = T1S0,
      T1S1 = T1S1_cov,
      S0S1 = S0S1,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    if (nrow(combins) == 0) {
      return(list(total = 0, valid = combins[0, , drop = FALSE]))
    }

    keep <- logical(nrow(combins))

    for (i in seq_len(nrow(combins))) {
      Sigma_c <- build_covariance_matrix(
        T0T1 = combins$T0T1[i],
        T0S0 = combins$T0S0[i],
        T0S1 = combins$T0S1[i],
        T1S0 = combins$T1S0[i],
        T1S1 = combins$T1S1[i],
        S0S1 = combins$S0S1[i],
        T0T0 = T0T0,
        T1T1 = T1T1,
        S0S0 = S0S0,
        S1S1 = S1S1
      )

      if (is_positive_definite(Sigma_c)) {
        keep[i] <- TRUE
      }
    }

    list(total = nrow(combins), valid = combins[keep, , drop = FALSE])
  }

  p <- as.numeric(p)
  AR1_rho <- as.numeric(AR1_rho)

  if (!isTRUE(AR1)) {
    stop("Currently only the AR(1) structure is implemented, so `AR1` must be TRUE.")
  }

  if (length(p) != 1 || !is.finite(p) || p < 1 || p != round(p)) {
    stop("`p` must be a single positive integer.")
  }

  if (length(AR1_rho) != 1 || !is.finite(AR1_rho) || AR1_rho <= -1 || AR1_rho >= 1) {
    stop("`AR1_rho` must be a single finite value in (-1, 1).")
  }

  scalar_values <- list(
    VT0S0 = VT0S0,
    VT1S1 = VT1S1,
    VT0T0 = VT0T0,
    VT1T1 = VT1T1,
    VS0S0 = VS0S0,
    VS1S1 = VS1S1,
    DT0S0 = DT0S0,
    DT1S1 = DT1S1,
    DT0T0 = DT0T0,
    DT1T1 = DT1T1,
    DS0S0 = DS0S0,
    DS1S1 = DS1S1
  )

  variance_names <- c("VT0T0", "VT1T1", "VS0S0", "VS1S1", "DT0T0", "DT1T1", "DS0S0", "DS1S1")

  for (name in names(scalar_values)) {
    validate_scalar(scalar_values[[name]], name, positive = name %in% variance_names)
  }

  validate_covariance_scalar(VT0S0, VT0T0, VS0S0, "VT0S0")
  validate_covariance_scalar(VT1S1, VT1T1, VS1S1, "VT1S1")
  validate_covariance_scalar(DT0S0, DT0T0, DS0S0, "DT0S0")
  validate_covariance_scalar(DT1S1, DT1T1, DS1S1, "DT1S1")

  correlation_vectors <- list(
    VT0T1 = VT0T1,
    VT0S1 = VT0S1,
    VT1S0 = VT1S0,
    VS0S1 = VS0S1,
    DT0T1 = DT0T1,
    DT0S1 = DT0S1,
    DT1S0 = DT1S0,
    DS0S1 = DS0S1
  )

  for (name in names(correlation_vectors)) {
    validate_correlation_vector(correlation_vectors[[name]], name)
  }

  alpha <- (p - (p * AR1_rho) + (2 * AR1_rho)) / (1 + AR1_rho)

  v_candidates <- enumerate_valid_matrices(
    T0S0_cov = VT0S0,
    T1S1_cov = VT1S1,
    T0T0 = VT0T0,
    T1T1 = VT1T1,
    S0S0 = VS0S0,
    S1S1 = VS1S1,
    T0T1 = VT0T1,
    T0S1 = VT0S1,
    T1S0 = VT1S0,
    S0S1 = VS0S1
  )

  d_candidates <- enumerate_valid_matrices(
    T0S0_cov = DT0S0,
    T1S1_cov = DT1S1,
    T0T0 = DT0T0,
    T1T1 = DT1T1,
    S0S0 = DS0S0,
    S1S1 = DS1S1,
    T0T1 = DT0T1,
    T0S1 = DT0S1,
    T1S0 = DT1S0,
    S0S1 = DS0S1
  )

  num_pos_def_v <- nrow(v_candidates$valid)
  num_pos_def_d <- nrow(d_candidates$valid)
  v_components <- summarize_delta_components(
    valid_candidates = v_candidates$valid,
    T0T0 = VT0T0,
    T1T1 = VT1T1,
    S0S0 = VS0S0,
    S1S1 = VS1S1,
    T0S0_cov = VT0S0,
    T1S1_cov = VT1S1
  )

  d_components <- summarize_delta_components(
    valid_candidates = d_candidates$valid,
    T0T0 = DT0T0,
    T1T1 = DT1T1,
    S0S0 = DS0S0,
    S1S1 = DS1S1,
    T0S0_cov = DT0S0,
    T1S1_cov = DT1S1
  )

  max_pairs <- nrow(v_components) * nrow(d_components)
  r2_lambda <- numeric(max_pairs)
  out_index <- 0

  if (max_pairs > 0) {
    d11_all <- d_components$c11
    d12_all <- d_components$c12
    d22_all <- d_components$c22

    for (i in seq_len(nrow(v_components))) {
      u11 <- v_components$c11[i]
      u12 <- v_components$c12[i]
      u22 <- v_components$c22[i]
      rho_u2 <- v_components$rho2[i]

      denom_ud <- (u11 + (alpha * d11_all)) * (u22 + (alpha * d22_all))
      keep_pairs <- is.finite(denom_ud) & denom_ud > 0

      if (!any(keep_pairs)) {
        next
      }

      rho_ud2 <- ((u12 + (alpha * d12_all[keep_pairs]))^2) / denom_ud[keep_pairs]
      pair_values <- 1 - ((1 - rho_u2)^(p - 1) * (1 - rho_ud2))
      pair_values <- pair_values[is.finite(pair_values)]

      if (length(pair_values) == 0) {
        next
      }

      insert_at <- out_index + seq_along(pair_values)
      r2_lambda[insert_at] <- pair_values
      out_index <- out_index + length(pair_values)
    }
  }

  if (out_index == 0) {
    r2_lambda <- numeric(0)
  } else {
    r2_lambda <- r2_lambda[seq_len(out_index)]
  }

  num_pos_def_pairs <- length(r2_lambda)

  fit <- list(
    Num.Pos.Def.V = num_pos_def_v,
    Num.Pos.Def.D = num_pos_def_d,
    Num.Pos.Def.Pairs = num_pos_def_pairs,
    R2_Lambda = r2_lambda,
    Call = match.call()
  )

  class(fit) <- "ICA_contcont_long_ri"
  fit
}
