#' Assess surrogacy in the information-theoretic causal-inference framework (Individual Causal Association, ICA)
#' for a continuous surrogate and true endpoint measured repeatedly over time in a single-trial setting
#'
#' This function examines how the ICA changes as surrogate endpoint measurements are added
#' sequentially when the true-endpoint treatment-effect trajectory is fixed. See Details below.
#'
#' @param V0,V1 Treatment-specific residual 2 x 2 covariance matrices for T,S.
#' @param D0,D1 Treatment-specific random-intercept 2 x 2 covariance matrices for T,S.
#' @param times_T Clinical evaluation time points for T.
#' @param times_S Clinical evaluation time points for S.
#' @param correlation Temporal correlation structures such as exponential, compound_symmetry, and independence.
#' @param temporal_rho Temporal-correlation parameter. Ignored when \code{correlation = "independence"}
#' @param grid_V Sensitivity grid for unidentified residual correlations.
#' @param grid_D Sensitivity grid for unidentified random-intercept correlations.
#' @param tol Numerical tolerance used when checking positive definiteness and validating ICA values.
#' @param return_candidates Logical. If \code{TRUE}, the retained admissible covariance
#'  completion candidates for the residual and random-intercept components are returned.
#' @param progress Logical. If \code{TRUE}, progress messages are printed while
#'  processing residual covariance completion candidates.
#'
#' @return
#' * Num.Pos.Def.V: An object of class numeric that contains the number of
#'   positive definite `V` matrices.
#' * Num.Pos.Def.D: An object of class numeric that contains the number of
#'   positive definite `D` matrices.
#' * Num.Pos.Def.Pairs: An object of class numeric that contains the number of
#'   valid `V-D` pairs used in the final calculations.
#' * R2_Lambda: A scalar or vector that contains the individual causal association \eqn{R_{\Lambda}^2}.
#' * V_candidates: retained residual covariance candidates, if requested.
#' * D_candidates: retained random-intercept covariance candidates, if requested.
#'
#' @details
#' ICA for progressively enlarged surrogate grids under endpoint-specific RI:
#' It keeps the complete true endpoint treatment-effect trajectory fixed on `times_T` and examine how the
#' ICA changes as surrogate endpoint measurements are added sequentially for specified surrogate grid in `times_S`.
#' The objective is to determine how much information the first k surrogate endpoint treatment effects
#' provide about the complete true endpoint treatment-effect trajectory.
#'
#' Consider the vector \eqn{(\Delta T1,...,\Delta Tp, \Delta S1,...\Delta Sk )} for  \eqn{k=1,...,p}.
#' Its covariance matrix, denoted by \eqn{\Sigma_{\Delta,p}(k)} is obtained by retaining the corresponding rows
#' and columns of the original \eqn{2p \times 2p} covariance matrix. Thus, \eqn{\Sigma_{\Delta,p}(k)} is the \eqn{5+k \times 5+k}
#' principal submatrix associated with the complete PANSS trajectory and the first k surrogate treatment effects.
#' For any value of \eqn{k}, partition the selected covariance matrix as
#'
#' \deqn{
#' \Sigma_{\Delta,p}(k) =
#' \begin{pmatrix}
#' \Sigma_{\Delta_{TT,p} } & \Sigma_{\Delta_{TS,p}(k)} \\
#' \Sigma_{\Delta_{TS,p}(k)} & \Sigma_{\Delta_{SS,p}(k) }
#' \end{pmatrix}
#' }
#'
#' The ICA for covariance completion \eqn{p} and surrogate endpoint grid size \eqn{k} is
#' \deqn{ R_{\Lambda, n}^{2}=1- \frac{\Sigma_{\Delta,p}}{\Delta_{TT,p}\Delta_{SS,p}(k)} }
#'
#' @references
#' Deliorman, G., Pardo, M. C., Van der Elst W., Molenberghs G., & Alonso, A. (submitted in 2026).
#' Assessing Surrogate Endpoints in Longitudinal Studies: An Information-Theoretic
#' Causal Inference Approach.
#'
#'
#' @examples
#' #Example outputs
#' V0 <- matrix(c(0.02579325, 0.02514841,
#'               0.02514841, 0.02872256),
#'               nrow = 2, byrow = TRUE,
#'               dimnames = list(c("T", "S"), c("T", "S")))
#'
#'
#' V1 <- matrix(c(0.02907611, 0.02962312,
#'                0.02962312, 0.03400179),
#'                nrow = 2, byrow = TRUE,
#'                dimnames = list(c("T", "S"), c("T", "S")))
#'
#' D0 <- matrix(c(0.02587092, 0.02458997,
#'                0.02458997, 0.02405581),
#'                nrow = 2, byrow = TRUE,
#'                dimnames = list(c("T", "S"), c("T", "S")))
#'
#' D1 <- matrix(c(0.03649246, 0.03225987,
#'               0.03225987, 0.02897932),
#'               nrow = 2, byrow = TRUE,
#'               dimnames = list(c("T", "S"), c("T", "S")))
#'
#' times_S <- list(week1 = c(1),
#'                         week1_2 = c(1, 2),
#'                         week1_2_4 = c(1, 2, 4),
#'                         week1_2_4_6 = c(1, 2, 4, 6),
#'                         week1_2_4_6_8 = c(1, 2, 4, 6, 8))
#'
#'
#' # Example run
#' ica_nested<- ICA_contcont_long_ri_nested(V0 = V0,  V1 = V1,
#'                             D0 = D0,  D1 = D1,
#'                             times_T = c(1,2,4,6,8),
#'                             times_S = times_S,
#'                             correlation  = "exponential",
#'                             temporal_rho = 0.7117544,
#'                             grid_V = seq(-1, 1, by = 0.5),
#'                             grid_D = seq(-1, 1, by = 0.5),
#'                             return_candidates = FALSE,
#'                             tol = 1e-12)
#'
#' summary_longitudinal_ica(list(R2_Lambda = ica_nested$R2_Lambda[, "week1_2_4"]))
#'
#' @export
ICA_contcont_long_ri_nested <- function(V0,
                                        V1,
                                        D0,
                                        D1,
                                        times_T,
                                        times_S,
                                        correlation = c("exponential", "compound_symmetry", "independence"),
                                        temporal_rho = NULL,
                                        grid_V = seq(-1, 1, by = 0.1),
                                        grid_D = seq(-1, 1, by = 0.1),
                                        return_candidates = FALSE,
                                        tol = 1e-12,
                                        progress = FALSE) {
  times_T <- as.numeric(times_T)
  correlation <- match.arg(correlation)

  if (length(times_T) < 1L || any(!is.finite(times_T)) || anyDuplicated(times_T)) {
    stop("`times_T` must contain distinct finite numeric values.")
  }
  if (!is.list(times_S) || length(times_S) < 1L) stop("`times_S` must be a non-empty list of surrogate grids.")
  times_S <- lapply(times_S, as.numeric)
  if (any(vapply(times_S, function(x) length(x) < 1L || any(!is.finite(x)) || anyDuplicated(x), logical(1)))) {
    stop("Each element of `times_S` must contain distinct finite numeric values.")
  }
  if (correlation %in% c("exponential", "compound_symmetry")) {
    if (is.null(temporal_rho) ||
        length(temporal_rho) != 1L ||
        !is.finite(temporal_rho)) {
      stop("`temporal_rho` must be supplied as a single finite value when correlation is 'exponential' or 'compound_symmetry'.")
    }
  }

  if (correlation == "exponential") {
    if (temporal_rho <= 0 || temporal_rho >= 1) {
      stop("For exponential correlation, `temporal_rho` must lie in (0,1).")
    }
  }

  if (correlation == "compound_symmetry") {
    all_times <- unique(c(times_T, unlist(times_S)))
    lower <- if (length(all_times) <= 1L) -1 else -1 / (length(all_times) - 1L)
    if (temporal_rho <= lower || temporal_rho >= 1) {
      stop("For compound symmetry, `temporal_rho` must lie in (", signif(lower, 5), ", 1).")
    }
  }

  if (correlation == "independence") {
    temporal_rho <- NA_real_
  }

  Vc <- .surr_completion_components(V0, V1, grid_V, tol = tol, return_candidates = TRUE)
  Dc <- .surr_completion_components(D0, D1, grid_D, tol = tol, return_candidates = TRUE)
  U <- Vc$components
  DD <- Dc$components


  nV <- nrow(U); nD <- nrow(DD); n_pairs <- nV * nD

  d11 <- DD$c11; d12 <- DD$c12; d22 <- DD$c22
  pT <- length(times_T)

  make_cross_R <- function(t1, t2, correlation, temporal_rho) {
    if (correlation == "exponential") {
      return(temporal_rho ^ abs(outer(t1, t2, "-")))
    }

    if (correlation == "compound_symmetry") {
      out <- matrix(temporal_rho, length(t1), length(t2))
      out[outer(t1, t2, "==")] <- 1
      return(out)
    }

    if (correlation == "independence") {
      return(1 * outer(t1, t2, "=="))
    }
  }

  if (correlation == "exponential") {
    R_TT <- make_longitudinal_R(times_T, rho = temporal_rho, type = "exponential")
  }

  if (correlation == "compound_symmetry") {
    R_TT <- make_longitudinal_R(times_T, rho = temporal_rho, type = "cs")
  }

  if (correlation == "independence") {
    R_TT <- make_longitudinal_R(times_T, rho = NULL, type = "independence")
  }

  logdet_R_TT <- as.numeric(determinant(R_TT, logarithm = TRUE)$modulus)
  oneT <- rep(1, pT)
  alpha_T <- as.numeric(crossprod(oneT, solve(R_TT, oneT)))

  R2 <- matrix(NA_real_, n_pairs, length(times_S))
  colnames(R2) <- names(times_S)
  if (is.null(colnames(R2)) || any(!nzchar(colnames(R2)))) {
    colnames(R2) <- paste0("grid_", seq_along(times_S))
  }

  for (ii in seq_len(nV)) {
    u11 <- U$c11[ii]; u12 <- U$c12[ii]; u22 <- U$c22[ii]
    rows <- ((ii - 1L) * nD + 1L):(ii * nD)

    T_update <- 1 + (d11 / u11) * alpha_T
    if (any(!is.finite(T_update) | T_update <= tol)) {
      stop("A non-positive true-endpoint determinant update was encountered.")
    }
    logdet_T <- pT * log(u11) + logdet_R_TT + log(T_update)

    for (kk in seq_along(times_S)) {
      s_times <- times_S[[kk]]
      pS <- length(s_times)

      if (correlation == "exponential") {
        R_SS <- make_longitudinal_R(s_times, rho = temporal_rho, type = "exponential")
      }

      if (correlation == "compound_symmetry") {
        R_SS <- make_longitudinal_R(s_times, rho = temporal_rho, type = "cs")
      }

      if (correlation == "independence") {
        R_SS <- make_longitudinal_R(s_times, rho = NULL, type = "independence")
      }

      R_TS <- make_cross_R(times_T, s_times, correlation, temporal_rho)


      Sigma0 <- rbind(cbind(u11 * R_TT, u12 * R_TS),
                      cbind(u12 * t(R_TS), u22 * R_SS))
      ch0 <- tryCatch(chol(Sigma0), error = function(e) NULL)
      if (is.null(ch0)) stop("The residual treatment-effect covariance is not positive definite.")
      logdet_Sigma0 <- 2 * sum(log(diag(ch0)))
      Sigma0_inv <- chol2inv(ch0)

      Z <- matrix(0, pT + pS, 2L)
      Z[seq_len(pT), 1] <- 1
      Z[pT + seq_len(pS), 2] <- 1
      M <- crossprod(Z, Sigma0_inv %*% Z)
      M <- 0.5 * (M + t(M))
      m11 <- M[1,1]; m12 <- M[1,2]; m22 <- M[2,2]

      det_update <- (1 + d11 * m11 + d12 * m12) *
        (1 + d12 * m12 + d22 * m22) -
        (d11 * m12 + d12 * m22) *
        (d12 * m11 + d22 * m12)
      if (any(!is.finite(det_update) | det_update <= tol)) {
        stop("A non-positive full-covariance determinant update was encountered.")
      }

      oneS <- rep(1, pS)
      alpha_S <- as.numeric(crossprod(oneS, solve(R_SS, oneS)))
      S_update <- 1 + (d22 / u22) * alpha_S
      if (any(!is.finite(S_update) | S_update <= tol)) {
        stop("A non-positive surrogate determinant update was encountered.")
      }
      logdet_R_SS <- as.numeric(determinant(R_SS, logarithm = TRUE)$modulus)
      logdet_full <- logdet_Sigma0 + log(det_update)
      logdet_S <- pS * log(u22) + logdet_R_SS + log(S_update)

      vals <- 1 - exp(logdet_full - logdet_T - logdet_S)
      R2[rows, kk] <- .surr_validate_R2(vals)
    }

    if (isTRUE(progress) && (ii %% 500L == 0L || ii == nV)) {
      message("Processed V completion ", ii, " of ", nV)
    }
  }

  out <- list(
    Num.Pos.Def.V = nV,
    Num.Pos.Def.D = nD,
    Num.Pos.Def.Pairs = n_pairs,
    times_T = times_T,
    times_S = times_S,
    correlation = correlation,
    temporal_rho = temporal_rho,
    R2_Lambda = R2,
    V_candidates = if (return_candidates) U else NULL,
    D_candidates = if (return_candidates) DD else NULL,
    Call = match.call()
  )
  class(out) <- c("ICA_contcont_long_ri_nested", "surrogate_longitudinal_ica_nested")
  out
}
