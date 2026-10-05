## ============================================================================
## Surrogate package staging code: utilities for longitudinal ICA sensitivity
## ============================================================================

.surr_assert_2x2_cov <- function(M, name = "matrix") {
  M <- as.matrix(M)
  if (!all(dim(M) == c(2L, 2L))) stop("`", name, "` must be a 2 x 2 matrix.")
  if (any(!is.finite(M))) stop("`", name, "` contains non-finite values.")
  if (!isTRUE(all.equal(M, t(M), tolerance = 1e-10))) stop("`", name, "` must be symmetric.")
  if (inherits(try(chol(M), silent = TRUE), "try-error")) stop("`", name, "` must be positive definite.")
  M
}

.surr_as_TS_cov <- function(M, name = "matrix") {
  M <- .surr_assert_2x2_cov(M, name)
  rn <- rownames(M)
  cn <- colnames(M)

  if (!is.null(rn) && !is.null(cn)) {
    if (all(c("T", "S") %in% rn) && all(c("T", "S") %in% cn)) {
      return(M[c("T", "S"), c("T", "S"), drop = FALSE])
    }
    map_name <- function(x) {
      out <- rep(NA_character_, length(x))
      out[grepl("(^|_)T$", x)] <- "T"
      out[grepl("(^|_)S$", x)] <- "S"
      out
    }
    rr <- map_name(rn); cc <- map_name(cn)
    if (all(c("T", "S") %in% rr) && all(c("T", "S") %in% cc)) {
      i <- match(c("T", "S"), rr)
      j <- match(c("T", "S"), cc)
      M <- M[i, j, drop = FALSE]
      dimnames(M) <- list(c("T", "S"), c("T", "S"))
      return(M)
    }
  }

  stop("`", name, "` must have endpoint names identifying T and S, e.g. T/S or ri_T/ri_S.")
}

.surr_validate_grid <- function(grid, name = "grid") {
  grid <- as.numeric(grid)
  if (length(grid) < 1L || any(!is.finite(grid))) stop("`", name, "` must contain finite values.")
  if (any(grid < -1 | grid > 1)) stop("`", name, "` must lie in [-1,1].")
  grid
}

.surr_validate_R2 <- function(x, tol = 1e-10) {
  if (any(!is.finite(x))) stop("Non-finite ICA values were produced.")
  if (any(x < -tol | x > 1 + tol)) {
    stop("ICA values outside [0,1] were produced; inspect the covariance completion and numerical stability.")
  }
  pmin(1, pmax(0, x))
}

## Enumerate all admissible 4 x 4 covariance completions for potential outcomes
## ordered as (T0,T1,S0,S1), given treatment-specific observed 2 x 2
## covariance matrices V0 = Var(T0,S0) and V1 = Var(T1,S1).
.surr_completion_components <- function(M0, M1, grid, tol = 1e-12,
                                        return_candidates = TRUE) {
  M0 <- .surr_as_TS_cov(M0, "M0")
  M1 <- .surr_as_TS_cov(M1, "M1")
  grid <- .surr_validate_grid(grid)

  T0T0 <- M0["T", "T"]; S0S0 <- M0["S", "S"]; T0S0 <- M0["T", "S"]
  T1T1 <- M1["T", "T"]; S1S1 <- M1["S", "S"]; T1S1 <- M1["T", "S"]

  comb <- expand.grid(T0T1 = grid,
                      T0S1 = grid,
                      T1S0 = grid,
                      S0S1 = grid,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)

  cov_T0T1 <- comb$T0T1 * sqrt(T0T0 * T1T1)
  cov_T0S1 <- comb$T0S1 * sqrt(T0T0 * S1S1)
  cov_T1S0 <- comb$T1S0 * sqrt(T1T1 * S0S0)
  cov_S0S1 <- comb$S0S1 * sqrt(S0S0 * S1S1)

  ## Schur-complement PD check, partitioning the 4 x 4 matrix into T and S blocks.
  det_TT <- T0T0 * T1T1 - cov_T0T1^2
  good_TT <- is.finite(det_TT) & det_TT > tol

  cond_S0S0 <- rep(NA_real_, nrow(comb))
  cond_S0S1 <- rep(NA_real_, nrow(comb))
  cond_S1S1 <- rep(NA_real_, nrow(comb))

  idx <- which(good_TT)
  if (length(idx) > 0L) {
    den <- det_TT[idx]
    a <- cov_T0T1[idx]
    b01 <- cov_T0S1[idx]
    b10 <- cov_T1S0[idx]

    cond_S0S0[idx] <- S0S0 -
      (T1T1 * T0S0^2 - 2 * a * T0S0 * b10 + T0T0 * b10^2) / den

    cond_S0S1[idx] <- cov_S0S1[idx] -
      (T1T1 * T0S0 * b01 - a * T0S0 * T1S1 - a * b10 * b01 + T0T0 * b10 * T1S1) / den

    cond_S1S1[idx] <- S1S1 -
      (T1T1 * b01^2 - 2 * a * b01 * T1S1 + T0T0 * T1S1^2) / den
  }

  det_cond_SS <- cond_S0S0 * cond_S1S1 - cond_S0S1^2
  keep <- good_TT & is.finite(cond_S0S0) & is.finite(cond_S1S1) &
    is.finite(det_cond_SS) & cond_S0S0 > tol & det_cond_SS > tol

  comb <- comb[keep, , drop = FALSE]
  cov_T0T1 <- cov_T0T1[keep]
  cov_T0S1 <- cov_T0S1[keep]
  cov_T1S0 <- cov_T1S0[keep]
  cov_S0S1 <- cov_S0S1[keep]

  c11 <- T0T0 + T1T1 - 2 * cov_T0T1
  c12 <- T0S0 + T1S1 - cov_T0S1 - cov_T1S0
  c22 <- S0S0 + S1S1 - 2 * cov_S0S1
  det_delta <- c11 * c22 - c12^2
  keep_delta <- is.finite(c11) & is.finite(c12) & is.finite(c22) &
    c11 > tol & c22 > tol & det_delta > tol

  out <- data.frame(c11 = c11[keep_delta],
                    c12 = c12[keep_delta],
                    c22 = c22[keep_delta],
                    rho2 = c12[keep_delta]^2 / (c11[keep_delta] * c22[keep_delta]))

  if (return_candidates) {
    comb <- comb[keep_delta, , drop = FALSE]
    out <- cbind(comb,
                 cov_T0T1 = cov_T0T1[keep_delta],
                 cov_T0S1 = cov_T0S1[keep_delta],
                 cov_T1S0 = cov_T1S0[keep_delta],
                 cov_S0S1 = cov_S0S1[keep_delta],
                 out,
                 det_delta = det_delta[keep_delta])
  }

  list(total = length(grid)^4, components = out)
}


#' Compute ICA from covariance blocks
#'
#' @param Sigma_TT Covariance matrix of the true-endpoint treatment-effect vector.
#' @param Sigma_SS Covariance matrix of the surrogate treatment-effect vector.
#' @param Sigma_TS Cross-covariance matrix Cov(Delta T, Delta S).
#' @param tol Numerical tolerance used when validating the ICA value.
#'
#' @return The determinant-based longitudinal ICA.
#'
#' @export
R2Lambda_from_covariance <- function(Sigma_TT, Sigma_SS, Sigma_TS, tol = 1e-10) {
  Sigma_TT <- as.matrix(Sigma_TT)
  Sigma_SS <- as.matrix(Sigma_SS)
  Sigma_TS <- as.matrix(Sigma_TS)
  if (nrow(Sigma_TT) != ncol(Sigma_TT) || nrow(Sigma_SS) != ncol(Sigma_SS)) {
    stop("`Sigma_TT` and `Sigma_SS` must be square.")
  }
  if (!all(dim(Sigma_TS) == c(nrow(Sigma_TT), nrow(Sigma_SS)))) {
    stop("`Sigma_TS` has incompatible dimensions.")
  }
  Sigma <- rbind(cbind(Sigma_TT, Sigma_TS),
                 cbind(t(Sigma_TS), Sigma_SS))
  ch_full <- tryCatch(chol(Sigma), error = function(e) NULL)
  ch_T <- tryCatch(chol(Sigma_TT), error = function(e) NULL)
  ch_S <- tryCatch(chol(Sigma_SS), error = function(e) NULL)
  if (is.null(ch_full) || is.null(ch_T) || is.null(ch_S)) {
    stop("All covariance matrices must be positive definite.")
  }
  ld_full <- 2 * sum(log(diag(ch_full)))
  ld_T <- 2 * sum(log(diag(ch_T)))
  ld_S <- 2 * sum(log(diag(ch_S)))
  .surr_validate_R2(1 - exp(ld_full - ld_T - ld_S), tol = tol)
}

##  Construct a 2 x 2 endpoint covariance matrix in T,S ordering
as_TS_covariance <- function(var_T, var_S, cov_TS) {
  if (any(!is.finite(c(var_T, var_S, cov_TS))) || var_T <= 0 || var_S <= 0) {
    stop("Variances must be positive and all inputs finite.")
  }
  M <- matrix(c(var_T, cov_TS, cov_TS, var_S), 2, 2,
              dimnames = list(c("T", "S"), c("T", "S")))
  .surr_assert_2x2_cov(M)
}

#' Summary of a longitudinal ICA sensitivity result
#'
#' This function summarizes the distribution of longitudinal ICA values obtained
#' from a sensitivity analysis.
#'
#' @param object An object returned by a longitudinal ICA sensitivity function.
#' @param threshold Numeric threshold used to compute the fraction of ICA values
#' below the threshold.
#'
#' @return A data.frame with summary statistics for ICA.
#'
#' @export
summary_longitudinal_ica <- function(object, threshold = 0.95) {
  if (is.null(object$R2_Lambda)) stop("Object does not contain `R2_Lambda`.")
  x <- object$R2_Lambda
  qs <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
  data.frame(minimum = min(x),
             q1 = qs[1],
             median = qs[2],
             q3 = qs[3],
             maximum = max(x),
             IQR = qs[3] - qs[1],
             fraction_below_threshold = mean(x < threshold),
             threshold = threshold)
}

