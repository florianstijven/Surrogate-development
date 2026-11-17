## ============================================================================
## Unified longitudinal individual causal association
## ============================================================================

.surr_sensitivity_value <- function(sensitivity, name, default) {
  x <- sensitivity[[name]]
  if (is.null(x)) default else x
}


.surr_galecki_ica <- function(U, R, evaluation_times, tol) {
  if (nrow(U) == 0L) return(numeric(0))

  if (isTRUE(evaluation_times$common)) {
    return(.surr_validate_R2(
      1 - (1 - U$rho2)^length(evaluation_times$T)
    ))
  }

  logdet_T_R <- .surr_logdet_pd(R$TT, "R_TT")
  logdet_S_R <- .surr_logdet_pd(R$SS, "R_SS")
  pT <- length(evaluation_times$T)
  pS <- length(evaluation_times$S)
  out <- numeric(nrow(U))

  for (i in seq_len(nrow(U))) {
    full <- rbind(
      cbind(U$c11[i] * R$TT, U$c12[i] * R$TS),
      cbind(U$c12[i] * t(R$TS), U$c22[i] * R$SS)
    )
    logdet_full <- .surr_logdet_pd(full, "completed residual covariance")
    logdet_T <- pT * log(U$c11[i]) + logdet_T_R
    logdet_S <- pS * log(U$c22[i]) + logdet_S_R
    out[i] <- 1 - exp(logdet_full - logdet_T - logdet_S)
  }
  .surr_validate_R2(out, tol)
}


.surr_endpoint_ri_ica <- function(U, DD, R, evaluation_times,
                                  tol, progress = FALSE) {
  nU <- nrow(U)
  nD <- nrow(DD)
  if (nU == 0L || nD == 0L) return(numeric(0))

  pT <- length(evaluation_times$T)
  pS <- length(evaluation_times$S)
  logdet_R_TT <- .surr_logdet_pd(R$TT, "R_TT")
  logdet_R_SS <- .surr_logdet_pd(R$SS, "R_SS")
  oneT <- rep(1, pT)
  oneS <- rep(1, pS)
  alpha_T <- as.numeric(crossprod(oneT, solve(R$TT, oneT)))
  alpha_S <- as.numeric(crossprod(oneS, solve(R$SS, oneS)))

  d11 <- DD$c11
  d12 <- DD$c12
  d22 <- DD$c22
  out <- numeric(nU * nD)

  for (i in seq_len(nU)) {
    u11 <- U$c11[i]
    u12 <- U$c12[i]
    u22 <- U$c22[i]

    Sigma0 <- rbind(
      cbind(u11 * R$TT, u12 * R$TS),
      cbind(u12 * t(R$TS), u22 * R$SS)
    )
    ch0 <- tryCatch(chol(Sigma0), error = function(e) NULL)
    if (is.null(ch0)) {
      stop("The completed residual treatment-effect covariance is not ",
           "positive definite on the evaluation grid.")
    }
    logdet_Sigma0 <- 2 * sum(log(diag(ch0)))
    Sigma0_inv <- chol2inv(ch0)

    Z <- matrix(0, pT + pS, 2L)
    Z[seq_len(pT), 1L] <- 1
    Z[pT + seq_len(pS), 2L] <- 1
    M <- crossprod(Z, Sigma0_inv %*% Z)
    M <- 0.5 * (M + t(M))
    m11 <- M[1L, 1L]
    m12 <- M[1L, 2L]
    m22 <- M[2L, 2L]

    det_update <- (1 + d11 * m11 + d12 * m12) *
      (1 + d12 * m12 + d22 * m22) -
      (d11 * m12 + d12 * m22) *
      (d12 * m11 + d22 * m12)

    T_update <- 1 + (d11 / u11) * alpha_T
    S_update <- 1 + (d22 / u22) * alpha_S
    keep <- is.finite(det_update) & det_update > tol &
      is.finite(T_update) & T_update > tol &
      is.finite(S_update) & S_update > tol

    vals <- rep(NA_real_, nD)
    if (any(keep)) {
      logdet_full <- logdet_Sigma0 + log(det_update[keep])
      logdet_T <- pT * log(u11) + logdet_R_TT + log(T_update[keep])
      logdet_S <- pS * log(u22) + logdet_R_SS + log(S_update[keep])
      vals[keep] <- 1 - exp(logdet_full - logdet_T - logdet_S)
    }
    if (any(!keep)) {
      stop("A non-positive determinant update was encountered after a ",
           "positive-definite parameter completion.")
    }
    rows <- (i - 1L) * nD + seq_len(nD)
    out[rows] <- .surr_validate_R2(vals, tol)

    if (isTRUE(progress) && (i %% 500L == 0L || i == nU)) {
      message("Processed residual completion ", i, " of ", nU)
    }
  }
  out
}


.surr_shared_ri_candidates <- function(D0, D1, grid, tol) {
  d0 <- as.numeric(D0[1L, 1L])
  d1 <- as.numeric(D1[1L, 1L])
  grid <- .surr_validate_grid(grid, "sensitivity$random_effects")
  cov01 <- grid * sqrt(d0 * d1)
  keep <- is.finite(cov01) & (d0 * d1 - cov01^2) > tol
  data.frame(
    correlation = grid[keep],
    d_delta = d0 + d1 - 2 * cov01[keep]
  )
}


.surr_shared_ri_ica <- function(U, DD, R, evaluation_times, tol,
                                progress = FALSE) {
  nU <- nrow(U)
  nD <- nrow(DD)
  if (nU == 0L || nD == 0L) return(numeric(0))

  pT <- length(evaluation_times$T)
  pS <- length(evaluation_times$S)
  oneT <- rep(1, pT)
  oneS <- rep(1, pS)
  oneFull <- rep(1, pT + pS)
  alpha_T <- as.numeric(crossprod(oneT, solve(R$TT, oneT)))
  alpha_S <- as.numeric(crossprod(oneS, solve(R$SS, oneS)))
  logdet_R_TT <- .surr_logdet_pd(R$TT, "R_TT")
  logdet_R_SS <- .surr_logdet_pd(R$SS, "R_SS")
  d <- DD$d_delta
  out <- numeric(nU * nD)

  for (i in seq_len(nU)) {
    u11 <- U$c11[i]
    u12 <- U$c12[i]
    u22 <- U$c22[i]
    Sigma0 <- rbind(
      cbind(u11 * R$TT, u12 * R$TS),
      cbind(u12 * t(R$TS), u22 * R$SS)
    )
    ch0 <- tryCatch(chol(Sigma0), error = function(e) NULL)
    if (is.null(ch0)) {
      stop("The completed residual treatment-effect covariance is not ",
           "positive definite on the evaluation grid.")
    }
    logdet0 <- 2 * sum(log(diag(ch0)))
    alpha_full <- as.numeric(crossprod(
      oneFull, chol2inv(ch0) %*% oneFull
    ))

    full_update <- 1 + d * alpha_full
    T_update <- 1 + (d / u11) * alpha_T
    S_update <- 1 + (d / u22) * alpha_S
    keep <- is.finite(full_update) & full_update > tol &
      is.finite(T_update) & T_update > tol &
      is.finite(S_update) & S_update > tol
    if (any(!keep)) {
      stop("A non-positive shared-random-intercept determinant update ",
           "was encountered.")
    }

    vals <- 1 - exp(
      logdet0 + log(full_update) -
        (pT * log(u11) + logdet_R_TT + log(T_update)) -
        (pS * log(u22) + logdet_R_SS + log(S_update))
    )
    rows <- (i - 1L) * nD + seq_len(nD)
    out[rows] <- .surr_validate_R2(vals, tol)

    if (isTRUE(progress) && (i %% 500L == 0L || i == nU)) {
      message("Processed residual completion ", i, " of ", nU)
    }
  }
  out
}


.surr_validate_random_completions <- function(D0, D1, completions,
                                               tol = 1e-12) {
  if (!is.list(completions) || length(completions) < 1L) {
    stop("Random-intercept-and-slope ICA requires ",
         "sensitivity$random_completions, a non-empty list of cross-arm ",
         "Cov(b0,b1) matrices.")
  }
  q0 <- nrow(D0)
  q1 <- nrow(D1)

  out <- vector("list", length(completions))
  keep <- logical(length(completions))
  for (i in seq_along(completions)) {
    C01 <- as.matrix(completions[[i]])
    if (!all(dim(C01) == c(q0, q1)) || any(!is.finite(C01))) next
    if (!is.null(rownames(C01)) && !is.null(rownames(D0)) &&
        all(rownames(D0) %in% rownames(C01))) {
      C01 <- C01[rownames(D0), , drop = FALSE]
    }
    if (!is.null(colnames(C01)) && !is.null(rownames(D1)) &&
        all(rownames(D1) %in% colnames(C01))) {
      C01 <- C01[, rownames(D1), drop = FALSE]
    }
    joint <- rbind(cbind(D0, C01), cbind(t(C01), D1))
    eig_min <- min(eigen(0.5 * (joint + t(joint)),
                         symmetric = TRUE, only.values = TRUE)$values)
    if (is.finite(eig_min) && eig_min > tol) {
      out[[i]] <- C01
      keep[i] <- TRUE
    }
  }
  out[keep]
}


.surr_random_delta_from_completion <- function(object, evaluation_times,
                                                D0, D1, C01) {
  design0 <- .surr_evaluation_design(object, evaluation_times, 0L)$Z
  design1 <- .surr_evaluation_design(object, evaluation_times, 1L)$Z
  if (is.null(design0) || is.null(design1)) {
    stop("Random-effect designs are missing for a random-effects model.")
  }
  design0 %*% D0 %*% t(design0) +
    design1 %*% D1 %*% t(design1) -
    design0 %*% C01 %*% t(design1) -
    design1 %*% t(C01) %*% t(design0)
}


.surr_general_random_ica <- function(object, U, completions, R,
                                     evaluation_times, tol,
                                     progress = FALSE) {
  pars <- extract_longitudinal_parameters(object)
  C <- .surr_validate_random_completions(
    pars$D0, pars$D1, completions, tol
  )
  if (length(C) == 0L) return(list(values = numeric(0), n_random = 0L))

  random_delta <- lapply(C, function(x) {
    .surr_random_delta_from_completion(
      object, evaluation_times, pars$D0, pars$D1, x
    )
  })

  out <- numeric(nrow(U) * length(random_delta))
  pT <- length(evaluation_times$T)
  pS <- length(evaluation_times$S)
  iT <- seq_len(pT)
  iS <- pT + seq_len(pS)

  for (i in seq_len(nrow(U))) {
    base <- rbind(
      cbind(U$c11[i] * R$TT, U$c12[i] * R$TS),
      cbind(U$c12[i] * t(R$TS), U$c22[i] * R$SS)
    )
    vals <- vapply(random_delta, function(D) {
      full <- base + D
      R2Lambda_from_covariance(
        full[iT, iT, drop = FALSE],
        full[iS, iS, drop = FALSE],
        full[iT, iS, drop = FALSE],
        tol = tol
      )
    }, numeric(1))
    rows <- (i - 1L) * length(random_delta) + seq_along(random_delta)
    out[rows] <- vals

    if (isTRUE(progress) && (i %% 100L == 0L || i == nrow(U))) {
      message("Processed residual completion ", i, " of ", nrow(U))
    }
  }
  list(values = out, n_random = length(random_delta))
}


#' Calculate the longitudinal ICA from a fitted model
#'
#' The evaluation grid is mandatory because it defines the treatment-effect
#' trajectories and therefore the longitudinal ICA estimand.
#'
#' @param object A fitted longitudinal surrogate model.
#' @param evaluation_times A common numeric grid, endpoint-specific T/S grids,
#'   or \code{"common_observed"}.
#' @param sensitivity A list. Element \code{residual} is the grid for the four
#'   unidentified residual cross-world correlations. For RI models,
#'   \code{random_effects} is the corresponding random-intercept grid. For
#'   RI/RS models, \code{random_completions} must be a list of cross-arm
#'   covariance matrices \code{Cov(b0,b1)}.
#' @param extrapolation Extrapolation action; default \code{"warn"}.
#' @param return_candidates Retain admissible completion components.
#' @param progress Display progress for large calculations.
#' @param tol Numerical tolerance.
#' @return An object of class \code{surrogate_longitudinal_ica}.
#' @export
longitudinal_ica <- function(
    object,
    evaluation_times,
    sensitivity = list(
      residual = seq(-1, 1, by = 0.1),
      random_effects = seq(-1, 1, by = 0.1)
    ),
    extrapolation = c("warn", "allow", "error"),
    return_candidates = FALSE,
    progress = FALSE,
    tol = 1e-12) {
  if (!inherits(object, "surrogate_longitudinal_fit")) {
    stop("object must inherit from 'surrogate_longitudinal_fit'.")
  }
  extrapolation <- match.arg(extrapolation)
  eval <- validate_longitudinal_evaluation_times(
    object, evaluation_times, extrapolation
  )
  model <- .surr_model_key(object)
  pars <- extract_longitudinal_parameters(object)
  R <- .surr_evaluation_correlations(object, eval)

  grid_V <- .surr_sensitivity_value(
    sensitivity, "residual", seq(-1, 1, by = 0.1)
  )
  Vc <- .surr_completion_components(
    pars$V0, pars$V1, grid_V, tol = tol,
    return_candidates = return_candidates
  )
  U <- Vc$components

  D_candidates <- NULL
  nD <- NA_integer_

  if (model == "galecki") {
    values <- .surr_galecki_ica(U, R, eval, tol)
  } else if (model == "endpoint_ri") {
    grid_D <- .surr_sensitivity_value(
      sensitivity, "random_effects", seq(-1, 1, by = 0.1)
    )
    Dc <- .surr_completion_components(
      pars$D0, pars$D1, grid_D, tol = tol,
      return_candidates = return_candidates
    )
    D_candidates <- Dc$components
    nD <- nrow(D_candidates)
    values <- .surr_endpoint_ri_ica(
      U, D_candidates, R, eval, tol, progress
    )
  } else if (model == "shared_ri") {
    grid_D <- .surr_sensitivity_value(
      sensitivity, "random_effects", seq(-1, 1, by = 0.1)
    )
    D_candidates <- .surr_shared_ri_candidates(
      pars$D0, pars$D1, grid_D, tol
    )
    nD <- nrow(D_candidates)
    values <- .surr_shared_ri_ica(
      U, D_candidates, R, eval, tol, progress
    )
  } else if (model %in% c("endpoint_ri_rs", "shared_ri_rs", "custom")) {
    general <- .surr_general_random_ica(
      object = object,
      U = U,
      completions = sensitivity$random_completions,
      R = R,
      evaluation_times = eval,
      tol = tol,
      progress = progress
    )
    values <- general$values
    nD <- general$n_random
  } else {
    stop("Unsupported fitted model: ", model)
  }

  out <- list(
    model = model,
    evaluation_times = eval,
    R2_Lambda = .surr_validate_R2(values, tol),
    Num.Pos.Def.V = nrow(U),
    Num.Pos.Def.D = nD,
    Num.Pos.Def.Pairs = length(values),
    V_candidates = if (return_candidates) U else NULL,
    D_candidates = if (return_candidates) D_candidates else NULL,
    sensitivity = sensitivity,
    fit = object,
    Call = match.call()
  )
  class(out) <- c("surrogate_longitudinal_ica", "list")
  out
}


#' Fit a longitudinal model and calculate its ICA
#'
#' This is the single-call user interface. Model fitting uses all available
#' observations at their actual measurement times. The mandatory evaluation
#' grid is used only after fitting to construct the causal covariance matrices.
#' Baseline responses must be rows in \code{data} with \code{time = 0}; they
#' contribute to fitting but are excluded from the ICA evaluation grid.
#'
#' @inheritParams fit_longitudinal_surrogate
#' @inheritParams longitudinal_ica
#' @return An object of class \code{surrogate_longitudinal_ica}.
#' @export
ICA_contcont_long <- function(
    data,
    model,
    formula,
    evaluation_times,
    id,
    treatment,
    time,
    endpoint,
    outcome,
    surrogate_label,
    true_label,
    treatment_levels = c(0, 1),
    covariates = character(0),
    residual_type = "exponential",
    sensitivity = list(
      residual = seq(-1, 1, by = 0.1),
      random_effects = seq(-1, 1, by = 0.1)
    ),
    extrapolation = c("warn", "allow", "error"),
    start_fit = NULL,
    control = list(rel.tol = 1e-10,
                   x.tol = 1e-10,
                   eval.max = 5000,
                   iter.max = 3000),
    return_candidates = FALSE,
    progress = FALSE,
    tol = 1e-12) {
  fit <- fit_longitudinal_surrogate(
    data = data,
    model = model,
    formula = formula,
    id = id,
    treatment = treatment,
    time = time,
    endpoint = endpoint,
    outcome = outcome,
    surrogate_label = surrogate_label,
    true_label = true_label,
    treatment_levels = treatment_levels,
    covariates = covariates,
    residual_type = residual_type,
    start_fit = start_fit,
    control = control
  )

  longitudinal_ica(
    object = fit,
    evaluation_times = evaluation_times,
    sensitivity = sensitivity,
    extrapolation = match.arg(extrapolation),
    return_candidates = return_candidates,
    progress = progress,
    tol = tol
  )
}


#' Evaluate one fitted model on several evaluation grids
#'
#' @param object A fitted longitudinal surrogate model.
#' @param evaluation_sets A named list of evaluation-time specifications.
#' @param sensitivity,extrapolation,return_candidates,progress,tol Arguments
#'   with the same meaning as in longitudinal_ica().
#' @return A named list of longitudinal ICA objects.
#' @export
longitudinal_ica_evaluation_sets <- function(
    object,
    evaluation_sets,
    sensitivity = list(
      residual = seq(-1, 1, by = 0.1),
      random_effects = seq(-1, 1, by = 0.1)
    ),
    extrapolation = c("warn", "allow", "error"),
    return_candidates = FALSE,
    progress = FALSE,
    tol = 1e-12) {
  if (!is.list(evaluation_sets) || length(evaluation_sets) < 1L) {
    stop("evaluation_sets must be a non-empty named list.")
  }
  if (is.null(names(evaluation_sets)) || any(!nzchar(names(evaluation_sets)))) {
    stop("Every evaluation set must have a non-empty name.")
  }
  extrapolation <- match.arg(extrapolation)
  model <- .surr_model_key(object)

  ## Reuse the expensive completion enumeration for the principal multi-grid
  ## use case rather than repeating it for every evaluation set.
  if (model == "endpoint_ri") {
    pars <- extract_longitudinal_parameters(object)
    grid_V <- .surr_sensitivity_value(
      sensitivity, "residual", seq(-1, 1, by = 0.1)
    )
    grid_D <- .surr_sensitivity_value(
      sensitivity, "random_effects", seq(-1, 1, by = 0.1)
    )
    U <- .surr_completion_components(
      pars$V0, pars$V1, grid_V, tol = tol,
      return_candidates = return_candidates
    )$components
    DD <- .surr_completion_components(
      pars$D0, pars$D1, grid_D, tol = tol,
      return_candidates = return_candidates
    )$components

    return(lapply(evaluation_sets, function(x) {
      eval <- validate_longitudinal_evaluation_times(
        object, x, extrapolation
      )
      R <- .surr_evaluation_correlations(object, eval)
      values <- .surr_endpoint_ri_ica(
        U, DD, R, eval, tol, progress
      )
      out <- list(
        model = model,
        evaluation_times = eval,
        R2_Lambda = values,
        Num.Pos.Def.V = nrow(U),
        Num.Pos.Def.D = nrow(DD),
        Num.Pos.Def.Pairs = length(values),
        V_candidates = if (return_candidates) U else NULL,
        D_candidates = if (return_candidates) DD else NULL,
        sensitivity = sensitivity,
        fit = object,
        Call = match.call()
      )
      class(out) <- c("surrogate_longitudinal_ica", "list")
      out
    }))
  }

  lapply(evaluation_sets, function(x) {
    longitudinal_ica(
      object = object,
      evaluation_times = x,
      sensitivity = sensitivity,
      extrapolation = extrapolation,
      return_candidates = return_candidates,
      progress = progress,
      tol = tol
    )
  })
}
