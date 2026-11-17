## ============================================================================
## Evaluation-grid handling and covariance construction
## ============================================================================

.surr_model_key <- function(object) {
  model <- object[["model", exact = TRUE]]
  if (!is.null(model)) return(model)
  spec <- object$model_spec
  if (identical(spec$re_mode, "none")) return("galecki")

  if (identical(spec$re_mode, "endpoint_specific")) {
    x <- unlist(spec$re_struct, use.names = FALSE)
    if (length(x) == 4L && all(x == "ri")) return("endpoint_ri")
    if (length(x) == 4L && all(x == "rislope")) return("endpoint_ri_rs")
  }

  if (identical(spec$re_mode, "shared_by_arm")) {
    x <- unlist(spec$shared_re, use.names = FALSE)
    if (length(x) == 2L && all(x == "ri")) return("shared_ri")
    if (length(x) == 2L && all(x == "rislope")) return("shared_ri_rs")
  }
  "custom"
}


.surr_validate_time_vector <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) < 1L || any(!is.finite(x))) {
    stop(name, " must contain at least one finite numeric time.")
  }
  if (anyDuplicated(x)) stop(name, " must not contain duplicated times.")
  if (any(x <= 0)) {
    stop(name, " must contain strictly positive post-randomization times; ",
         "baseline at time 0 is used for fitting but is not part of the ICA.")
  }
  if (is.unsorted(x, strictly = TRUE)) {
    stop(name, " must be strictly increasing.")
  }
  x
}


#' Validate a longitudinal ICA evaluation grid
#'
#' @param object A fitted longitudinal surrogate model.
#' @param evaluation_times A common numeric grid, a list with elements
#'   \code{T} and \code{S}, or \code{"common_observed"}.
#' @param extrapolation Action when requested times lie outside the endpoint's
#'   observed range: \code{"warn"} (default), \code{"allow"}, or
#'   \code{"error"}.
#' @return A list containing endpoint-specific evaluation times.
#' @export
validate_longitudinal_evaluation_times <- function(
    object,
    evaluation_times,
    extrapolation = c("warn", "allow", "error")) {
  if (!inherits(object, "surrogate_longitudinal_fit")) {
    stop("object must inherit from 'surrogate_longitudinal_fit'.")
  }
  extrapolation <- match.arg(extrapolation)

  if (is.character(evaluation_times)) {
    if (!identical(evaluation_times, "common_observed")) {
      stop("The only character specification is ",
           "evaluation_times = 'common_observed'.")
    }
    obs <- object$observed_times
    if (is.null(obs) || is.null(obs$T) || is.null(obs$S)) {
      stop("The fitted object does not store endpoint-specific observed times.")
    }
    if (!isTRUE(all.equal(obs$T, obs$S, tolerance = 0))) {
      stop("The observed schedules differ between T and S. ",
           "Specify evaluation_times = list(T = ..., S = ...) explicitly.")
    }
    out <- list(T = obs$T, S = obs$S)
  } else if (is.numeric(evaluation_times)) {
    tt <- .surr_validate_time_vector(evaluation_times, "evaluation_times")
    out <- list(T = tt, S = tt)
  } else if (is.list(evaluation_times) &&
             all(c("T", "S") %in% names(evaluation_times))) {
    out <- list(
      T = .surr_validate_time_vector(evaluation_times$T, "evaluation_times$T"),
      S = .surr_validate_time_vector(evaluation_times$S, "evaluation_times$S")
    )
  } else {
    stop("evaluation_times must be a numeric vector, a named list with ",
         "elements T and S, or 'common_observed'.")
  }

  ranges <- object$observed_time_ranges
  if (!is.null(ranges)) {
    messages <- character(0)
    for (ep in c("T", "S")) {
      rr <- ranges[[ep]]
      outside <- out[[ep]][out[[ep]] < rr[1L] | out[[ep]] > rr[2L]]
      if (length(outside) > 0L) {
        messages <- c(
          messages,
          paste0(ep, ": ", paste(outside, collapse = ", "),
                 " outside [", rr[1L], ", ", rr[2L], "]")
        )
      }
    }
    if (length(messages) > 0L) {
      msg <- paste0(
        "Evaluation times outside the observed endpoint-specific range were ",
        "requested (", paste(messages, collapse = "; "),
        "). The fitted covariance model will be extrapolated."
      )
      if (extrapolation == "warn") warning(msg, call. = FALSE)
      if (extrapolation == "error") stop(msg, call. = FALSE)
    }
  }

  out$common <- identical(out$T, out$S)
  class(out) <- c("longitudinal_evaluation_times", "list")
  out
}


.surr_temporal_cross <- function(times_a, times_b, rho, type) {
  times_a <- as.numeric(times_a)
  times_b <- as.numeric(times_b)
  delta <- abs(outer(times_a, times_b, "-"))

  if (type == "independence") return(1 * (delta == 0))
  if (length(rho) != 1L || !is.finite(rho)) {
    stop("A finite temporal-correlation parameter is required.")
  }
  if (type == "exponential") return(rho ^ delta)
  if (type == "cs") {
    out <- matrix(rho, length(times_a), length(times_b))
    out[delta == 0] <- 1
    return(out)
  }
  stop("Unsupported residual-correlation type: ", type)
}


.surr_evaluation_correlations <- function(object, evaluation_times) {
  pars <- extract_longitudinal_parameters(object)
  list(
    TT = .surr_temporal_cross(
      evaluation_times$T, evaluation_times$T,
      pars$rho, object$residual_type
    ),
    SS = .surr_temporal_cross(
      evaluation_times$S, evaluation_times$S,
      pars$rho, object$residual_type
    ),
    TS = .surr_temporal_cross(
      evaluation_times$T, evaluation_times$S,
      pars$rho, object$residual_type
    )
  )
}


.surr_evaluation_design <- function(object, evaluation_times, arm) {
  ep <- c(rep("T", length(evaluation_times$T)),
          rep("S", length(evaluation_times$S)))
  tt <- c(evaluation_times$T, evaluation_times$S)
  center <- object$time_center
  if (is.null(center)) center <- 0
  di <- data.frame(
    Z = rep(as.integer(arm), length(tt)),
    endpoint = factor(ep, levels = c("S", "T")),
    ctime = tt - center
  )
  Z <- .surr_make_Zrand(di, object$model_spec, object$endpoint_order)
  list(
    Z = Z,
    index_T = seq_len(length(evaluation_times$T)),
    index_S = length(evaluation_times$T) + seq_len(length(evaluation_times$S))
  )
}


.surr_logdet_pd <- function(M, name = "matrix") {
  M <- 0.5 * (M + t(M))
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(ch)) stop(name, " is not positive definite.")
  2 * sum(log(diag(ch)))
}


#' Construct covariance blocks on an evaluation grid
#'
#' This low-level function constructs the covariance blocks of the true- and
#' surrogate-endpoint treatment-effect trajectories for one completed set of
#' unidentified covariance parameters. It is primarily intended for validation
#' and advanced analyses.
#'
#' @param object A fitted longitudinal surrogate model.
#' @param evaluation_times Valid evaluation times.
#' @param residual_delta A 2 by 2 covariance matrix for the residual
#'   treatment-effect components, ordered T,S.
#' @param random_delta Optional covariance contribution on the complete
#'   evaluation vector, ordered as all T times followed by all S times.
#' @return A list with \code{Sigma_TT}, \code{Sigma_SS}, and \code{Sigma_TS}.
#' @export
construct_longitudinal_ica_covariance <- function(
    object,
    evaluation_times,
    residual_delta,
    random_delta = NULL) {
  evaluation_times <- validate_longitudinal_evaluation_times(
    object, evaluation_times, extrapolation = "allow"
  )
  U <- .surr_as_TS_cov(residual_delta, "residual_delta")
  R <- .surr_evaluation_correlations(object, evaluation_times)
  pT <- length(evaluation_times$T)
  pS <- length(evaluation_times$S)

  Sigma_TT <- U["T", "T"] * R$TT
  Sigma_SS <- U["S", "S"] * R$SS
  Sigma_TS <- U["T", "S"] * R$TS

  if (!is.null(random_delta)) {
    random_delta <- as.matrix(random_delta)
    if (!all(dim(random_delta) == c(pT + pS, pT + pS))) {
      stop("random_delta has incompatible dimensions.")
    }
    iT <- seq_len(pT)
    iS <- pT + seq_len(pS)
    Sigma_TT <- Sigma_TT + random_delta[iT, iT, drop = FALSE]
    Sigma_SS <- Sigma_SS + random_delta[iS, iS, drop = FALSE]
    Sigma_TS <- Sigma_TS + random_delta[iT, iS, drop = FALSE]
  }

  list(
    Sigma_TT = Sigma_TT,
    Sigma_SS = Sigma_SS,
    Sigma_TS = Sigma_TS,
    evaluation_times = evaluation_times
  )
}
