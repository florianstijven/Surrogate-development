#' Likelihood and AIC summary for a longitudinal fit
#'
#' This function returns basic likelihood-based summaries for a fitted
#' longitudinal surrogate model.
#'
#' @param object An object of class \code{surrogate_longitudinal_fit}.
#'
#' @return A data.frame containing the model name, number of parameters,
#' log-likelihood, negative log-likelihood, AIC, and convergence code.
#'
#' @seealso \code{\link{fit_longitudinal_surrogate_model}} for a complete example
#' of fitting a longitudinal surrogate model.
#'
#' @export
longitudinal_fit_stats <- function(object) {
  if (!inherits(object, "surrogate_longitudinal_fit")) {
    stop("`object` must inherit from 'surrogate_longitudinal_fit'.")
  }

  logLik <- -object$fit$objective
  k <- length(object$fit$par)

  data.frame(
    model = object$model_spec$name,
    npar = k,
    logLik = logLik,
    neglogLik = object$fit$objective,
    AIC = 2 * k - 2 * logLik,
    convergence = object$fit$convergence,
    stringsAsFactors = FALSE
  )
}


#' Matrix diagnostics for fitted covariance components
#'
#' This function returns numerical diagnostics for the fitted covariance
#' components of a longitudinal surrogate model.
#'
#' @param object An object of class \code{surrogate_longitudinal_fit}.
#'
#' @return A data.frame containing determinant, condition number, reciprocal
#' condition number, minimum eigenvalue, Cholesky status, and diagnostic status
#' for each fitted covariance matrix.
#'
#' @seealso \code{\link{fit_longitudinal_surrogate_model}} for a complete example
#' of fitting a longitudinal surrogate model.
#'
#' @export
longitudinal_matrix_diagnostics <- function(object) {
  if (!inherits(object, "surrogate_longitudinal_fit")) {
    stop("`object` must inherit from 'surrogate_longitudinal_fit'.")
  }

  pars <- extract_longitudinal_parameters(object)
  mats <- list(V0 = pars$V0, V1 = pars$V1)
  if (!is.null(pars$D0)) mats$D0 <- pars$D0
  if (!is.null(pars$D1)) mats$D1 <- pars$D1

  one <- function(M, nm) {
    M <- 0.5 * (as.matrix(M) + t(as.matrix(M)))
    eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    cond <- tryCatch(kappa(M, exact = TRUE), error = function(e) NA_real_)
    chol_ok <- !inherits(try(chol(M), silent = TRUE), "try-error")
    rcond <- if (is.finite(cond)) 1 / cond else NA_real_
    status <- if (!chol_ok) {
      "not PD"
    } else if (is.finite(rcond) && rcond < 1e-8) {
      "PD but numerically singular"
    } else if (is.finite(rcond) && rcond < 1e-4) {
      "PD but ill-conditioned"
    } else {
      "PD"
    }
    data.frame(
      model = object$model_spec$name,
      matrix = nm,
      dimension = paste0(nrow(M), "x", ncol(M)),
      determinant = det(M),
      cond_num = cond,
      rcond_num = rcond,
      min_eig = min(eig),
      chol_ok = chol_ok,
      status = status,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, Map(one, mats, names(mats)))
}


#' Extract fitted covariance parameters
#'
#' This function extracts the fitted residual covariance matrices, temporal
#' correlation parameter, and, when present, fitted random-effects covariance
#' matrices from a longitudinal surrogate model.
#'
#' @param object An object of class \code{surrogate_longitudinal_fit}.
#'
#' @return A list containing fitted covariance parameters, including
#' \code{rho}, \code{V0}, \code{V1}, and, when present, \code{D0} and \code{D1}.
#'
#' @seealso \code{\link{fit_longitudinal_surrogate_model}} for a complete example
#' of fitting a longitudinal surrogate model.
#'
#' @export
extract_longitudinal_parameters <- function(object) {
  if (!inherits(object, "surrogate_longitudinal_fit")) {
    stop("`object` must inherit from 'surrogate_longitudinal_fit'.")
  }

  .surr_unpack_theta(
    object$fit$par,
    n_beta = object$n_beta,
    model_spec = object$model_spec,
    residual_type = object$residual_type,
    n_visits = length(object$visit_times),
    endpoint_order = object$endpoint_order
  )
}
