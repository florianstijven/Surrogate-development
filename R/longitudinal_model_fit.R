#' Longitudinal observed-data model fitting
#'
#'
#' This function fits a bivariate longitudinal surrogate/true-endpoint model. General-purpose functions for two continuous longitudinal
#' endpoints measured in two treatment arms. Residual covariance is separable at the observed-data
#' level, with treatment-specific 2 x 2 endpoint covariance matrices and a common
#' temporal correlation matrix. Random effects may be endpoint-specific or shared
#' within treatment arm.
#'
#' @param data Standardized data returned by prepare_longitudinal_surrogate_data().
#' @param formula Fixed-effects formula using the standardized response `Y`.
#' @param model_spec Random-effects specification.  \code{\link{make_longitudinal_re_none}},
#' \code{\link{make_longitudinal_re_endpoint_specific}}, or \code{\link{make_longitudinal_re_shared_by_arm}}.
#' @param residual_type One of "exponential", "independence", or "cs".
#' @param start_fit Optional previously fitted surrogate_longitudinal_fit used to
#'        initialize common parameters.
#' @param control Control list passed to nlminb().
#'
#'
#' @return An object of class surrogate_longitudinal_fit.
#'
#'
#' @details
#' The data should first be standardized using
#' \code{\link{prepare_longitudinal_surrogate_data}}. The mean structure is
#' specified through \code{mean_formula}, while the random-effects structure is
#' specified through \code{model_spec}.
#'
#' The argument \code{model_spec} can be created using one of the helper
#' functions:
#' \itemize{
#'   \item \code{\link{make_longitudinal_re_none}} for Galecki's marginal model,
#'   without subject-specific random effects;
#'   \item \code{\link{make_longitudinal_re_endpoint_specific}} for
#'   endpoint-specific random-effects models;
#'   \item \code{\link{make_longitudinal_re_shared_by_arm}} for shared
#'   random-effects models by treatment arm.
#' }
#'
#' The residual temporal correlation structure is selected using
#' \code{residual_type}. The currently supported options are
#' \code{"exponential"}, \code{"cs"} denotes compound symmetry, and \code{"independence"}.
#' The exponential structure uses the actual temporal distances between
#' measurements; for equally spaced visits, it reduces to the usual AR(1)
#' structure.
#'
#'
#'
#' @references
#' Deliorman, G., Pardo, M. C., Van der Elst W., Molenberghs G., & Alonso, A. (submitted in 2026).
#' Assessing Surrogate Endpoints in Longitudinal Studies: An Information-Theoretic
#' Causal Inference Approach.
#'
#' @examples
#' \dontrun{
#' data(Schizo_PANSS_BPRS)
#'
#' dat <- prepare_longitudinal_surrogate_data(
#'        data = Schizo_PANSS_BPRS,
#'        id = "Id",
#'        treatment = "Treat",
#'        time = "Time_cont",
#'        endpoint = "Endpoint",
#'        outcome = "Response_log",
#'        baseline = "Baseline_log",
#'        surrogate_label = "BPRS",
#'        true_label = "PANSS",
#'        treatment_levels = c(0, 1))
#'
#' # endpoint- and treatment-specific baseline adjustment.
#'
#' mean_formula <- Y ~ 0 + endpoint:visitF:treatF + endpoint:treatF:Ybase
#'
#' spec_galecki <- make_longitudinal_re_none(name = "Galecki's marginal model")
#'
#' #do not run
#' fit_galecki <- fit_longitudinal_surrogate_model(dat,
#'                                                 mean_formula,
#'                                                 model_spec = spec_galecki,
#'                                                 residual_type = "exponential")
#'
#' longitudinal_fit_stats(fit_galecki)
#' extract_longitudinal_parameters(fit_galecki)
#' longitudinal_matrix_diagnostics(fit_galecki)}
#'
#' @export
fit_longitudinal_surrogate_model <- function(data,
                                             formula,
                                             model_spec,
                                             residual_type = c("exponential", "independence", "cs"),
                                             start_fit = NULL,
                                             control = list(rel.tol = 1e-10,
                                                            x.tol = 1e-10,
                                                            eval.max = 5000,
                                                            iter.max = 3000)) {
  residual_type <- match.arg(residual_type)
  needed <- c("id", "Z", "time_cont", "ctime", "endpoint", "Y")
  miss <- setdiff(needed, names(data))
  if (length(miss) > 0L) stop("Standardized data are missing: ", paste(miss, collapse = ", "))

  endpoint_order <- c("S", "T")
  if (!all(endpoint_order %in% levels(data$endpoint))) {
    stop("The standardized endpoint factor must contain levels 'S' and 'T'.")
  }

  Xfull <- model.matrix(formula, data = data)
  n_beta <- ncol(Xfull)
  if (qr(Xfull)$rank < n_beta) warning("The fixed-effects design matrix is rank deficient.")

  ols <- lm(formula, data = data, y = TRUE)
  beta0 <- coef(ols)
  beta0[is.na(beta0)] <- 0
  subjects <- split(data, data$id)
  visit_times <- sort(unique(data$time_cont))
  n_visits <- length(visit_times)

  make_D_start <- function(arm, previous = NULL) {
    nm <- .surr_random_effect_names(model_spec, arm, endpoint_order)
    q <- length(nm)
    if (q == 0L) return(numeric(0))

    s0 <- sd(data$Y, na.rm = TRUE)
    diag_vals <- ifelse(grepl("^rs_", nm) | grepl("shared_rs", nm),
                        0.01 * s0^2,
                        0.20 * s0^2)

    D_start <- diag(diag_vals, q)
    dimnames(D_start) <- list(nm, nm)

    if (!is.null(previous)) {
      prev_pars <- extract_longitudinal_parameters(previous)
      prev_D <- prev_pars[[paste0("D", arm)]]
      if (!is.null(prev_D)) {
        common <- intersect(nm, rownames(prev_D))
        if (length(common) > 0L) {
          D_start[common, common] <- prev_D[common, common, drop = FALSE]
        }
      }
    }
    .surr_safe_pack_chol(D_start)
  }

  make_start <- function(previous = NULL) {
    s0 <- sd(data$Y, na.rm = TRUE)
    if (is.null(previous)) {
      beta_start <- beta0
      rho_start <- if (residual_type == "independence") numeric(0) else
        .surr_eta_from_rho(0.6, residual_type, n_visits)
      V0_start <- .surr_pack_chol(diag(rep(s0^2, 2L)))
      V1_start <- .surr_pack_chol(diag(rep(s0^2, 2L)))
    } else {
      prev <- extract_longitudinal_parameters(previous)
      beta_start <- prev$beta
      if (length(beta_start) != n_beta) {
        stop("`start_fit` has an incompatible fixed-effects design.")
      }
      rho_start <- if (residual_type == "independence") numeric(0) else {
        rho0 <- if (is.finite(prev$rho)) prev$rho else 0.6
        bounds <- .surr_rho_bounds(residual_type, n_visits)
        rho0 <- min(max(rho0, bounds[1] + 1e-6), bounds[2] - 1e-6)
        .surr_eta_from_rho(rho0, residual_type, n_visits)
      }
      V0_start <- .surr_safe_pack_chol(prev$V0)
      V1_start <- .surr_safe_pack_chol(prev$V1)
    }
    c(beta_start, rho_start, V0_start, V1_start,
      make_D_start(0, previous), make_D_start(1, previous))
  }

  nll_subject <- function(di, pars) {
    Xi <- model.matrix(formula, data = di)
    mu <- drop(Xi %*% pars$beta)
    arm <- di$Z[1]
    V <- if (arm == 0L) pars$V0 else pars$V1

    Rtime <- make_longitudinal_R(di$time_cont,
                                 rho = if (residual_type == "independence") NULL else pars$rho,
                                 type = residual_type)
    Vend <- V[as.character(di$endpoint), as.character(di$endpoint), drop = FALSE]
    Sigmae <- Rtime * Vend

    Zb <- .surr_make_Zrand(di, model_spec, endpoint_order)
    if (is.null(Zb)) {
      Sigmab <- matrix(0, nrow(di), nrow(di))
    } else {
      D <- pars[[paste0("D", arm)]]
      Sigmab <- Zb %*% D %*% t(Zb)
    }
    Sigma <- Sigmae + Sigmab

    U <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(U)) return(1e12)
    res <- di$Y - mu
    tmp <- backsolve(U, res, transpose = TRUE)
    0.5 * (length(res) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(tmp^2))
  }

  objective <- function(theta) {
    pars <- .surr_unpack_theta(theta, n_beta, model_spec,
                               residual_type, n_visits, endpoint_order)
    sum(vapply(subjects, nll_subject, numeric(1), pars = pars))
  }

  start <- make_start(start_fit)
  opt <- nlminb(start = start, objective = objective, control = control)

  out <- list(
    fit = opt,
    formula = formula,
    n_beta = n_beta,
    beta_names = colnames(Xfull),
    model_spec = model_spec,
    residual_type = residual_type,
    endpoint_order = endpoint_order,
    visit_times = visit_times,
    n_subjects = length(subjects),
    n_observations = nrow(data),
    call = match.call()
  )
  class(out) <- "surrogate_longitudinal_fit"
  out
}
## ---- Model specification helpers -------------------------------------------

#' No random effects
#' @param name Model name.
#'
#' @return A random-effects model specification.
#'
#' @export
make_longitudinal_re_none <- function(name = "Marginal model") {
  list(name = name, re_mode = "none")
}

#' Endpoint-specific random effects within treatment arm
#' @param name Model name.
#' @param S_0,T_0,S_1,T_1 Random-effects structure for the surrogate and true
#' endpoints in the control and treatment arms. One of \code{"none"},
#' \code{"ri"}, or \code{"rislope"}.
#'
#' @return A random-effects model specification.
#'
#' @export
make_longitudinal_re_endpoint_specific <- function(name,
                                                   S_0 = "none",
                                                   T_0 = "none",
                                                   S_1 = "none",
                                                   T_1 = "none") {
  allowed <- c("none", "ri", "rislope")
  vals <- c(S_0 = S_0, T_0 = T_0, S_1 = S_1, T_1 = T_1)
  if (any(!vals %in% allowed)) {
    stop("Random-effects types must be one of: ", paste(allowed, collapse = ", "))
  }
  list(name = name,
       re_mode = "endpoint_specific",
       re_struct = as.list(vals))
}

#' Random effects shared by the two endpoints within each treatment arm
#' @param name Model name.
#' @param arm0,arm1 Random-effects structure in the control and treatment arms.
#' One of \code{"none"}, \code{"ri"}, or \code{"rislope"}.
#'
#' @return A random-effects model specification.
#'
#' @export
make_longitudinal_re_shared_by_arm <- function(name,
                                               arm0 = "none",
                                               arm1 = "none") {
  allowed <- c("none", "ri", "rislope")
  vals <- c(`0` = arm0, `1` = arm1)
  if (any(!vals %in% allowed)) {
    stop("Random-effects types must be one of: ", paste(allowed, collapse = ", "))
  }
  list(name = name,
       re_mode = "shared_by_arm",
       shared_re = as.list(vals))
}


## ---- Internal covariance utilities -----------------------------------------

.surr_make_chol_from_par <- function(par, q) {
  L <- matrix(0, q, q)
  idx <- 1L
  for (r in seq_len(q)) {
    for (c in seq_len(r)) {
      L[r, c] <- if (r == c) exp(par[idx]) else par[idx]
      idx <- idx + 1L
    }
  }
  L
}

.surr_pack_chol <- function(S) {
  S <- 0.5 * (S + t(S))
  L <- t(chol(S))
  out <- numeric(nrow(S) * (nrow(S) + 1L) / 2L)
  idx <- 1L
  for (r in seq_len(nrow(L))) {
    for (c in seq_len(r)) {
      out[idx] <- if (r == c) log(L[r, c]) else L[r, c]
      idx <- idx + 1L
    }
  }
  out
}

.surr_safe_pack_chol <- function(S, jitter = 1e-8) {
  S <- 0.5 * (S + t(S))
  out <- tryCatch(.surr_pack_chol(S), error = function(e) NULL)
  if (!is.null(out)) return(out)

  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, jitter)
  S_pd <- ee$vectors %*% diag(vals, nrow = length(vals)) %*% t(ee$vectors)
  .surr_pack_chol(0.5 * (S_pd + t(S_pd)))
}

.surr_n_chol_par <- function(q) q * (q + 1L) / 2L

.surr_re_type_endpoint <- function(model_spec, endpoint, arm) {
  if (model_spec$re_mode == "none") return("none")
  if (model_spec$re_mode != "endpoint_specific") {
    stop("Endpoint-specific random-effect lookup requested for incompatible model specification.")
  }
  key <- paste0(endpoint, "_", arm)
  out <- model_spec$re_struct[[key]]
  if (is.null(out)) stop("Missing random-effects structure for ", key, ".")
  out
}

.surr_re_type_shared <- function(model_spec, arm) {
  if (model_spec$re_mode == "none") return("none")
  if (model_spec$re_mode != "shared_by_arm") {
    stop("Shared random-effect lookup requested for incompatible model specification.")
  }
  out <- model_spec$shared_re[[as.character(arm)]]
  if (is.null(out)) stop("Missing shared random-effects structure for arm ", arm, ".")
  out
}

.surr_random_effect_names <- function(model_spec, arm, endpoint_order = c("S", "T")) {
  if (model_spec$re_mode == "none") return(character(0))

  if (model_spec$re_mode == "endpoint_specific") {
    out <- character(0)
    for (ep in endpoint_order) {
      re_type <- .surr_re_type_endpoint(model_spec, ep, arm)
      if (re_type == "ri") out <- c(out, paste0("ri_", ep))
      if (re_type == "rislope") out <- c(out, paste0("ri_", ep), paste0("rs_", ep))
    }
    return(out)
  }

  if (model_spec$re_mode == "shared_by_arm") {
    re_type <- .surr_re_type_shared(model_spec, arm)
    if (re_type == "none") return(character(0))
    if (re_type == "ri") return("shared_ri")
    if (re_type == "rislope") return(c("shared_ri", "shared_rs"))
  }

  stop("Unknown `re_mode` in model specification.")
}

.surr_random_dim <- function(model_spec, arm, endpoint_order = c("S", "T")) {
  length(.surr_random_effect_names(model_spec, arm, endpoint_order))
}

.surr_make_Zrand <- function(di, model_spec, endpoint_order = c("S", "T")) {
  arm <- di$Z[1]
  n <- nrow(di)

  if (model_spec$re_mode == "none") return(NULL)

  if (model_spec$re_mode == "endpoint_specific") {
    nm <- .surr_random_effect_names(model_spec, arm, endpoint_order)
    if (length(nm) == 0L) return(NULL)
    Z <- matrix(0, n, length(nm), dimnames = list(NULL, nm))

    for (ep in endpoint_order) {
      re_type <- .surr_re_type_endpoint(model_spec, ep, arm)
      rows <- which(di$endpoint == ep)
      if (length(rows) == 0L) next
      if (re_type == "ri") {
        Z[rows, paste0("ri_", ep)] <- 1
      } else if (re_type == "rislope") {
        Z[rows, paste0("ri_", ep)] <- 1
        Z[rows, paste0("rs_", ep)] <- di$ctime[rows]
      }
    }
    return(Z)
  }

  if (model_spec$re_mode == "shared_by_arm") {
    re_type <- .surr_re_type_shared(model_spec, arm)
    if (re_type == "none") return(NULL)
    if (re_type == "ri") {
      Z <- matrix(1, n, 1, dimnames = list(NULL, "shared_ri"))
      return(Z)
    }
    if (re_type == "rislope") {
      Z <- cbind(shared_ri = 1, shared_rs = di$ctime)
      return(Z)
    }
  }

  stop("Unknown `re_mode` in model specification.")
}

.surr_rho_bounds <- function(residual_type, n_visits) {
  residual_type <- match.arg(residual_type, c("exponential", "independence", "cs"))
  if (residual_type == "independence") return(c(NA_real_, NA_real_))
  if (residual_type == "exponential") return(c(0, 1))
  if (n_visits <= 1L) return(c(-1, 1))
  c(-1 / (n_visits - 1), 1)
}

.surr_rho_from_eta <- function(eta, residual_type, n_visits) {
  bounds <- .surr_rho_bounds(residual_type, n_visits)
  if (residual_type == "independence") return(NA_real_)
  bounds[1] + (bounds[2] - bounds[1]) * plogis(eta)
}

.surr_eta_from_rho <- function(rho, residual_type, n_visits) {
  if (residual_type == "independence") return(numeric(0))
  bounds <- .surr_rho_bounds(residual_type, n_visits)
  if (!is.finite(rho) || rho <= bounds[1] || rho >= bounds[2]) {
    stop("Starting temporal correlation lies outside the admissible interval.")
  }
  qlogis((rho - bounds[1]) / (bounds[2] - bounds[1]))
}

## Construct a temporal residual correlation matrix
make_longitudinal_R <- function(times,
                                rho = NULL,
                                type = c("exponential", "independence", "cs")) {
  type <- match.arg(type)
  times <- as.numeric(times)
  if (length(times) < 1L || any(!is.finite(times))) stop("`times` must be finite numeric values.")

  if (type == "independence") return(1 * outer(times, times, "=="))

  if (length(rho) != 1L || !is.finite(rho)) stop("`rho` must be a single finite value.")

  if (type == "exponential") {
    if (rho <= 0 || rho >= 1) stop("For exponential correlation, `rho` must lie in (0,1).")
    return(rho ^ abs(outer(times, times, "-")))
  }

  m <- length(unique(times))
  lower <- if (m <= 1L) -1 else -1 / (m - 1)
  if (rho <= lower || rho >= 1) {
    stop("For compound symmetry, `rho` must lie in (", signif(lower, 5), ", 1).")
  }
  R <- matrix(rho, length(times), length(times))
  R[outer(times, times, "==")] <- 1
  R
}


## ---- Parameter unpacking ----------------------------------------------------

.surr_unpack_theta <- function(theta,
                               n_beta,
                               model_spec,
                               residual_type,
                               n_visits,
                               endpoint_order = c("S", "T")) {
  beta <- theta[seq_len(n_beta)]
  k <- n_beta

  rho <- NA_real_
  if (residual_type != "independence") {
    rho <- .surr_rho_from_eta(theta[k + 1L], residual_type, n_visits)
    k <- k + 1L
  }

  V0_par <- theta[(k + 1L):(k + 3L)]; k <- k + 3L
  V1_par <- theta[(k + 1L):(k + 3L)]; k <- k + 3L
  L0 <- .surr_make_chol_from_par(V0_par, 2L)
  L1 <- .surr_make_chol_from_par(V1_par, 2L)
  V0 <- L0 %*% t(L0)
  V1 <- L1 %*% t(L1)
  dimnames(V0) <- list(endpoint_order, endpoint_order)
  dimnames(V1) <- list(endpoint_order, endpoint_order)

  out <- list(beta = beta, rho = rho, V0 = V0, V1 = V1)

  for (arm in 0:1) {
    q <- .surr_random_dim(model_spec, arm, endpoint_order)
    if (q == 0L) next
    np <- .surr_n_chol_par(q)
    dpar <- theta[(k + 1L):(k + np)]; k <- k + np
    LD <- .surr_make_chol_from_par(dpar, q)
    D <- LD %*% t(LD)
    nm <- .surr_random_effect_names(model_spec, arm, endpoint_order)
    dimnames(D) <- list(nm, nm)
    out[[paste0("D", arm)]] <- D
  }

  if (k != length(theta)) stop("Internal parameter-length mismatch.")
  out
}




