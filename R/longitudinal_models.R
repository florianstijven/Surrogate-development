## ============================================================================
## User-facing model specifications for longitudinal surrogate evaluation
## ============================================================================

#' Construct a standard longitudinal surrogate-model specification
#'
#' @param model Character model name. Supported values are \code{"galecki"},
#'   \code{"endpoint_ri"}, \code{"endpoint_ri_rs"}, \code{"shared_ri"}, and
#'   \code{"shared_ri_rs"}.
#' @return A random-effects specification accepted by
#'   \code{fit_longitudinal_surrogate_model()}.
#' @export
longitudinal_model_spec <- function(model) {
  model <- match.arg(
    model,
    c("galecki", "endpoint_ri", "endpoint_ri_rs",
      "shared_ri", "shared_ri_rs")
  )

  switch(
    model,
    galecki = make_longitudinal_re_none(
      name = "Galecki's marginal model"
    ),
    endpoint_ri = make_longitudinal_re_endpoint_specific(
      name = "Endpoint-specific random-intercept model",
      S_0 = "ri", T_0 = "ri", S_1 = "ri", T_1 = "ri"
    ),
    endpoint_ri_rs = make_longitudinal_re_endpoint_specific(
      name = "Endpoint-specific random-intercept-and-slope model",
      S_0 = "rislope", T_0 = "rislope",
      S_1 = "rislope", T_1 = "rislope"
    ),
    shared_ri = make_longitudinal_re_shared_by_arm(
      name = "Shared random-intercept model",
      arm0 = "ri", arm1 = "ri"
    ),
    shared_ri_rs = make_longitudinal_re_shared_by_arm(
      name = "Shared random-intercept-and-slope model",
      arm0 = "rislope", arm1 = "rislope"
    )
  )
}


#' Fit a standard bivariate longitudinal surrogate model
#'
#' This convenience function standardizes a long-format data set and fits one
#' of the standard covariance models. The surrogate and true endpoints may be
#' observed at different times.
#'
#' @param data A data.frame in long format.
#' @param model One of the model names accepted by
#'   \code{longitudinal_model_spec()}.
#' @param formula Fixed-effects formula written in terms of the standardized
#'   response \code{Y}.
#' @param id,treatment,time,endpoint,outcome Column names in \code{data}.
#' @param surrogate_label,true_label Values identifying the two endpoints.
#' @param treatment_levels Control and active treatment values, in that order.
#' @param covariates Optional additional covariate columns.
#' @param residual_type Temporal residual-correlation structure.
#' @param start_fit Optional fitted model used for initialization.
#' @param control Control list passed to \code{nlminb()}.
#' @details Baseline responses are supplied in the same long data set as the
#'   repeated outcomes, using \code{time = 0}; no separate baseline argument is
#'   accepted. The baseline responses contribute to estimation, but time 0 is
#'   never part of the longitudinal ICA evaluation grid.
#' @return An object of class \code{surrogate_longitudinal_fit}.
#' @export
fit_longitudinal_surrogate <- function(data,
                                       model,
                                       formula,
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
                                       start_fit = NULL,
                                       control = list(rel.tol = 1e-10,
                                                      x.tol = 1e-10,
                                                      eval.max = 5000,
                                                      iter.max = 3000)) {
  dat <- prepare_longitudinal_surrogate_data(
    data = data,
    id = id,
    treatment = treatment,
    time = time,
    endpoint = endpoint,
    outcome = outcome,
    surrogate_label = surrogate_label,
    true_label = true_label,
    treatment_levels = treatment_levels,
    covariates = covariates
  )

  spec <- longitudinal_model_spec(model)
  fit <- fit_longitudinal_surrogate_model(
    data = dat,
    formula = formula,
    model_spec = spec,
    residual_type = residual_type,
    start_fit = start_fit,
    control = control
  )
  fit$model <- model
  fit$data_roles <- list(
    id = id,
    treatment = treatment,
    time = time,
    endpoint = endpoint,
    outcome = outcome,
    surrogate_label = surrogate_label,
    true_label = true_label,
    treatment_levels = treatment_levels,
    covariates = covariates
  )
  fit
}
