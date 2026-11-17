#' Prepare data for longitudinal surrogate modeling
#'
#'
#' This function standardizes a long-format data set to the internal variables used by the
#' longitudinal model-fitting functions.
#'
#' @param data A data.frame in long format.
#' @param id Character name of the subject identifier column.
#' @param treatment Character name of the treatment column.
#' @param time Character name of the measurement-time column.
#' @param endpoint Character name of the endpoint-label column.
#' @param outcome Character name of the observed outcome column.
#' @param surrogate_label Label identifying the surrogate endpoint.
#' @param true_label Label identifying the true endpoint.
#' @param treatment_levels Two treatment values, ordered as control then active.
#' @param covariates Optional character vector of additional baseline covariates.
#'
#' @return A standardized data.frame.
#'
#'
#' @details Baseline responses must be supplied as ordinary outcome rows with
#' \code{time = 0}. Post-randomization observations must have strictly
#'  positive times. This single convention allows baseline outcomes to
#'  contribute to estimation of the augmented treatment-arm marginal models,
#'  while the longitudinal ICA is evaluated only at positive times.
#'
#' @examples
#' data(Schizo_PANSS_BPRS)
#' dat <- prepare_longitudinal_surrogate_data(
#'        data = Schizo_PANSS_BPRS,
#'        id = "Id",
#'        treatment = "Treat",
#'        time = "Time_cont",
#'        endpoint = "Endpoint",
#'        outcome = "Response_log",
#'        surrogate_label = "BPRS",
#'        true_label = "PANSS",
#'        treatment_levels = c(0, 1))
#'
#' @export
prepare_longitudinal_surrogate_data <- function(data,
                                                id,
                                                treatment,
                                                time,
                                                endpoint,
                                                outcome,
                                                surrogate_label,
                                                true_label,
                                                treatment_levels = c(0, 1),
                                                covariates = character(0)) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (length(treatment_levels) != 2L) stop("`treatment_levels` must have length 2.")

  required <- c(id, treatment, time, endpoint, outcome, covariates)
  required <- unique(required[!is.na(required) & nzchar(required)])
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  out <- data.frame(
    id = data[[id]],
    time = data[[time]],
    endpoint_raw = data[[endpoint]],
    Y = data[[outcome]],
    stringsAsFactors = FALSE
  )

  out$treatF <- factor(data[[treatment]],
                       levels = treatment_levels,
                       labels = c("0", "1"))
  if (anyNA(out$treatF)) {
    stop("Some treatment values could not be mapped to `treatment_levels`.")
  }
  out$Z <- as.integer(out$treatF) - 1L

  out$endpoint <- ifelse(out$endpoint_raw == surrogate_label, "S",
                         ifelse(out$endpoint_raw == true_label, "T", NA_character_))
  if (anyNA(out$endpoint)) {
    stop("Some endpoint labels were not recognized. Check `surrogate_label` and `true_label`.")
  }
  out$endpoint <- factor(out$endpoint, levels = c("S", "T"))
  out$epIdx <- ifelse(out$endpoint == "S", 1L, 2L)

  if (is.numeric(out$time)) {
    out$time_cont <- as.numeric(out$time)
  } else {
    suppressWarnings(time_num <- as.numeric(out$time))
    if (all(!is.na(time_num))) {
      out$time_cont <- time_num
    } else {
      stop("`time` must be numeric or coercible to numeric so that baseline is ",
           "represented explicitly by time = 0.")
    }
  }

  if (any(!is.finite(out$time_cont))) stop("Measurement times must be finite.")
  if (any(out$time_cont < 0)) {
    stop("Measurement times must be nonnegative: use time = 0 for baseline ",
         "and strictly positive times after randomization.")
  }
  out$visitF <- factor(out$time_cont, levels = sort(unique(out$time_cont)))
  out$is_baseline <- out$time_cont == 0
  out$is_postbaseline <- out$time_cont > 0

  if (length(covariates) > 0L) {
    for (cc in covariates) {
      out[[cc]] <- data[[cc]]
      if (is.character(out[[cc]])) out[[cc]] <- factor(out[[cc]])
    }
  }

  duplicated_rows <- duplicated(out[c("id", "time_cont", "endpoint")])
  if (any(duplicated_rows)) {
    stop("Each subject may have at most one observation per endpoint and time.")
  }
  if (!any(out$is_baseline)) {
    warning("No time-0 baseline responses were found. The model will be fitted ",
            "to post-randomization responses only.", call. = FALSE)
  }
  if (anyNA(out$id)) stop("Subject identifiers must not be missing.")
  treatment_counts <- tapply(out$Z, out$id, function(x) length(unique(x)))
  if (any(treatment_counts != 1L)) {
    stop("Treatment assignment must be constant within subject.")
  }

  ## Time is anchored at the explicitly defined baseline. This makes random
  ## intercepts interpretable at time 0 and avoids a hidden centering rule.
  out$ctime <- out$time_cont

  out <- out[order(out$id, out$time_cont, out$epIdx), , drop = FALSE]
  rownames(out) <- NULL
  out
}
