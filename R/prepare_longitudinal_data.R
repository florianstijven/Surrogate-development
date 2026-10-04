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
#' @param baseline Optional character name of an endpoint-specific baseline value.
#' @param covariates Optional character vector of additional baseline covariates.
#'
#' @return A standardized data.frame.
#'
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
#'        baseline = "Baseline_log",
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
                                                baseline = NULL,
                                                covariates = character(0)) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (length(treatment_levels) != 2L) stop("`treatment_levels` must have length 2.")

  required <- c(id, treatment, time, endpoint, outcome, baseline, covariates)
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

  time_raw <- data[[time]]

  if (is.numeric(time_raw)) {
    out$time_cont <- as.numeric(time_raw)
  } else {
    time_chr <- if (is.factor(time_raw)) as.character(time_raw) else time_raw
    suppressWarnings(time_num <- as.numeric(time_chr))

    if (all(!is.na(time_num))) {
      out$time_cont <- time_num
    } else {
      warning("`time` is not numeric; ordered integer scores are being used. ",
              "For continuous-time correlation, numeric visit times are preferable.")
      time_levels <- if (is.factor(time_raw)) levels(time_raw) else sort(unique(as.character(time_raw)))
      out$time_cont <- as.numeric(factor(as.character(time_raw), levels = time_levels))
    }
  }

  if (any(!is.finite(out$time_cont))) stop("Measurement times must be finite.")
  out$visitF <- factor(out$time_cont, levels = sort(unique(out$time_cont)))

  if (!is.null(baseline)) out$Ybase <- data[[baseline]]

  if (length(covariates) > 0L) {
    for (cc in covariates) {
      out[[cc]] <- data[[cc]]
      if (is.character(out[[cc]])) out[[cc]] <- factor(out[[cc]])
    }
  }

  visit_times <- sort(unique(out$time_cont))
  out$ctime <- out$time_cont - mean(visit_times)

  out <- out[order(out$id, out$time_cont, out$epIdx), , drop = FALSE]
  rownames(out) <- NULL
  out
}
