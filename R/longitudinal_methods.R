## ============================================================================
## Print, summary, and plot methods
## ============================================================================

#' @export
print.surrogate_longitudinal_fit <- function(x, ...) {
  stats <- longitudinal_fit_stats(x)
  cat("Bivariate longitudinal surrogate model\n")
  cat("  Model:       ", stats$model, "\n", sep = "")
  cat("  Subjects:    ", x$n_subjects, "\n", sep = "")
  cat("  Observations:", x$n_observations, "\n")
  cat("  LogLik:      ", format(stats$logLik, digits = 7), "\n", sep = "")
  cat("  AIC:         ", format(stats$AIC, digits = 7), "\n", sep = "")
  cat("  Convergence: ", stats$convergence, "\n", sep = "")
  invisible(x)
}


#' @export
print.surrogate_longitudinal_ica <- function(x, ...) {
  cat("Longitudinal individual causal association\n")
  cat("  Model: ", x$model, "\n", sep = "")
  cat("  T evaluation times: ",
      paste(x$evaluation_times$T, collapse = ", "), "\n", sep = "")
  cat("  S evaluation times: ",
      paste(x$evaluation_times$S, collapse = ", "), "\n", sep = "")
  cat("  Admissible ICA values: ", length(x$R2_Lambda), "\n", sep = "")
  if (length(x$R2_Lambda) > 0L) {
    ss <- summary_longitudinal_ica(x)
    cat("  Range: [", format(ss$minimum, digits = 6), ", ",
        format(ss$maximum, digits = 6), "]\n", sep = "")
    cat("  Median: ", format(ss$median, digits = 6), "\n", sep = "")
  }
  invisible(x)
}


#' @export
summary.surrogate_longitudinal_ica <- function(object,
                                               threshold = 0.95,
                                               ...) {
  out <- summary_longitudinal_ica(object, threshold = threshold)
  attr(out, "model") <- object$model
  attr(out, "evaluation_times") <- object$evaluation_times
  class(out) <- c("summary.surrogate_longitudinal_ica", class(out))
  out
}


#' @export
print.summary.surrogate_longitudinal_ica <- function(x, ...) {
  cat("Longitudinal ICA sensitivity summary\n")
  cat("  Model: ", attr(x, "model"), "\n", sep = "")
  NextMethod("print")
  invisible(x)
}


#' @export
plot.surrogate_longitudinal_ica <- function(x,
                                            xlim = c(0, 1),
                                            main = "ECDF of longitudinal ICA",
                                            xlab = expression(R[Lambda]^2),
                                            ylab = "ECDF",
                                            ...) {
  if (length(x$R2_Lambda) < 1L) stop("The ICA result contains no values.")
  plot(
    stats::ecdf(x$R2_Lambda),
    do.points = FALSE,
    verticals = TRUE,
    xlim = xlim,
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  invisible(x)
}
