#' @export
summary.ICA.BinCont <- function(object, ..., Object){

  options(digits = 4)

  if (missing(Object)){Object <- object}

  Object$R2_H <- na.exclude(Object$R2_H)

  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }

  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n# Total number of valid R2_H values")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(length(Object$R2_H))

  cat("\n\n\n# R2_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat("Mean (SD) R2_H:", format(round(mean(Object$R2_H), 4), nsmall = 4), " (", format(round(sd(Object$R2_H), 4), nsmall = 4), ")")
  cat("\nMode R2_H:", format(round(mode(Object$R2_H)$mode_val, 4), nsmall = 4))
  cat("\nMedian R2_H:", round(median(Object$R2_H), 4))
  cat("\nRange R2_H: [", round(range(Object$R2_H)[1], 4), ",", round(range(Object$R2_H)[2], 4), "]")
  cat("\nSDI0.95 R2_H: [", round(quantile(Object$R2_H, probs=0.025), 4), ",", round(quantile(Object$R2_H, probs=0.975), 4), "]")
  cat("\nSDI0.90 R2_H: [", round(quantile(Object$R2_H, probs=0.05), 4), ",", round(quantile(Object$R2_H, probs=0.95), 4), "]")
  cat("\nSDI0.85 R2_H: [", round(quantile(Object$R2_H, probs=0.075), 4), ",", round(quantile(Object$R2_H, probs=0.925), 4), "]")
  cat("\nSDI0.80 R2_H: [", round(quantile(Object$R2_H, probs=0.1), 4), ",", round(quantile(Object$R2_H, probs=0.9), 4), "]")
  cat("\nSDI0.75 R2_H: [", round(quantile(Object$R2_H, probs=0.125), 4), ",", round(quantile(Object$R2_H, probs=0.875), 4), "]")
}


