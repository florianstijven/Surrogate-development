summary.ICA.BinBin <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  
  Object$R2_H <- na.exclude(Object$R2_H)
  Object$Theta_T <- na.exclude(Object$Theta_T)
  Object$Theta_S <- na.exclude(Object$Theta_S)
  
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n# Total number of valid Pi vectors")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(dim(Object$Pi.Vectors)[1])
  cat("\n\n\n# R2_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) R2_H: ", format(round(mean(Object$R2_H), 4), nsmall = 4), " (", format(round(sd(Object$R2_H), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$R2_H), 4), nsmall = 4), "; max: ",  format(round(max(Object$R2_H), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(Object$R2_H)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the R2_H distribution: \n\n")
  quant <- quantile(Object$R2_H, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  cat("\n\n\n# R_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) R_H: ", format(round(mean(sqrt(Object$R2_H)), 4), nsmall = 4), " (", format(round(sd(sqrt(Object$R2_H)), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(sqrt(Object$R2_H)), 4), nsmall = 4), "; max: ",  format(round(max(sqrt(Object$R2_H)), 4), nsmall = 4), "]", sep="")
  cat("\nMode R_H: ", format(round(mode(sqrt(Object$R2_H))$mode_val, 4), nsmall = 4))
  
  cat("\n\nQuantiles of the R_H distribution: \n\n")
  quant <- quantile(sqrt(Object$R2_H), probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  
  cat("\n\n# Theta_T results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) Theta_T: ", format(round(mean(Object$Theta_T), 4), nsmall = 4), " (", format(round(sd(Object$Theta_T), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$Theta_T), 4), nsmall = 4), "; max: ",  format(round(max(Object$Theta_T), 4), nsmall = 4), "]", sep="")
  if (Object$Theta_T[1]!=Inf) {
    cat("\nMode Theta_T: ", format(round(mode(Object$Theta_T)$mode_val, 4), nsmall = 4))}
  cat("\n\nQuantiles of the Theta_T distribution: \n\n")
  quant <- quantile(Object$Theta_T, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  cat("\n\n# Theta_S results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) Theta_S: ", format(round(mean(Object$Theta_S), 4), nsmall = 4), " (", format(round(sd(Object$Theta_S), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$Theta_S), 4), nsmall = 4), "; max: ",  format(round(max(Object$Theta_S), 4), nsmall = 4), "]", sep="")
  if (Object$Theta_S[1]!=Inf) {
    cat("\nMode Theta_S: ", format(round(mode(Object$Theta_S)$mode_val, 4), nsmall = 4))}
  cat("\n\nQuantiles of the Theta_S distribution: \n\n")
  quant <- quantile(Object$Theta_S, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
}
