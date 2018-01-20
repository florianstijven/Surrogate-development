summary.SPF.BinCont <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  
  mode <- function(data) {
    x <- data
    if (unique(x[1])!=0){
      z <- density(x)
      mode_val <- z$x[which.max(z$y)]
      if (mode_val < 0){mode_val <- c(0)}
    }
    if (unique(x[1])==0){
      mode_val <- c(0)
    }  
    fit <- list(mode_val= mode_val)  
  }
  
  suppressWarnings(warning("mode"))
  
  options(digits=5)
  
  cat("Requested interval Delta S: [",Object$a, ", ", Object$b, "]", sep="")
  cat("\n\nSPF Descriptives")
  cat("\n~~~~~~~~~~~~~~~~")
  
  DT <- paste("Delta S in interval [",Object$a, ", ", Object$b, "]", sep="")
  
  cat("\n\nP[Delta T = -1 | ", DT, ":", sep="")
  cat("\n----------------------------------------------")
  cat("\nMean: ", mean(na.exclude(Object$P_Delta_T_min1)), ";  Median: ", median((na.exclude(Object$P_Delta_T_min1))), 
          ";  Mode: ", mode((na.exclude(Object$P_Delta_T_min1)))$mode_val,
          ";  SD: ", sd((na.exclude(Object$P_Delta_T_min1))), "\nMin: ", min(na.exclude(Object$P_Delta_T_min1)), "; Max: ", 
      max(na.exclude(Object$P_Delta_T_min1)), 
          "; 95% CI = [", quantile((na.exclude(Object$P_Delta_T_min1)), probs = c(.025)), "; ",  quantile((na.exclude(Object$P_Delta_T_min1)), probs = c(.975)), "]\n",
          sep="")
  
  cat("\n\nP[Delta T = 0 | ", DT, ":", sep="")
  cat("\n----------------------------------------------")
  cat("\nMean: ", mean(na.exclude(Object$P_Delta_T_0)), ";  Median: ", median((na.exclude(Object$P_Delta_T_0))), 
      ";  Mode: ", mode((na.exclude(Object$P_Delta_T_0)))$mode_val,
      ";  SD: ", sd((na.exclude(Object$P_Delta_T_0))), "\nMin: ", min(na.exclude(Object$P_Delta_T_0)), "; Max: ", 
      max(na.exclude(Object$P_Delta_T_0)), 
      "; 95% CI = [", quantile((na.exclude(Object$P_Delta_T_0)), probs = c(.025)), "; ",  quantile((na.exclude(Object$P_Delta_T_0)), probs = c(.975)), "]\n",
      sep="")
  
  cat("\n\nP[Delta T = 1 | ", DT, ":", sep="")
  cat("\n----------------------------------------------")
  cat("\nMean: ", mean(na.exclude(Object$P_Delta_T_1)), ";  Median: ", median((na.exclude(Object$P_Delta_T_1))), 
      ";  Mode: ", mode((na.exclude(Object$P_Delta_T_1)))$mode_val,
      ";  SD: ", sd((na.exclude(Object$P_Delta_T_1))), "\nMin: ", min(na.exclude(Object$P_Delta_T_1)), "; Max: ", 
      max(na.exclude(Object$P_Delta_T_1)), 
      "; 95% CI = [", quantile((na.exclude(Object$P_Delta_T_1)), probs = c(.025)), "; ",  quantile((na.exclude(Object$P_Delta_T_1)), probs = c(.975)), "]\n",
      sep="")
  
  }
