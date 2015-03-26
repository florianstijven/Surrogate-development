summ.ICA.By.Mon <- function(x){
  
  Object <- x 
  
  Object$R2_H <- na.exclude(Object$R2_H)
  Object$Theta_T <- na.exclude(Object$Theta_T)
  Object$Theta_S <- na.exclude(Object$Theta_S)
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n# Total number of valid Pi vectors")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(dim(Object$Pi.Vectors)[1])
  
    results <- cbind.data.frame(Object$Pi.Vectors$Monotonicity, Object$R2_H)
    
    cat("\n\n# Summary of results obtained in different monotonicity scenarios")
    cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    
    cat("# R2_H results summary")
    cat("\n~~~~~~~~~~~~~~~~~~~~~\n\nMean:\n")
    print(tapply(results[,2], list(results[,1]), mean))
    
    pc50 <- function(x=x){
      quantile(x = x, probs = .5, na.rm = T)}
    cat("\nMedian:\n")
    print(tapply(results[,2], list(results[,1]), pc50))
    
    mode <- function(data) {
      x <- data
      z <- density(x)
      mode_val <- z$x[which.max(z$y)]
      fit <- list(mode_val= mode_val)
    }
    cat("\nMode:\n")
    print(data.frame(t(tapply(results[,2], list(results[,1]), mode)), row.names=""))
    
    cat("\nSD:\n")
    print(tapply(results[,2], list(results[,1]), sd))
    
    cat("\nMin:\n")
    print(tapply(results[,2], list(results[,1]), min))
    cat("\nMax:\n")
    print(tapply(results[,2], list(results[,1]), max))
    
#  results <- cbind.data.frame(Object$Pi.Vectors$Monotonicity, Object$C3)
  
#  cat("\n\n# C3 results summary")
#  cat("\n#~~~~~~~~~~~~~~~~~~~\n\nMean:\n")
#  print(tapply(results[,2], list(results[,1]), mean))
  
#  cat("\nMedian:\n")
#  print(tapply(results[,2], list(results[,1]), pc50))
  
#  cat("\nMode:\n")
#  print(data.frame(t(tapply(results[,2], list(results[,1]), mode)), row.names=""))
  
#  cat("\nSD:\n")
#  print(tapply(results[,2], list(results[,1]), sd))
  
#  cat("\nMin:\n")
#  print(tapply(results[,2], list(results[,1]), min))
#  cat("\nMax:\n")
#  print(tapply(results[,2], list(results[,1]), max))
  
  
  results <- cbind.data.frame(Object$Pi.Vectors$Monotonicity, Object$Theta_T)
  
  cat("\n\n# Theta_T results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~\n\nMean:\n")
  print(tapply(results[,2], list(results[,1]), mean))
  
  cat("\nMedian:\n")
  print(tapply(results[,2], list(results[,1]), pc50))
  
  cat("\nSD:\n")
  print(tapply(results[,2], list(results[,1]), sd))
  
  cat("\nMin:\n")
  print(tapply(results[,2], list(results[,1]), min))
  cat("\nMax:\n")
  print(tapply(results[,2], list(results[,1]), max))
  
  
  results <- cbind.data.frame(Object$Pi.Vectors$Monotonicity, Object$Theta_S)
  
  cat("\n\n# Theta_S results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~\n\nMean:\n")
  print(tapply(results[,2], list(results[,1]), mean))
  
  cat("\nMedian:\n")
  print(tapply(results[,2], list(results[,1]), pc50))
  
  cat("\nSD:\n")
  print(tapply(results[,2], list(results[,1]), sd))
  
  cat("\nMin:\n")
  print(tapply(results[,2], list(results[,1]), min))
  cat("\nMax:\n")
  print(tapply(results[,2], list(results[,1]), max))
  
}

