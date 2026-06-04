#' @export
summary.Single.Trial.ContCont <- function(object, ..., Object){

  if (missing(Object)){Object <- object}
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1])
  cat("\nTotal number of patients in experimental treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat==1]),
      "\nTotal number of patients in control treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat!=1]))

  means_table <- rbind(tapply(Object$Data.Analyze$Surr, list(Object$Data.Analyze$Treat), mean), tapply(Object$Data.Analyze$True, list(Object$Data.Analyze$Treat), mean))
  colnames(means_table) <- c("Control Treatment", "Experimental treatment")
  rownames(means_table) <- c("Surrogate", "True endpoint")
  cat("\n\nMean surrogate and true endpoint values in each treatment group: \n\n")
  print(format(round(data.frame(means_table, stringsAsFactors = TRUE), 4), nsmall = 4))
  Var_table <- rbind(tapply(Object$Data.Analyze$Surr, list(Object$Data.Analyze$Treat), var), tapply(Object$Data.Analyze$True, list(Object$Data.Analyze$Treat), var))
  colnames(Var_table) <- c("Control Treatment", "Experimental treatment")
  rownames(Var_table) <- c("Surrogate", "True endpoint")
  cat("\n\nVar surrogate and true endpoint values in each treatment group: \n\n")
  print(format(round(data.frame(Var_table, stringsAsFactors = TRUE), 4), nsmall = 4))


  cat("\n\nCorrelations between the true and surrogate endpoints in the control (r_T0S0)")
  cat("\nand the experimental treatment groups (r_T1S1):\n\n")
  print(round(Object$Cor.Endpoints, 4), nsmall = 4)

  cat("\n\n\n# Prentice criteria")
  cat("\n#~~~~~~~~~~~~~~~~~~\n")
  
  cat("# Estimated parameters models 1-4")
  cat("\nalpha = ", round(x = Object$Prentice.Model.1[2,1], digits=4), ", se = ", round(x = Object$Prentice.Model.1[2,2], digits=4),  
      ", p-value = ", Object$Prentice.Model.1[2,4], sep="")
  cat("\nbeta = ", round(x = Object$Prentice.Model.2[2,1], digits=4), ", se = ", round(x = Object$Prentice.Model.2[2,2], digits=4),
      ", p-value = ", Object$Prentice.Model.2[2,4], sep="")
  cat("\ngamma = ", round(x = Object$Prentice.Model.3[2,1], digits=4), ", se = ", round(x = Object$Prentice.Model.3[2,2], digits=4),
      ", p-value = ", Object$Prentice.Model.3[2,4], sep="")
  cat("\nbeta_S = ", round(x = Object$Prentice.Model.4[2,1], digits=4), ", se = ", round(x = Object$Prentice.Model.4[2,2], digits=4),
      ", p-value = ", Object$Prentice.Model.4[2,4], sep="")
  
  cat("\n\nPrentice criteria passed: ")
  if (Object$Prentice.Passed == FALSE){cat("No")}
  if (Object$Prentice.Passed == TRUE){cat("Yes")}
  
  
  cat("\n\n\n# Proportion Explained")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~\n")
  cat("\nDelta method-based confidence interval:\n\n")
  print(format(round(Object$PE.Delta, 4), nsmall = 4))
  
  cat("\nFieller method-based confidence interval:\n\n")
  print(format(round(Object$PE.Fieller, 4), nsmall = 4))
  
  cat("\nBootstrap-based confidence interval:\n\n")
  print(format(round(Object$PE.Boot, 4), nsmall = 4))
  
  
  cat("\n\n# Relative Effect")
  cat("\n#~~~~~~~~~~~~~~~~")
  cat("\n\nDelta method-based confidence interval:\n\n")
  print(format(round(Object$RE.Delta, 4), nsmall = 4))
  cat("\nFieller theorem-based confidence interval:\n\n")
  print(format(round(Object$RE.Fieller, 4), nsmall = 4))
  cat("\nBootstrap-based confidence interval:\n\n")
  print(format(round(Object$RE.Boot, 4), nsmall = 4))
  
  
  cat("\n\n# Adjusted Association")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nFisher Z-based confidence interval:\n\n")
  print(format(round(Object$AA, 4), nsmall = 4))
  cat("\nBootstrap-based confidence interval:\n\n")
  print(format(round(Object$AA.Boot, 4), nsmall = 4))
}
