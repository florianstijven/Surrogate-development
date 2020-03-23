
summary.UnifixedContCont <- summary.UnimixedContCont <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1])
  cat("\nM(SD) patients per trial: ", format(round(mean((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), " (", format(round(sd((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), ")", 
      "  [min: ", min((Object$Obs.Per.Trial$Obs.per.trial)), "; max: ",  max((Object$Obs.Per.Trial$Obs.per.trial)), "]", sep="")
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
  cat("\n\n\n# Meta-analytic results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  print(format(round(Object$Trial.R2, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Indiv.R2, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Trial.R, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Indiv.R, 4), nsmall = 4))
  
  #ICA
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  cat("\n\n\n# ICA results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat("Mean (SD) ICA: ", format(round(mean(Object$ICA$ICA), 4), nsmall = 4), " (", format(round(sd(Object$ICA$ICA), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$ICA$ICA), 4), nsmall = 4), "; max: ",  format(round(max(Object$ICA$ICA), 4), nsmall = 4), "]", sep="")
  cat("\nMode ICA: ", format(round(mode(Object$ICA$ICA)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the ICA distribution: \n\n")
  quant <- quantile(Object$ICA$ICA, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  # Max ent results
  Fit <- MaxEntContCont(x=Object$ICA, T0T0 = Object$T0T0, T1T1 = Object$T1T1, S0S0 = Object$S0S0, S1S1 = Object$S1S1)
  cat("\n\nMaximum entropy ICA: ")
  cat(Fit$ICA.Max.Ent) #, " (with entropy = ", Object$Max.Ent, ")\n", sep = "")
  cat("\n\nObtained under correlation structure:\n")
  print(Fit$Table.ICA.Entropy[1,3:8], row.names = "")
  cat("\n")
}

summary.BifixedContCont <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1]) 
  cat("\nM(SD) patients per trial: ", format(round(mean((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), " (", format(round(sd((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), ")", 
      "  [min: ", min((Object$Obs.Per.Trial$Obs.per.trial)), "; max: ",  max((Object$Obs.Per.Trial$Obs.per.trial)), "]", sep="")
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
  cat("\n\n\n# Meta-analytic results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  print(format(round(Object$Trial.R2, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Indiv.R2, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Trial.R, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Indiv.R, 4), nsmall = 4))  
  
  #ICA
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  cat("\n\n\n# ICA results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat("Mean (SD) ICA: ", format(round(mean(Object$ICA$ICA), 4), nsmall = 4), " (", format(round(sd(Object$ICA$ICA), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$ICA$ICA), 4), nsmall = 4), "; max: ",  format(round(max(Object$ICA$ICA), 4), nsmall = 4), "]", sep="")
  cat("\nMode ICA: ", format(round(mode(Object$ICA$ICA)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the ICA distribution: \n\n")
  quant <- quantile(Object$ICA$ICA, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  # Max ent results
  Fit <- MaxEntContCont(x=Object$ICA, T0T0 = Object$T0T0, T1T1 = Object$T1T1, S0S0 = Object$S0S0, S1S1 = Object$S1S1)
  cat("\n\nMaximum entropy ICA: ")
  cat(Fit$ICA.Max.Ent) #, " (with entropy = ", Object$Max.Ent, ")\n", sep = "")
  cat("\n\nObtained under correlation structure:\n")
  print(Fit$Table.ICA.Entropy[1,3:8], row.names = "")
  cat("\n")
}


summary.BimixedContCont <- function(object, ..., Object){

  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1])
  cat("\nM(SD) patients per trial: ", format(round(mean((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), " (", format(round(sd((Object$Obs.Per.Trial$Obs.per.trial)), 4), nsmall = 4), ")", 
      "  [min: ", min((Object$Obs.Per.Trial$Obs.per.trial)), "; max: ",  max((Object$Obs.Per.Trial$Obs.per.trial)), "]", sep="")
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
  cat("\n\n\n# Meta-analytic results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  print(format(round(Object$Trial.R2, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Indiv.R2, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Trial.R, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Indiv.R, 4), nsmall = 4))  
  
  #ICA
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
  cat("\n\n\n# ICA results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat("Mean (SD) ICA: ", format(round(mean(Object$ICA$ICA), 4), nsmall = 4), " (", format(round(sd(Object$ICA$ICA), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$ICA$ICA), 4), nsmall = 4), "; max: ",  format(round(max(Object$ICA$ICA), 4), nsmall = 4), "]", sep="")
  cat("\nMode ICA: ", format(round(mode(Object$ICA$ICA)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the ICA distribution: \n\n")
  quant <- quantile(Object$ICA$ICA, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  # Max ent results
  Fit <- MaxEntContCont(x=Object$ICA, T0T0 = Object$T0T0, T1T1 = Object$T1T1, S0S0 = Object$S0S0, S1S1 = Object$S1S1)
  cat("\n\nMaximum entropy ICA: ")
  cat(Fit$ICA.Max.Ent) #, " (with entropy = ", Object$Max.Ent, ")\n", sep = "")
  cat("\n\nObtained under correlation structure:\n")
  print(Fit$Table.ICA.Entropy[1,3:8], row.names = "")
  cat("\n")
}
