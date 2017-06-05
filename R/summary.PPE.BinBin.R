summary.PPE.BinBin<-function(object, ..., Object){

if (missing(Object)){Object <- object} 
x <- Object  

mono=as.character(unique(x$Monotonicity))

mode <- function(data) {
    y <- data
    z <- density(y)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
 

for (j in 1:length(mono)){

valid<-length(x$index[x$Monotonicity==mono[j]])

cat("\n\n\n# Number of valid vectors",": Monotonicity=",mono[j],sep=" ")
cat("\n#---------------------------------------------------------------------------------------\n\n")  

cat("\n# N= ", valid,sep=" ")

cat("\n\n\n# Probability of a Prediction Error (PPE)",": Monotonicity=",mono[j],sep=" ")

 cat("\n#---------------------------------------------------------------------------------------\n\n")  
  cat("Mean (SD) PPE: ", format(round(mean(x$PPE[x$Monotonicity==mono[j]]), 4), nsmall = 4), " (", format(round(sd(x$PPE[x$Monotonicity==mono[j]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$PPE[x$Monotonicity==mono[j]]), 4), nsmall = 4), "; max: ",  format(round(max(x$PPE[x$Monotonicity==mono[j]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode PPE: ", format(round(mode(x$PPE[x$Monotonicity==mono[j]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(x$PPE[x$Monotonicity==mono[j]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# Reduction in Predicton Error attributed to the Surrogate (RPE)",": Monotonicity=",mono[j],sep=" ")

 cat("\n#---------------------------------------------------------------------------------------------------------\n\n")  
  cat("Mean (SD) RPE: ", format(round(mean(x$RPE[x$Monotonicity==mono[j]]), 4), nsmall = 4), " (", format(round(sd(x$RPE[x$Monotonicity==mono[j]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$RPE[x$Monotonicity==mono[j]]), 4), nsmall = 4), "; max: ",  format(round(max(x$RPE[x$Monotonicity==mono[j]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode Pe: ", format(round(mode(x$RPE[x$Monotonicity==mono[j]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(x$RPE[x$Monotonicity==mono[j]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# Probability of a Prediction Error based on T only (PPE_T)",": Monotonicity=",mono[j],sep=" ")

 cat("\n#---------------------------------------------------------------------------------------------------------\n\n")  
  cat("Mean (SD) PPE_T: ", format(round(mean(x$PPE_T[x$Monotonicity==mono[j]]), 4), nsmall = 4), " (", format(round(sd(x$PPE_T[x$Monotonicity==mono[j]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$PPE_T[x$Monotonicity==mono[j]]), 4), nsmall = 4), "; max: ",  format(round(max(x$PPE_T[x$Monotonicity==mono[j]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode PPE_T: ", format(round(mode(x$PPE_T[x$Monotonicity==mono[j]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(x$PPE_T[x$Monotonicity==mono[j]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# Individual Causal Association (R2_H)",": Monotonicity=",mono[j],sep=" ")

 cat("\n#---------------------------------------------------------------------------------------------------------\n\n")  
  cat("Mean (SD) ICA: ", format(round(mean(x$R2_H[x$Monotonicity==mono[j]]), 4), nsmall = 4), " (", format(round(sd(x$R2_H[x$Monotonicity==mono[j]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$R2_H[x$Monotonicity==mono[j]]), 4), nsmall = 4), "; max: ",  format(round(max(x$R2_H[x$Monotonicity==mono[j]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode ICA	: ", format(round(mode(x$R2_H[x$Monotonicity==mono[j]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(x$R2_H[x$Monotonicity==mono[j]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)

}
}

