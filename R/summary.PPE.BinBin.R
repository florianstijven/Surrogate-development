
summary.PPE.BinBin<-function(object, ..., Object){

if (missing(Object)){Object <- object} 
x <- Object
  
mode <- function(data) {
    y <- data
    z <- density(y)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }
 

valid<-length(x$RPE)

cat("\n\n\n# Number of valid vectors",sep=" ")
cat("\n#---------------------------------------------------------------------------------------\n\n")  

cat("\n# N= ", valid,sep=" ")

cat("\n\n\n# Probability of a Prediction Error (PPE)",sep=" ")

 cat("\n#---------------------------------------------------------------------------------------\n\n")  
  cat("Mean (SD) PPE: ", format(round(mean(x$PPE), 4), nsmall = 4), " (", format(round(sd(x$PPE), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$PPE), 4), nsmall = 4), "; max: ",  format(round(max(x$PPE), 4), nsmall = 4), "]", sep="")
  cat("\nMode PPE: ", format(round(mode(x$PPE)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(x$PPE, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# Reduction in Predicton Error attributed to the Surrogate (RPE)",sep=" ")

 cat("\n#---------------------------------------------------------------------------------------------------------\n\n")  
  cat("Mean (SD) RPE: ", format(round(mean(x$RPE), 4), nsmall = 4), " (", format(round(sd(x$RPE), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$RPE), 4), nsmall = 4), "; max: ",  format(round(max(x$RPE), 4), nsmall = 4), "]", sep="")
  cat("\nMode Pe: ", format(round(mode(x$RPE)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(x$RPE, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# Probability of a Prediction Error based on T only (PPE_T)",sep=" ")

 cat("\n#---------------------------------------------------------------------------------------------------------\n\n")  
  cat("Mean (SD) PPE_T: ", format(round(mean(x$PPE_T), 4), nsmall = 4), " (", format(round(sd(x$PPE_T), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$PPE_T), 4), nsmall = 4), "; max: ",  format(round(max(x$PPE_T), 4), nsmall = 4), "]", sep="")
  cat("\nMode PPE_T: ", format(round(mode(x$PPE_T)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(x$PPE_T, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# Individual Causal Association (R2_H)",sep=" ")

 cat("\n#---------------------------------------------------------------------------------------------------------\n\n")  
  cat("Mean (SD) ICA: ", format(round(mean(x$R2_H), 4), nsmall = 4), " (", format(round(sd(x$R2_H), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$R2_H), 4), nsmall = 4), "; max: ",  format(round(max(x$R2_H), 4), nsmall = 4), "]", sep="")
  cat("\nMode ICA	: ", format(round(mode(x$R2_H)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(x$R2_H, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)

}

