summary.Fano.BinBin<-function(object,...,Object, Type="Overall"){

if (missing(Object)){Object<-object}
x <- Object
  
mode <- function(data) {
    y <- data
    z <- density(y)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }

Delta=as.numeric(unique(x$delta))


if (Type=="Overall"){

for (j in 1:length(Delta)){

bydelta<-x$R2_HL[x$delta==Delta[j]]
  

cat("R2_HL results summary for","Delta=",Delta[j],sep=" ")


  cat("\n#-------------------------------------------------------\n\n\n")  
  cat("Mean (SD) R2_HL: ", format(round(mean(x$R2_HL[x$delta==Delta[j]]), 4), nsmall = 4), " (", format(round(sd(x$R2_HL[x$delta==Delta[j]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$R2_HL[x$delta==Delta[j]]), 4), nsmall = 4), "; max: ",  format(round(max(x$R2_HL[x$delta==Delta[j]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(x$R2_HL)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the R2_HL distribution: \n\n")
  quant <- quantile(x$R2_HL[x$delta==Delta[j]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  cat("\n\n\n")  

}


cat("\n\n\n# H_Delta_T summary")

 cat("\n#--------------------------------------------------------\n\n")  
  cat("Mean (SD) H_Delta_T: ", format(round(mean(x$H_Delta_T[x$delta==Delta[1]]), 4), nsmall = 4), " (", format(round(sd(x$H_Delta_T[x$delta==Delta[1]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$H_Delta_T[x$delta==Delta[1]]), 4), nsmall = 4), "; max: ",  format(round(max(x$H_Delta_T[x$delta==Delta[1]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(x$H_Delta_T[x$delta==Delta[1]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the H_Delta_T distribution: \n\n")
  quant <- quantile(x$H_Delta_T[x$delta==Delta[1]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# pi_00 summary")

 cat("\n#--------------------------------------------------------\n\n")  
  cat("Mean (SD) pi_00: ", format(round(mean(x$pi_00[x$delta==Delta[1]]), 4), nsmall = 4), " (", format(round(sd(x$pi_00[x$delta==Delta[1]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$pi_00[x$delta==Delta[1]]), 4), nsmall = 4), "; max: ",  format(round(max(x$pi_00[x$delta==Delta[1]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(x$pi_00[x$delta==Delta[1]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the pi_00 distribution: \n\n")
  quant <- quantile(x$pi_00[x$delta==Delta[1]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# pi_11 summary")

 cat("\n#--------------------------------------------------------\n\n")  
  cat("Mean (SD) pi_11: ", format(round(mean(x$pi_11[x$delta==Delta[1]]), 4), nsmall = 4), " (", format(round(sd(x$pi_11[x$delta==Delta[1]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$pi_11[x$delta==Delta[1]]), 4), nsmall = 4), "; max: ",  format(round(max(x$pi_11[x$delta==Delta[1]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(x$pi_11[x$delta==Delta[1]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the pi_11 distribution: \n\n")
  quant <- quantile(x$pi_11[x$delta==Delta[1]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)

cat("\n\n\n# pi_01 summary")

 cat("\n#--------------------------------------------------------\n\n")  
  cat("Mean (SD) pi_01: ", format(round(mean(x$pi_01[x$delta==Delta[1]]), 4), nsmall = 4), " (", format(round(sd(x$pi_01[x$delta==Delta[1]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$pi_01[x$delta==Delta[1]]), 4), nsmall = 4), "; max: ",  format(round(max(x$pi_01[x$delta==Delta[1]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(x$pi_01[x$delta==Delta[1]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the pi_01 distribution: \n\n")
  quant <- quantile(x$pi_01[x$delta==Delta[1]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)


cat("\n\n\n# pi_10 summary")

 cat("\n#--------------------------------------------------------\n\n")  
  cat("Mean (SD) pi_10: ", format(round(mean(x$pi_10[x$delta==Delta[1]]), 4), nsmall = 4), " (", format(round(sd(x$pi_10[x$delta==Delta[1]]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(x$pi_10[x$delta==Delta[1]]), 4), nsmall = 4), "; max: ",  format(round(max(x$pi_10[x$delta==Delta[1]]), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(x$pi_10[x$delta==Delta[1]])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the pi_10 distribution: \n\n")
  quant <- quantile(x$pi_10[x$delta==Delta[1]], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
}


if (Type=="Mean"){


Delta=as.numeric(unique(x$delta))
pi10=as.numeric(unique(x$samplepi10))
meanvector <- matrix(, ncol = length(Delta), nrow = length(pi10))


  for (i in 1:length(pi10)){
   for (j in 1:length(Delta)){
     meanvector[i,j]=mean(x$R2_HL[x$delta==Delta[j] & x$samplepi10==pi10[i]])
      }
   }

colnames(meanvector)=Delta

for (j in 1:length(Delta)){


cat("R2_HL results summary for","Delta=",Delta[j],sep=" ")


  cat("\n#-------------------------------------------------------\n\n\n")  
  cat("Mean (SD) R2_HL: ", format(round(mean(meanvector[,j]), 4), nsmall = 4), " (", format(round(sd(meanvector[,j]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(meanvector[,j]), 4), nsmall = 4), "; max: ",  format(round(max(meanvector[,j]), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(meanvector[,j])$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the R2_HL distribution: \n\n")
  quant <- quantile(meanvector[,j], probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  cat("\n\n\n")  

}
}
}

