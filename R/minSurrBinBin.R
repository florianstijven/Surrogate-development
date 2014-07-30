MinSurrBinBin <- function(Object, Monotonicity=c("No"), delta=0.5){
    
  if (Monotonicity=="No") {
    Delta_T=c(3)}
  if (Monotonicity=="True.Endp") {
    Delta_T=c(2)}
  if (Monotonicity=="Surr.Endp") {
    Delta_T=c(3)}
  if (Monotonicity=="Surr.True.Endp") {
    Delta_T=c(2)}
  
  R2_HL_all <- as.numeric(1 - ((delta * log2(Delta_T) + 1)/ Object$H_Delta_T))
  
  fit <- 
    list(R2_HL = as.numeric(R2_HL_all), Call=match.call())   
  
  class(fit) <- "MinSurrBinBin"
  fit
}

summary.MinSurrBinBin <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n\n# R2_HL results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) R2_HL: ", format(round(mean(Object$R2_HL[which(is.finite(Object$R2_HL))]), 4), nsmall = 4), " (", format(round(sd(Object$R2_HL[which(is.finite(Object$R2_HL))]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$R2_HL[which(is.finite(Object$R2_HL))]), 4), nsmall = 4), "; max: ",  format(round(max(Object$R2_HL[which(is.finite(Object$R2_HL))]), 4), nsmall = 4), "]", sep="")
  cat("\n\nQuantiles of the R2_HL distribution: \n\n")
  quant <- quantile(Object$R2_HL, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
}
