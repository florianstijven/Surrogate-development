AA.MultS <- function(Sigma_gamma, N, Alpha=0.05){
  dims <- dim(Sigma_gamma)[1] 
  K <- dims-1
  Alpha_orig <- Alpha
  sigma_TT <- matrix(data = Sigma_gamma[1,1], nrow = 1)
  sigma_ST <- matrix(data = Sigma_gamma[2:dims,1], nrow = (dims-1))
  sigma_SS <- matrix(data = Sigma_gamma[2:dims,2:dims], nrow = (dims-1))
  
  Gamma.Delta <- as.numeric((t(sigma_ST) %*% solve(sigma_SS) %*% sigma_ST) / sigma_TT)
  
  #Z <- .5*log((1+Gamma.Delta)/(1-Gamma.Delta)) 
  #rho_lb <- max((exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N-3)))))+1))
  #rho_ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N-3)))))+1))
  #rho_sd <- sqrt((1-Gamma.Delta**2)/(N-2))
  #rho_results_FishZ <- data.frame(cbind(Gamma.Delta, rho_sd , rho_lb, rho_ub))
  #colnames(rho_results_FishZ) <- c("Multiple-surrogate AA", "Standard Error", "CI lower limit", "CI upper limit")
  #rownames(rho_results_FishZ) <- c(" ")
  
  sd_val <- sqrt((4*Gamma.Delta*(1-Gamma.Delta)^2)/(N-3))
  lb <- max(0, as.numeric(as.numeric(Gamma.Delta) + qnorm(Alpha/2) * (sd_val)))
  ub <- min(1, as.numeric(Gamma.Delta + qnorm(1-Alpha/2)*(sd_val)))
  Gamma.Delta_Results <- data.frame(cbind(Gamma.Delta, sd_val, lb, ub))
  colnames(Gamma.Delta_Results) <- c("Multivariate AA", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Gamma.Delta_Results) <- c(" ") 
  
  N <- data.frame(N); rownames(N) <- c(" ")
  Alpha <- data.frame(Alpha); rownames(Alpha) <- c(" ")
  
  # Adjusted
  Alpha <- Alpha_orig
  Adj.Gamma.Delta2 <- 
    1 - (1 - Gamma.Delta) * ((N - 1) / (N - K - 1))  
  
  sd_val <- sqrt((4*Adj.Gamma.Delta2*(1-Adj.Gamma.Delta2)^2)/(N-3))
  lb <- max(0, as.numeric(as.numeric(Adj.Gamma.Delta2) + qnorm(Alpha/2) * (sd_val)))
  ub <- min(1, as.numeric(Adj.Gamma.Delta2 + qnorm(1-Alpha/2)*(sd_val)))
  Adj.Gamma.Delta2_Results <- data.frame(cbind(Adj.Gamma.Delta2, sd_val, lb, ub))
  colnames(Adj.Gamma.Delta2_Results) <- c("Adjusted multivariate AA", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(Adj.Gamma.Delta2_Results) <- c(" ") 
  
  #Adj.Gamma.Delta <- sqrt(Adj.Gamma.Delta2)
  
  #Z <- as.numeric(.5*log((1+Adj.Gamma.Delta)/(1-Adj.Gamma.Delta))) 
  #rho_lb <- max((exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N-3)))))+1))
  #rho_ub <- min(1, as.numeric(as.numeric(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N-3)))))+1)))
  #rho_sd <- sqrt((1-Gamma.Delta**2)/(N-2))
  #Adj.Gamma.Delta <- data.frame(cbind(Adj.Gamma.Delta, rho_sd , rho_lb, rho_ub))
  #colnames(Adj.Gamma.Delta) <- c("Corr. multiple-surrogate AA", "Standard Error", "CI lower limit", "CI upper limit")
  #rownames(Adj.Gamma.Delta) <- c(" ")
  
  
  
  
  fit <-   
    list(Gamma.Delta = Gamma.Delta_Results, Corr.Gamma.Delta=Adj.Gamma.Delta2_Results, Sigma_gamma=Sigma_gamma, N=N, Alpha=Alpha, Call=match.call())  
  
  class(fit) <- "AA.MultS"
  fit
}


summary.AA.MultS <- function(object, ..., Object) {
  if (missing(Object)) {
    Object <- object}
  cat("\nFunction call:\n\n")
  print(Object$Call)
  
  cat("\n\n# Uncorrected multivariate Adjusted Association")
  cat("\n# Fisher Z ", (1-as.numeric(Object$Alpha))*100, "%-based confidence interval (N = ", as.numeric(Object$N), ")", sep="")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  print(format(round(Object$Gamma.Delta, 4), nsmall = 4))

  cat("\n\n# Bias-corrected multivariate Adjusted Association")
  cat("\n# Fisher Z ", (1-as.numeric(Object$Alpha))*100, "%-based confidence interval (N = ", as.numeric(Object$N), ")", sep="")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  print(format(round(Object$Corr.Gamma.Delta, 4), nsmall = 4))
    }
