# Bootstrap around SPF values
# library(Surrogate); data("Schizo_Bin")

# input
# Data = Schizo_Bin; Surr = "BPRS_Bin"; True = "PANSS_Bin"; Treat = "Treat";  M = 5; Seed=123  

# Voorbeeld
#MEP_CI <- Bootstrap.MEP.BinBin(Data = Schizo_Bin, Surr = "BPRS_Bin", True = "PANSS_Bin",
#                     Treat = "Treat", M = 3, Seed=123)
#summary(MEP_CI)



# start function
Bootstrap.MEP.BinBin <- function(Data, Surr, True, Treat, M=100, Seed=123){

All_r_1_1 <- All_r_min1_1 <- All_r_0_1 <- All_r_1_0 <- All_r_min1_0 <- All_r_0_0 <- All_r_1_min1 <- All_r_min1_min1 <- All_r_0_min1 <- NULL
All_R2H <- NULL
All_vector_p <- NULL

set.seed(Seed)
for (i in 1: M){
samples  <- sample(x=c(1:dim(Data)[1], n=dim(Data)[1]), replace = TRUE)  
Data_here <- Data[samples,]
Data_here <- Data_here[c(Surr, True, Treat)]
names(Data_here) <- c("Surr", "True", "Treat")

Probs <-  MarginalProbs(Dataset = Data_here, Surr = Surr, True = True,
              Treat = Treat)

MaxEntSPF <- MaxEntSPFBinBin(pi1_1_ = Probs$pi1_1_, pi1_0_ = Probs$pi1_0_, pi_1_1 = Probs$pi_1_1,
                             pi_1_0 = Probs$pi_1_0, pi0_1_ = Probs$pi0_1_, pi_0_1 = Probs$pi_0_1)
vector_p <- MaxEntSPF$Vector_p
All_vector_p <- rbind(All_vector_p, vector_p)


R2H <- MaxEntICABinBin(pi1_1_ = Probs$pi1_1_, pi1_0_ = Probs$pi1_0_, pi_1_1 = Probs$pi_1_1,
                pi_1_0 = Probs$pi_1_0, pi0_1_ = Probs$pi0_1_, pi_0_1 = Probs$pi_0_1)$R2_H
All_R2H <- c(All_R2H, R2H)


All_r_1_1 <- cbind(All_r_1_1, MaxEntSPF$r_1_1)
All_r_min1_1 <- cbind(All_r_min1_1, MaxEntSPF$r_min1_1)
All_r_0_1 <- cbind(All_r_0_1, MaxEntSPF$r_0_1)

All_r_1_0 <- cbind(All_r_1_0, MaxEntSPF$r_1_0)
All_r_min1_0 <- cbind(All_r_min1_0, MaxEntSPF$r_min1_0)
All_r_0_0 <- cbind(All_r_0_0, MaxEntSPF$r_0_0)

All_r_1_min1 <- cbind(All_r_1_min1, MaxEntSPF$r_1_min1)
All_r_min1_min1 <- cbind(All_r_min1_min1, MaxEntSPF$r_min1_min1)
All_r_0_min1 <- cbind(All_r_0_min1, MaxEntSPF$r_0_min1)

flush.console(); cat("\n", (i/M)*100, "% done", sep="")
}

fit <- list(R2H=All_R2H, r_1_1=All_r_1_1, r_min1_1=All_r_min1_1, r_0_1=All_r_0_1, 
            r_1_0=All_r_1_0, r_min1_0=All_r_min1_0, r_0_0=All_r_0_0, 
            r_1_min1=All_r_1_min1, r_min1_min1=All_r_min1_min1, r_0_min1=All_r_0_min1, 
            vector_p=All_vector_p, Call=match.call())
class(fit) <- "Bootstrap.MEP.BinBin"
fit

}






summary.Bootstrap.MEP.BinBin <- function(object, ..., Object){
  
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
      model_val <- c(0)
    }  
    fit <- list(mode_val= mode_val)  
  }
  
  suppressWarnings(warning("mode"))
  
  options(digits=5)
  
  try(Object$R2H <- na.exclude(Object$R2H), silent=TRUE)
  
  
  try(Object$r_1_1 <- na.exclude(Object$r_1_1), silent=TRUE)
  try(Object$r_min1_1 <- na.exclude(Object$r_min1_1), silent=TRUE)
  try(Object$r_0_1 <- na.exclude(Object$r_0_1), silent=TRUE)
  
  try(Object$r_1_0 <- na.exclude(Object$r_1_0), silent=TRUE)
  try(Object$r_min1_0 <- na.exclude(Object$r_min1_0), silent=TRUE)
  try(Object$r_0_0 <- na.exclude(Object$r_0_0), silent=TRUE)
  
  try(Object$r_1_min1 <- na.exclude(Object$r_1_min1), silent=TRUE)
  try(Object$r_min1_min1 <- na.exclude(Object$r_min1_min1), silent=TRUE)
  try(Object$r_0_min1 <- na.exclude(Object$r_0_min1), silent=TRUE)
  
  
  cat("\nResults bootstrap analyses:\n")
  cat("---------------------------\n\n")
  
  # R2H
  cat("\nR2_H results:\n")
  cat("-------------\n")
  
  cat("\nR2_H:\n", "Mean = ", mean(Object$R2H), ";  SD = ", sd(Object$R2H), 
          ";  95% CI = [", quantile(Object$R2H, probs = c(.025)), "; ",  
          quantile(Object$R2H, probs = c(.975)), "]\n",
          sep="")
  
  cat("\n\n\nSPF results:\n")
  cat("-------------\n")
  outcome_h <- Object$r_1_1
  cat("\nr_1_1:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
          ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
          quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  
  outcome_h <- Object$r_min1_1
  cat("\nr_min1_1:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
      ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
      quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  
  outcome_h <- Object$r_0_1
  cat("\nr_0_1:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
      ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
      quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  
  
  outcome_h <- Object$r_1_0
  cat("\nr_1_0:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
      ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
      quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  
  outcome_h <- Object$r_min1_0
  cat("\nr_min1_0:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
      ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
      quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  
  outcome_h <- Object$r_0_0
  cat("\nr_0_0:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
      ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
      quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  

  outcome_h <- Object$r_1_min1
  cat("\nr_1_min1:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
      ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
      quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  
  outcome_h <- Object$r_min1_min1
  cat("\nr_min1_min1:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
      ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
      quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  
  outcome_h <- Object$r_0_min1
  cat("\nr_0_min1:\n", "Mean: ", mean(outcome_h), ";  SD = ", sd(outcome_h), 
      ";  95% CI= [", quantile(outcome_h, probs = c(.025)), "; ",  
      quantile(outcome_h, probs = c(.975)), "]\n", sep="")
  

cat("\n\n\nMaximum entropy distribution vector of potential outcomes:\n")
cat("----------------------------------------------------------\n")
  
for (i in 1: dim(Object$vector_p)[2]){
  
  cat("\n", names(Object$vector_p)[i], ":\nMean: ", mean(Object$vector_p[,i]), ";  SD = ", sd(Object$vector_p[,i]), 
      ";  95% CI = [", quantile(Object$vector_p[,i], probs = c(.025)), "; ",  
      quantile(Object$vector_p[,i], probs = c(.975)), "]\n", sep="")

}

}
