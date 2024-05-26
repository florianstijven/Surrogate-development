summary.SPF.BinCont <- function(object, ..., Object){
 
  options(digits = 4)
  
  if (missing(Object)){Object <- object} 
  
  Object$P_DT_0_DS_0 <- na.exclude(Object$P_DT_0_DS_0)
  Object$P_DT_psi_DS_max <- na.exclude(Object$P_DT_psi_DS_max)
  
  Object$r_min1_min1 <- na.exclude(Object$r_min1_min1)
  Object$r_0_min1 <- na.exclude(Object$r_0_min1)
  Object$r_1_min1 <- na.exclude(Object$r_1_min1)
  
  Object$r_min1_0 <- na.exclude(Object$r_min1_0)
  Object$r_0_0 <- na.exclude(Object$r_0_0)
  Object$r_1_0 <- na.exclude(Object$r_1_0)
  
  Object$r_min1_1 <- na.exclude(Object$r_min1_1)
  Object$r_0_1 <- na.exclude(Object$r_0_1)
  Object$r_1_1 <- na.exclude(Object$r_1_1)
  
  cat("\nFunction call:\n\n")
  print(Object$Call)

  cat("Requested interval Delta S: [",Object$a,",", Object$b,"]", sep="")
  cat("\n\n\n# SPF results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat("P(Delta T=0|Delta S=0)")
  cat("\nMean:", round(mean(Object$P_DT_0_DS_0), 4))
  cat("\nRange: [",round(range(Object$P_DT_0_DS_0)[1], 4),",", round(range(Object$P_DT_0_DS_0)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$P_DT_0_DS_0, probs=0.1), 4),",", round(quantile(Object$P_DT_0_DS_0, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$P_DT_0_DS_0, probs=0.025), 4),",", round(quantile(Object$P_DT_0_DS_0, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
  cat("P(Delta T = psi_ab(Delta S)")
  cat("\nMean:", round(mean(Object$P_DT_psi_DS_max), 4))
  cat("\nRange: [",round(range(Object$P_DT_psi_DS_max)[1], 4),",", round(range(Object$P_DT_psi_DS_max)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$P_DT_psi_DS_max, probs=0.1), 4),",", round(quantile(Object$P_DT_psi_DS_max, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$P_DT_psi_DS_max, probs=0.025), 4),",", round(quantile(Object$P_DT_psi_DS_max, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
  cat("r(-1,-1)")
  cat("\nMean:", round(mean(Object$r_min1_min1), 4))
  cat("\nRange: [",round(range(Object$r_min1_min1)[1], 4),",", round(range(Object$r_min1_min1)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_min1_min1, probs=0.1), 4),",", round(quantile(Object$r_min1_min1, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_min1_min1, probs=0.025), 4),",", round(quantile(Object$r_min1_min1, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
  cat("r(0,-1)")
  cat("\nMean:", round(mean(Object$r_0_min1), 4))
  cat("\nRange: [",round(range(Object$r_0_min1)[1], 4),",", round(range(Object$r_0_min1)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_0_min1, probs=0.1), 4),",", round(quantile(Object$r_0_min1, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_0_min1, probs=0.025), 4),",", round(quantile(Object$r_0_min1, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
  cat("r(1,-1)")
  cat("\nMean:", round(mean(Object$r_1_min1), 4))
  cat("\nRange: [",round(range(Object$r_1_min1)[1], 4),",", round(range(Object$r_1_min1)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_1_min1, probs=0.1), 4),",", round(quantile(Object$r_1_min1, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_1_min1, probs=0.025), 4),",", round(quantile(Object$r_1_min1, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
  cat("r(-1,0)")
  cat("\nMean:", round(mean(Object$r_min1_0), 4))
  cat("\nRange: [",round(range(Object$r_min1_0)[1], 4),",", round(range(Object$r_min1_0)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_min1_0, probs=0.1), 4),",", round(quantile(Object$r_min1_0, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_min1_0, probs=0.025), 4),",", round(quantile(Object$r_min1_0, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
  cat("r(0,0)")
  cat("\nMean:", round(mean(Object$r_0_0), 4))
  cat("\nRange: [",round(range(Object$r_0_0)[1], 4),",", round(range(Object$r_0_0)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_0_0, probs=0.1), 4),",", round(quantile(Object$r_0_0, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_0_0, probs=0.025), 4),",", round(quantile(Object$r_0_0, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
  cat("r(1,0)")
  cat("\nMean:", round(mean(Object$r_1_0), 4))
  cat("\nRange: [",round(range(Object$r_1_0)[1], 4),",", round(range(Object$r_1_0)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_1_0, probs=0.1), 4),",", round(quantile(Object$r_1_0, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_1_0, probs=0.025), 4),",", round(quantile(Object$r_1_0, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
  cat("r(-1,1)")
  cat("\nMean:", round(mean(Object$r_min1_1), 4))
  cat("\nRange: [",round(range(Object$r_min1_1)[1], 4),",", round(range(Object$r_min1_1)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_min1_1, probs=0.1), 4),",", round(quantile(Object$r_min1_1, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_min1_1, probs=0.025), 4),",", round(quantile(Object$r_min1_1, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")

  cat("r(0,1)")
  cat("\nMean:", round(mean(Object$r_0_1), 4))
  cat("\nRange: [",round(range(Object$r_0_1)[1], 4),",", round(range(Object$r_0_0)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_0_1, probs=0.1), 4),",", round(quantile(Object$r_0_1, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_0_1, probs=0.025), 4),",", round(quantile(Object$r_0_1, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")

  cat("r(1,1)")
  cat("\nMean:", round(mean(Object$r_1_1), 4))
  cat("\nRange: [",round(range(Object$r_1_1)[1], 4),",", round(range(Object$r_1_1)[2], 4),"]")
  cat("\nSDI0.80: [",round(quantile(Object$r_1_1, probs=0.1), 4),",", round(quantile(Object$r_1_1, probs=0.9), 4),"]")
  cat("\nSDI0.95: [", round(quantile(Object$r_1_1, probs=0.025), 4),",", round(quantile(Object$r_1_1, probs=0.975), 4),"]")
  cat("\n--------------------------------------------------\n\n")
  
}
  

