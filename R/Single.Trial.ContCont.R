
# Dataset = ARMD_S
# attach(ARMD_S)
# Surr = Diff24
# True = Diff52
# Treat = Treat
# Pat.ID = Id
# Trial.ID <- rep(1, nrow(Dataset))
# Alpha=0.05

# data(Schizo)
# Dataset=Schizo
# attach(Schizo)
# attach(Dataset)
# Surr = BPRS
# True = PANSS
# Treat = Treat
# Pat.ID = Id
# Alpha=0.05



Single.Trial.ContCont <- function(Dataset, Surr, True, Treat, Pat.ID, Alpha=.05, 
                                  Number.Bootstraps=500, Seed=12345){
  
  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Trial.ID <- rep(1, nrow(Dataset))
  Pat.ID <- Dataset[,paste(substitute(Pat.ID))]
  
  Data.Proc <- .Data.Processing(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat, Trial.ID=Trial.ID, Pat.ID=Pat.ID, Min.Trial.Size=0)
  wide <- Data.Proc$wide
  dataS <- Data.Proc$dataS
  dataT <- Data.Proc$dataT
  Data.analyze <- Data.Proc$Data.analyze
  N.total <- Data.Proc$N.total
  
  
  # Prentice
  model12 <- lm(cbind(wide$Surr, wide$True)~ wide$Treat, data=wide)
  alpha <- model12$coefficients[2,1]
  beta <- model12$coefficients[2,2]
  model3 <- lm(wide$True ~ wide$Surr, data=wide)
  gamma <- model3$coefficients[2]
  model4 <- lm(wide$True ~ wide$Treat + wide$Surr, data=wide)
  beta_s <-  model4$coefficients[2]
  
  # p values for the relevant parameters
  P_model1 <- summary(model12)$"Response wide$Surr"$coefficients
  rownames(P_model1) <- c("Intercept", "Treatment")
  P_model2 <- summary(model12)$"Response wide$True"$coefficients
  rownames(P_model2) <- c("Intercept", "Treatment")
  P_model3 <- summary(model3)$coefficients
  rownames(P_model3) <- c("Intercept", "Surrogate")
  P_model4 <- summary(model4)$coefficients
  rownames(P_model4) <- c("Intercept", "Treatment", "Surrogate")
  
  P_crit1 <- summary(model12)$"Response wide$Surr"$coefficients[2, 4]
  P_crit2 <- summary(model12)$"Response wide$True"$coefficients[2, 4]
  P_crit3 <- summary(model3)$coefficients[2,4]
  P_crit4 <- data.frame(summary(model4)$coefficients, stringsAsFactors = TRUE)[2, 4]
  
  if ((P_crit1 < Alpha & P_crit2 < Alpha & P_crit3 < Alpha & P_crit4 > Alpha)==TRUE) {Prentice.Passed <- TRUE}
  if ((P_crit1 < Alpha & P_crit2 < Alpha & P_crit3 < Alpha & P_crit4 > Alpha)==FALSE) {Prentice.Passed <- FALSE}
  
  # Beta_s CI
  beta_s <-  model4$coefficients[2]
  beta_s_se <- summary(model4)$coefficients[2,2]
  beta_s_lb <- beta_s - (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*beta_s_se)
  beta_s_ub <- beta_s + (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*beta_s_se)
  beta_s_results <- data.frame(cbind(beta_s, beta_s_se , beta_s_lb, beta_s_ub), stringsAsFactors = TRUE)
  colnames(beta_s_results) <- c("beta_s", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(beta_s_results) <- c(" ")
  
  # alpha CI
  alpha <- model12$coefficients[2,1]
  alpha_se <- summary(model12)$"Response wide$Surr"$coefficients[2,2]
  alpha_lb <- alpha - (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*alpha_se)
  alpha_ub <- alpha + (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*alpha_se)
  alpha_results <- data.frame(cbind(alpha, alpha_se , alpha_lb, alpha_ub), stringsAsFactors = TRUE)
  colnames(alpha_results) <- c("Alpha", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(alpha_results) <- c(" ")
  
  # beta CI
  model12 <- lm(cbind(wide$Surr, wide$True) ~ wide$Treat, data=wide)
  Residuals <- data.frame(model12$residuals, stringsAsFactors = TRUE)
  colnames(Residuals) <- c("Surr", "True")
  beta <- model12$coefficients[2,2]
  beta_se <- summary(model12)$"Response wide$True"$coefficients[2,2]
  beta_lb <- beta - (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*beta_se)
  beta_ub <- beta + (qt(c(1-Alpha/2), df=N.total, lower.tail=TRUE)*beta_se)
  beta_results <- data.frame(cbind(beta, beta_se , beta_lb, beta_ub), stringsAsFactors = TRUE)
  colnames(beta_results) <- c("Beta", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(beta_results) <- c(" ")
  
  
  ## AA and RE

  # CI AA
  rho_z <- var(model12$residuals[,2], model12$residuals[,1])/ sqrt(var(model12$residuals[,2])*var(model12$residuals[,1]))
  Z <- .5*log((1+rho_z)/(1-rho_z))
  rho_lb <- max((exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_ub <- min(1, (exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_sd <- sqrt((1-rho_z**2)/(N.total-2))
  rho_results_FishZ <- data.frame(cbind(rho_z, rho_sd , rho_lb, rho_ub), stringsAsFactors = TRUE)
  colnames(rho_results_FishZ) <- c("AA (gamma)", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(rho_results_FishZ) <- c(" ")
  
  # RE 
  model12 <- lm(cbind(wide$Surr, wide$True) ~ wide$Treat, data=wide)
  alpha <- model12$coefficients[2,1]  
  beta <- model12$coefficients[2,2]
  RE <- beta/alpha
  
  # Delta method
  X <- model.matrix(~ wide$Treat)   
  df_res <- model12$df.residual     # = n - p
  resid_mat <- residuals(model12)   # matrix with two columns
      # Residual covariance matrix Sigma_hat (2 x 2)
  Sigma_hat <- crossprod(resid_mat) / df_res   # same as t(resid_mat) %*% resid_mat / df_res
  # (X'X)^{-1} and scalar for Treat row (row/col 2)
  XtX_inv <- solve(t(X) %*% X)
  v_scalar <- XtX_inv[2, 2]   # element for Treat coefficient
  # Variances and covariance of alpha and beta
  Var_alpha <- v_scalar * Sigma_hat[1, 1]
  Var_beta  <- v_scalar * Sigma_hat[2, 2]
  Cov_ab    <- v_scalar * Sigma_hat[1, 2]
  var_RE_delta <- Var_beta / (alpha^2) + (beta^2) * Var_alpha / (alpha^4) - 2 * beta * Cov_ab / (alpha^3)
  se_RE_delta  <- sqrt(var_RE_delta)
  tval <- qt(1-Alpha/2, df = df_res)
  CI_delta <- c(RE - tval * se_RE_delta, RE + tval * se_RE_delta)
  
  RE_results_Delta <- data.frame(cbind(RE, se_RE_delta, as.numeric(CI_delta[1]), as.numeric(CI_delta[2])), stringsAsFactors = TRUE)
  colnames(RE_results_Delta) <- c("RE", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(RE_results_Delta) <- c(" ")
  
  # Bootstrap CI
  d.size <- dim(wide)[1]
  obs <- c(1:d.size)
  k <- Number.Bootstraps
  RE_boot <- rho_z_boot <- alpha_boot <- beta_boot <- as.vector(NULL)
  for (i in 1:k){
    set.seed(Seed+i)
    index <- sample(obs, d.size, replace=TRUE)
    sample <- data.frame(wide[index,], stringsAsFactors = TRUE)
    sample <- na.exclude(sample[order(sample$Pat.ID),])
    model1 <- lm(cbind(sample$Surr, sample$True)~ sample$Treat, data=sample)
    alpha_boot[i] <- model1$coefficients[2,1]
    beta_boot[i] <- model1$coefficients[2,2]
    rho_z_boot[i] <- var(model1$residuals[,2], model1$residuals[,1])/ sqrt(var(model1$residuals[,2])*var(model1$residuals[,1]))
    RE_boot[i] <- beta_boot[i] / alpha_boot[i]
  }
  RE_CIs <- quantile(RE_boot, probs=c(Alpha/2, 1-Alpha/2), na.rm = TRUE)
  RE_results_Boot <- data.frame(cbind(RE, sd(RE_boot), RE_CIs[1], RE_CIs[2]), stringsAsFactors = TRUE)
  colnames(RE_results_Boot) <- c("RE", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(RE_results_Boot) <- c(" ")
  
  rho_CIs <- quantile(rho_z_boot, probs=c(Alpha/2, 1-Alpha/2), na.rm = TRUE)
  rho_results_Boot <- data.frame(cbind(rho_z, sd(rho_z_boot), rho_CIs[1], rho_CIs[2]), stringsAsFactors = TRUE)
  colnames(rho_results_Boot) <- c("AA (gamma)", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(rho_results_Boot) <- c(" ")
  
  NoTreat <- wide[wide$Treat!=1,]
  Treat <- wide[wide$Treat==1,]
  T0S0 <- cor(NoTreat$Surr, NoTreat$True)
  T1S1 <- cor(Treat$Surr, Treat$True)
  Z_T0S0 <- .5*log((1+T0S0)/(1-T0S0))
  rho_lb <- max(0, (exp(2*(Z_T0S0-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z_T0S0-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_ub <- min(1, (exp(2*(Z_T0S0+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z_T0S0+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_sd <- sqrt((1-T0S0**2)/(N.total-2))
  rho_results_T0S0 <- data.frame(cbind(T0S0, rho_sd , rho_lb, rho_ub), stringsAsFactors = TRUE)
  colnames(rho_results_T0S0) <- c("Estimate", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(rho_results_T0S0) <- c(" ")
  Z_T1S1 <- .5*log((1+T1S1)/(1-T1S1))
  rho_lb <- max(0, (exp(2*(Z_T1S1-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z_T1S1-(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_ub <- min(1, (exp(2*(Z_T1S1+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))-1)/(exp(2*(Z_T1S1+(qnorm(1-Alpha/2)*sqrt(1/(N.total-3)))))+1))
  rho_sd <- sqrt((1-T1S1**2)/(N.total-2))
  rho_results_T1S1 <- data.frame(cbind(T1S1, rho_sd , rho_lb, rho_ub), stringsAsFactors = TRUE)
  colnames(rho_results_T1S1) <- c("Estimate", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(rho_results_T1S1) <- c(" ")
  Cor.Endpoints <- data.frame(rbind(rho_results_T0S0, rho_results_T1S1), stringsAsFactors = TRUE)
  rownames(Cor.Endpoints) <- c("r_T0S0", "r_T1S1")
  colnames(Cor.Endpoints) <- c("Estimate", "Standard Error", "CI lower limit", "CI upper limit")
  
  # Fieler
  model12 <- lm(cbind(Surr, True) ~ Treat, data = wide)
      # Extract the full vector of coefficients and the VCV matrix
  beta_hat <- as.vector(coef(model12))
  Sigma <- crossprod(residuals(model12)) / model12$df.residual
  invXtX <- solve(crossprod(model.matrix(model12)))
  V_full <- kronecker(Sigma, invXtX)
  # Numerator contrast
  numC <- matrix(c(0, 0, 0, 1), nrow = 1)
  # Denominator contrast
  denC <- matrix(c(0, 1, 0, 0), nrow = 1)
  fieller_re <- gsci.ratio(est = beta_hat, vcmat =  V_full, Num.Contrast = numC, Den.Contrast = denC,
                           conf.level = 1-Alpha, adjusted = FALSE)
  
  RE_results_Fieller <- data.frame(cbind(fieller_re$estimate, fieller_re$conf.int[1], fieller_re$conf.int[2]), stringsAsFactors = TRUE)
  colnames(RE_results_Fieller) <- c("RE", "CI lower limit", "CI upper limit")
  rownames(RE_results_Fieller) <- c(" ")
  
  
  
  # Proportion explained
  # ~~~~~~~~~~~~~~~~~~~~
  model1 <- lm(wide$True ~ wide$Treat) 
  model2 <- lm(wide$True ~ wide$Treat + wide$Surr)
  beta <- model1$coefficients[2]
  beta_s <- model2$coefficients[2]
  PE <- 1-(beta_s / beta)
  
  # Delta method CI 
  X1 <- model.matrix(model1)   # design matrices n x p1
  X2 <- model.matrix(model2)   # n x p2
       # compute coefficient-linear-combination vectors:
  # coefficients = C %*% Y  with  C = (X'X)^{-1} X'
  C1 <- solve(t(X1) %*% X1) %*% t(X1)   # p1 x n
  C2 <- solve(t(X2) %*% X2) %*% t(X2)   # p2 x n
     # row for Treat is the 2nd row (intercept at row 1)
  a <- C1[2, ]   # length n, gives beta = a %*% Y
  b <- C2[2, ]   # length n, gives beta_s = b %*% Y
      # homoskedastic delta-method (use residual variance from full/adjusted model)
  sigma2_hat <- sum(resid(model2)^2) / model2$df.residual
  V <- matrix(0, 2, 2)
  V[1,1] <- sigma2_hat * sum(a * a)    # Var(beta)
  V[2,2] <- sigma2_hat * sum(b * b)    # Var(beta_s)
  V[1,2] <- V[2,1] <- sigma2_hat * sum(a * b)  # Cov(beta, beta_s)
  # gradient of g(beta, beta_s) = 1 - beta_s / beta
  # dg/dbeta   = beta_s / beta^2
  # dg/dbeta_s = - 1 / beta
  g <- c(beta_s / (beta^2), -1 / beta)
  var_PE <- as.numeric(t(g) %*% V %*% g)
  se_PE  <- sqrt(var_PE)
  z <- qnorm(1-Alpha/2)
  CI_lower <- PE - z * se_PE
  CI_upper <- PE + z * se_PE
  PE_results_Delta <- data.frame(cbind(PE, se_PE, CI_lower, CI_upper))   
  colnames(PE_results_Delta) <- c("PE", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(PE_results_Delta) <- c(" ")  
  
  
  # Bootstrap CI
  d.size <- dim(wide)[1]
  obs <- c(1:d.size)
  k <- Number.Bootstraps
  PE_boot <- as.vector(NULL)
  for (i in 1:k){
    set.seed(Seed+i)
    index <- sample(obs, d.size, replace=TRUE)
    sample <- data.frame(wide[index,], stringsAsFactors = TRUE)
    sample <- na.exclude(sample[order(sample$Pat.ID),])
    
    sample.model1 <- lm(sample$True ~ sample$Treat) 
    sample.model2 <- lm(sample$True ~ sample$Treat + sample$Surr)
    sample.beta <- sample.model1$coefficients[2]
    sample.beta_s <- sample.model2$coefficients[2]
    sample.PE <- 1-(sample.beta_s / sample.beta)
    PE_boot[i] <- sample.PE
  }
  PE_CIs <- quantile(PE_boot, probs=c(Alpha/2, 1-Alpha/2), na.rm = TRUE)
  PE_results_Boot <- data.frame(cbind(PE, sd(PE_boot), PE_CIs[1], PE_CIs[2]), stringsAsFactors = TRUE)
  colnames(PE_results_Boot) <- c("PE", "Standard Error", "CI lower limit", "CI upper limit")
  rownames(PE_results_Boot) <- c(" ")
  
  
  # Fieller 
    # Get coefficients and vcov from each model
  b1 <- coef(model1)
  b2 <- coef(model2)
  v1 <- vcov(model1)
  v2 <- vcov(model2)
    # Create a single vector of all coefficients
  beta_hat <- c(b1, b2)
  # Create a combined block-diagonal variance-covariance matrix
  V_combined <- as.matrix(Matrix::bdiag(v1, v2))
  colnames(V_combined) <- rownames(V_combined) <- names(beta_hat)
  # Define the ratio of interest (beta_s / beta) using contrast matrices
  # Numerator contrast: selects beta_s (model2$coefficients[2])
  numC <- matrix(c(0, 0, 0, 1, 0), nrow = 1)
  # Denominator contrast: selects beta (model1$coefficients[2])
  denC <- matrix(c(0, 1, 0, 0, 0), nrow = 1)
  # Calculate Fieller's CI for the ratio (theta = beta_s / beta)
  fieller_ratio <- gsci.ratio(
    est = beta_hat, vcmat =  V_combined, Num.Contrast = numC, Den.Contrast = denC, conf.level = 1-Alpha, adjusted = FALSE)
  # Transform the CI for PE = 1 - (beta_s / beta)
  PE <- 1 - fieller_ratio$estimate
  ci_ratio <- fieller_ratio$conf.int
  PE_ci_lower <- 1 - ci_ratio[2] # 1 - Upper bound of ratio CI
  PE_ci_upper <- 1 - ci_ratio[1] # 1 - Lower bound of ratio CI
  PE_results_Fieller <- data.frame(cbind(PE, PE_ci_lower, PE_ci_upper), stringsAsFactors = TRUE)
  colnames(PE_results_Fieller) <- c("PE", "CI lower limit", "CI upper limit")
  rownames(PE_results_Fieller) <- c(" ")
  
  
  
  fit <-
    list(Data.Analyze=wide, 
      Prentice.Model.1=P_model1, Prentice.Model.2=P_model2, Prentice.Model.3=P_model3, Prentice.Model.4=P_model4,
         Prentice.Passed=Prentice.Passed, 
         Alpha=alpha_results, Beta=beta_results, RE.Delta=RE_results_Delta, RE.Boot=RE_results_Boot, RE.Boot.Samples=RE_boot, AA=rho_results_FishZ, AA.Boot=rho_results_Boot, AA.Boot.Samples=rho_z_boot,
         Cor.Endpoints=Cor.Endpoints, Residuals=Residuals, 
      Beta_S = beta_s_results, PE.Delta = PE_results_Delta, PE.Boot=PE_results_Boot, PE.Boot.Samples=PE_boot,
      PE.Fieller = PE_results_Fieller, RE.Fieller = RE_results_Fieller,
      Call=match.call())
  
  class(fit) <- "Single.Trial.ContCont"
  fit
  
}

