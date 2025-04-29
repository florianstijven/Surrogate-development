ICA.Sample.ControlTreat <- function(T0S0, T1S1=seq(-1, 1, by = 0.001), T0T0 = 1, T1T1 = 1, S0S0 = 1, S1S1 = 1, 
          T0T1 = seq(-1, 1, by = 0.001), T0S1 = seq(-1, 1, by = 0.001), 
          T1S0 = seq(-1, 1, by = 0.001), S0S1 = seq(-1, 1, by = 0.001), 
          M = 50000, M.Target=NA) 
{
  T0S0_val <- T0S0
  T1S1_val <- T1S1
  T0T1_val <- T0T1
  T0S1_val <- T0S1
  T1S0_val <- T1S0
  S0S1_val <- S0S1
  # Variances
  T0T0_val <- T0T0
  T1T1_val <- T1T1
  S0S0_val <- S0S0
  S1S1_val <- S1S1
  Results <- na.exclude(matrix(NA, 1, 9))
  colnames(Results) <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", 
                         "S0S1", "ICA", "Sigma.Delta.T", "delta")
  
  All.Variances <- na.exclude(matrix(NA, 1, 4))
  colnames(All.Variances) <- c("T0T0", "T1T1", "S0S0", "S1S1")
  
  # Do M runs
  if (is.na(M.Target)==TRUE){
  for (i in 1:M) {
    T0T1 <- runif(n = 1, min = min(T0T1_val), max = max(T0T1_val))
    T0S0 <- runif(n = 1, min = min(T0S0_val), max = max(T0S0_val))
    T0S1 <- runif(n = 1, min = min(T0S1_val), max = max(T0S1_val))
    T1S0 <- runif(n = 1, min = min(T1S0_val), max = max(T1S0_val))
    T1S1 <- runif(n = 1, min = min(T1S1_val), max = max(T1S1_val))
    S0S1 <- runif(n = 1, min = min(S0S1_val), max = max(S0S1_val))
    Sigma_c <- diag(4)
    Sigma_c[2, 1] <- Sigma_c[1, 2] <- T0T1 * (sqrt(T0T0) * 
                                                sqrt(T1T1))
    Sigma_c[3, 1] <- Sigma_c[1, 3] <- T0S0 * (sqrt(T0T0) * 
                                                sqrt(S0S0))
    Sigma_c[4, 1] <- Sigma_c[1, 4] <- T0S1 * (sqrt(T0T0) * 
                                                sqrt(S1S1))
    Sigma_c[3, 2] <- Sigma_c[2, 3] <- T1S0 * (sqrt(T1T1) * 
                                                sqrt(S0S0))
    Sigma_c[4, 2] <- Sigma_c[2, 4] <- T1S1 * (sqrt(T1T1) * 
                                                sqrt(S1S1))
    Sigma_c[4, 3] <- Sigma_c[3, 4] <- S0S1 * (sqrt(S0S0) * 
                                                sqrt(S1S1))
    Sigma_c[1, 1] <- T0T0
    Sigma_c[2, 2] <- T1T1
    Sigma_c[3, 3] <- S0S0
    Sigma_c[4, 4] <- S1S1
    Cor_c <- cov2cor(Sigma_c)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)
    if (Min.Eigen.Cor > 0) {
      ICA <- ((sqrt(S0S0 * T0T0) * Cor_c[3, 1]) + (sqrt(S1S1 * 
                                                          T1T1) * Cor_c[4, 2]) - (sqrt(S0S0 * T1T1) * Cor_c[3, 
                                                                                                            2]) - (sqrt(S1S1 * T0T0) * Cor_c[4, 1]))/(sqrt((T0T0 + 
                                                                                                                                                              T1T1 - (2 * sqrt(T0T0 * T1T1) * Cor_c[2, 1])) * 
                                                                                                                                                             (S0S0 + S1S1 - (2 * sqrt(S0S0 * S1S1) * Cor_c[4, 
                                                                                                                                                                                                           3]))))
      if ((is.finite(ICA)) == TRUE) {
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0 * 
                                                   T1T1) * Cor_c[2, 1])
        delta <- sigma.delta.T * (1 - (ICA^2))
        results.part <- as.vector(cbind(T0T1, T0S0, T0S1, 
                                        T1S0, T1S1, S0S1, ICA, sigma.delta.T, delta))
        Results <- rbind(Results, results.part)
        rownames(Results) <- NULL
      }
    }
  }
  }
  
  # identify M ICA values  
  if (is.na(M.Target)==FALSE){
    count <- 0
    while (count < M.Target) {
      T0T1 <- runif(n = 1, min = min(T0T1_val), max = max(T0T1_val))
      T0S0 <- runif(n = 1, min = min(T0S0_val), max = max(T0S0_val))
      T0S1 <- runif(n = 1, min = min(T0S1_val), max = max(T0S1_val))
      T1S0 <- runif(n = 1, min = min(T1S0_val), max = max(T1S0_val))
      T1S1 <- runif(n = 1, min = min(T1S1_val), max = max(T1S1_val))
      S0S1 <- runif(n = 1, min = min(S0S1_val), max = max(S0S1_val))
      
      T0T0 <- runif(n = 1, min = min(T0T0_val), max = max(T0T0_val))
      T1T1 <- runif(n = 1, min = min(T1T1_val), max = max(T1T1_val))
      S0S0 <- runif(n = 1, min = min(S0S0_val), max = max(S0S0_val))
      S1S1 <- runif(n = 1, min = min(S1S1_val), max = max(S1S1_val))
      
      Sigma_c <- diag(4)
      Sigma_c[2, 1] <- Sigma_c[1, 2] <- T0T1 * (sqrt(T0T0) * 
                                                  sqrt(T1T1))
      Sigma_c[3, 1] <- Sigma_c[1, 3] <- T0S0 * (sqrt(T0T0) * 
                                                  sqrt(S0S0))
      Sigma_c[4, 1] <- Sigma_c[1, 4] <- T0S1 * (sqrt(T0T0) * 
                                                  sqrt(S1S1))
      Sigma_c[3, 2] <- Sigma_c[2, 3] <- T1S0 * (sqrt(T1T1) * 
                                                  sqrt(S0S0))
      Sigma_c[4, 2] <- Sigma_c[2, 4] <- T1S1 * (sqrt(T1T1) * 
                                                  sqrt(S1S1))
      Sigma_c[4, 3] <- Sigma_c[3, 4] <- S0S1 * (sqrt(S0S0) * 
                                                  sqrt(S1S1))
      Sigma_c[1, 1] <- T0T0
      Sigma_c[2, 2] <- T1T1
      Sigma_c[3, 3] <- S0S0
      Sigma_c[4, 4] <- S1S1
      Cor_c <- cov2cor(Sigma_c)
      Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE)
      if (Min.Eigen.Cor > 0) {
        ICA <- ((sqrt(S0S0 * T0T0) * Cor_c[3, 1]) + (sqrt(S1S1 * 
                                                            T1T1) * Cor_c[4, 2]) - (sqrt(S0S0 * T1T1) * Cor_c[3, 
                                                                                                              2]) - (sqrt(S1S1 * T0T0) * Cor_c[4, 1]))/(sqrt((T0T0 + 
                                                                                                                                                                T1T1 - (2 * sqrt(T0T0 * T1T1) * Cor_c[2, 1])) * 
                                                                                                                                                               (S0S0 + S1S1 - (2 * sqrt(S0S0 * S1S1) * Cor_c[4, 
                                                                                                                                                                                                             3]))))
        if ((is.finite(ICA)) == TRUE) {
          sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0 * 
                                                     T1T1) * Cor_c[2, 1])
          delta <- sigma.delta.T * (1 - (ICA^2))
          results.part <- as.vector(cbind(T0T1, T0S0, T0S1, 
                                          T1S0, T1S1, S0S1, ICA, sigma.delta.T, delta))
          Results <- rbind(Results, results.part)
          rownames(Results) <- NULL
          count <- count+1
          
          
          Variances.Here <- c(T0T0, T1T1, S0S0, S1S1)
          All.Variances <- rbind(All.Variances, Variances.Here)
          
        #  flush.console(); print(count)
        }
      }
    }
  }
  
  Results <- data.frame(Results, stringsAsFactors = TRUE)
  rownames(Results) <- NULL
  Total.Num.Matrices <- dim(Results)[1]
  fit <- list(Total.Num.Matrices = Total.Num.Matrices, Pos.Def = Results[, 
                                                                         1:6], ICA = Results$ICA, GoodSurr = Results[, 7:9], 
              Variances=data.frame(All.Variances),
              Call = match.call())
  class(fit) <- "ICA.ContCont"
  fit
}
