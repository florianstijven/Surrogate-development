ICA.ContCont.MultS_alt <- function(M = 500, N, Sigma,    
  G = seq(from=-1, to=1, by = .00001),
  Seed=c(123), Model = "Delta_T ~ Delta_S1 + Delta_S2", Show.Progress=FALSE){

SDs <- sqrt(diag(Sigma)); mu <- rep(0, times=length(SDs))
results_here <- results <- all_delta_S_T <- all_delta_S_T_here <- Lower.Dig.Corrs.All <- NULL

found <- 0
set.seed(Seed)
for (i in 1: M) {
  
  Sigma_c <- Sigma_c <- cov2cor(Sigma)
  num_elem <- dim(Sigma_c)[1] ** 2
  Sigma_c_orig <- Sigma_c
  size <- row_num <- col_num_ind <- 3; total_size <- dim(Sigma_c)[1]
  here <- Sigma_c
  
  while (size <= total_size) {
    here <- Sigma_c[(1:size), (1:size)] #LS!  
    here[is.na(here)] <- sample(x = G, size = length(here[is.na(here)]), replace = TRUE) #LS!
    here[upper.tri(here)] = t(here)[upper.tri(here)]  #LS!
    
      while (det(here) < 0){
      here <- Sigma_c[(1:size), (1:size)]  
      here[is.na(here)] <- sample(x = G, size = length(here[is.na(here)]), replace = TRUE)
      here[upper.tri(here)] = t(here)[upper.tri(here)]
    }
    
    Sigma_c[1:row_num, 1:col_num_ind] <- here 
    row_num <- row_num + 1; 
    col_num_ind <- col_num_ind + 1
    size <- size + 1
    }

  Sigma_c <- here
  
  Min.Eigen.Sigma <- try(min(eigen(Sigma_c)$values), TRUE)    # lowest eigenvalue

  if ((class(Min.Eigen.Sigma)!="try-error") & (Min.Eigen.Sigma >= 0.00000000001)) {
    found <- found + 1
    if (Show.Progress==TRUE){cat((i/M)*100, "% done... ", sep="")}
    corMat <- Sigma_c
    varVec <- SDs**2
    n = nrow(corMat)
    sdMat = diag(sqrt(varVec))
    rtn = sdMat %*% corMat %*% t(sdMat)
    Ideal <- MASS::mvrnorm(n = 10000, mu=mu, Sigma=rtn)  #T0, T1, S1_0, S1_1, etc
    
    # Make delta T and delta S1, S2 etc
    Delta_T <- Ideal[,2] - Ideal[,1]

    All_Delta_S <- NULL
    col_num <- 3
    for (k in 1: ((dim(Ideal)[2]-2)/2)){
      var_name <- c(1:10000)[k]
      assign(paste("Delta_S", var_name, sep=""), value =  Ideal[,(col_num + 1)] - Ideal[,col_num])
      Delta_S_hier <- assign(paste("Delta_S", var_name, sep=""), value =  Ideal[,(col_num + 1)] - Ideal[,col_num])
      col_num <- col_num + 2
      All_Delta_S <- cbind(All_Delta_S, Delta_S_hier)
    }

    Data_here <- cbind(Delta_T, All_Delta_S)
    Data_here <- data.frame(Data_here) 
    
    for (s in 1: (dim(Data_here)[2]-1)){
      names(Data_here)[s+1] <- paste("Delta_S", s, sep="")
    }   

    # Save lower diagonal Sigma, contains correlations
    Lower.Dig.Corrs.Here <- Sigma_c[lower.tri(Sigma_c)]
    
    
    # Fit models
    ICA <- sqrt(summary(lm(Model, data = Data_here))$r.squared)
    Adj.ICA2 <- (summary(lm(Model, data = Data_here))$adj.r.squared)
    if (Adj.ICA2>0){Adj.ICA <- sqrt(Adj.ICA2)}; if (Adj.ICA2<=0){Adj.ICA <- c(0)}
    res_error_Delta_T_given_delta_S <- (summary(lm(Model, data = Data_here)))$sigma # residual error delta T given delta S
    res_error_Delta_T <- (summary(lm(Delta_T~1, data = Data_here)))$sigma # residual error delta T 
    
    results_here <- 
      cbind(ICA, Adj.ICA, res_error_Delta_T_given_delta_S, res_error_Delta_T) 
    
    results <- rbind(results, results_here)
    
    Lower.Dig.Corrs.All <- data.frame(rbind(Lower.Dig.Corrs.All, Lower.Dig.Corrs.Here))
    row.names(Lower.Dig.Corrs.All) <- NULL
    
    all_delta_S_T_here <- cbind(i, ICA, Adj.ICA, Delta_T, All_Delta_S)
    all_delta_S_T <- rbind(all_delta_S_T, all_delta_S_T_here)  
  }  

  results <- data.frame(results)
  
  Sigma_c <- Sigma_c_orig  #LS!! 
  }

fit <- 
  list(R2_H=(as.numeric(results$ICA)**2), Corr.R2_H=(as.numeric(results$Adj.ICA)**2), Res_Err_Delta_T = as.numeric(results$res_error_Delta_T), 
       Res_Err_Delta_T_Given_S = as.numeric(results$res_error_Delta_T_given_delta_S), #All_Delta_S_T = all_delta_S_T, 
       Lower.Dig.Corrs.Sigma=Lower.Dig.Corrs.All, 
       Call=match.call()) 

class(fit) <- "ICA.ContCont.MultS"
fit

}

