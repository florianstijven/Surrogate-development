ISTE.ContCont <- function(Mean_T1, Mean_T0, Mean_S1, Mean_S0, N, Delta_S=c(-10, 0, 10), zeta.PI=0.05, 
  PI.Bound=0, PI.Lower=TRUE, Show.Prediction.Plots=TRUE, Save.Plots="No", 
  T0S0, T1S1, T0T0=1, T1T1=1, S0S0=1, S1S1=1, T0T1=seq(-1, 1, by=.001), T0S1=seq(-1, 1, by=.001), 
  T1S0=seq(-1, 1, by=.001), S0S1=seq(-1, 1, by=.001), M.PosDef=500, Seed=123) { 
  
  total <- 0  # For M.PosDef counter
  set.seed(Seed)
  T0S0_val <- T0S0
  T1S1_val <- T1S1 
  T0T1_val <- T0T1
  T0S1_val <- T0S1
  T1S0_val <- T1S0
  S0S1_val <- S0S1
  Mean_T1_val <- Mean_T1
  Mean_T0_val <- Mean_T0
  Mean_S1_val <- Mean_S1
  Mean_S0_val <- Mean_S0
  
  T0T0_val <- T0T0
  T1T1_val <- T1T1
  S0S0_val <- S0S0
  S1S1_val <- S1S1
  
  
  # To save ICA results
  Results.ICA <- na.exclude(matrix(NA, 1, 9))  
  colnames(Results.ICA) <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "ICA", "Sigma.Delta.T", "delta") 
  
  # Variables to save ISTE results
  STE_Low_PI_all <- STE_Up_PI_all <- MSE_all <- b0_all <- b1_all <- Exp_Delta_S_Delta_T_equal_O_all <- S_squared_pred_all <- 
    Predicted_Delta_T_all <- PI_Interval_Low_all <- PI_Interval_Up_all <- STE_Low_PI_0_all <- STE_Up_PI_0_all <- NULL
  T0T0_all <- T1T1_all <- S0S0_all <- S1S1_all <- Mean_DeltaT_all <- Mean_DeltaS_all <- NULL 
  
  
  while (total < M.PosDef){
    # vals for ISTE
    Mean_T1 <- sample(Mean_T1_val, size = 1)
    Mean_T0 <- sample(Mean_T0_val, size = 1)
    Mean_S1 <- sample(Mean_S1_val, size = 1)
    Mean_S0 <- sample(Mean_S0_val, size = 1)
   
    # vals for ISTE and ICA
    T0T0 <- sample(T0T0_val, size=1)
    T1T1 <- sample(T1T1_val, size=1)
    S0S0 <- sample(S0S0_val, size=1)
    S1S1 <- sample(S1S1_val, size=1)
    
    # computations for ICA
    T0T1 <- runif(n = 1, min = min(T0T1_val), max = max(T0T1_val))
    T0S0 <- runif(n = 1, min = min(T0S0_val), max = max(T0S0_val)) 
    T0S1 <- runif(n = 1, min = min(T0S1_val), max = max(T0S1_val))  
    T1S0 <- runif(n = 1, min = min(T1S0_val), max = max(T1S0_val))  
    T1S1 <- runif(n = 1, min = min(T1S1_val), max = max(T1S1_val))  
    S0S1 <- runif(n = 1, min = min(S0S1_val), max = max(S0S1_val)) 
    Sigma_c <- diag(4)         
    Sigma_c[2,1] <- Sigma_c[1,2] <- T0T1 * (sqrt(T0T0)*sqrt(T1T1))
    Sigma_c[3,1] <- Sigma_c[1,3] <- T0S0 * (sqrt(T0T0)*sqrt(S0S0))
    Sigma_c[4,1] <- Sigma_c[1,4] <- T0S1 * (sqrt(T0T0)*sqrt(S1S1))
    Sigma_c[3,2] <- Sigma_c[2,3] <- T1S0 * (sqrt(T1T1)*sqrt(S0S0))
    Sigma_c[4,2] <- Sigma_c[2,4] <- T1S1 * (sqrt(T1T1)*sqrt(S1S1))
    Sigma_c[4,3] <- Sigma_c[3,4] <- S0S1 * (sqrt(S0S0)*sqrt(S1S1))
    Sigma_c[1,1] <- T0T0
    Sigma_c[2,2] <- T1T1
    Sigma_c[3,3] <- S0S0
    Sigma_c[4,4] <- S1S1
    Cor_c <- cov2cor(Sigma_c)
    Min.Eigen.Cor <- try(min(eigen(Cor_c)$values), TRUE) 

    if (Min.Eigen.Cor > 0) {   # Pos Def Sigma found
      
      total <- total + 1  # increase counter pos def sigma found 
      cat(paste(total/M.PosDef*100, "% done...", sep=""))
      
      # Computations ICA
      ICA <- ((sqrt(S0S0*T0T0)*Cor_c[3,1])+(sqrt(S1S1*T1T1)*Cor_c[4,2])-(sqrt(S0S0*T1T1)*Cor_c[3,2])-(sqrt(S1S1*T0T0)*Cor_c[4,1]))/(sqrt((T0T0+T1T1-(2*sqrt(T0T0*T1T1)*Cor_c[2,1]))*(S0S0+S1S1-(2*sqrt(S0S0*S1S1)*Cor_c[4,3]))))
      if ((is.finite(ICA))==TRUE){
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * Cor_c[2,1])
        delta <- sigma.delta.T * (1-(ICA**2))
        Results.ICA.part <- as.vector(cbind(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1, ICA, sigma.delta.T, delta))
        Results.ICA <- rbind(Results.ICA, Results.ICA.part)
        rownames(Results.ICA) <- NULL}
      
      
      # Compuations ISTE
      Mean_DeltaT <- Mean_T1 - Mean_T0
      Mean_DeltaS <- Mean_S1 - Mean_S0
      
      Corr_mat_hier <- diag(c(1, 1, 1, 1))
      Corr_mat_hier[2,1] <- Corr_mat_hier[1,2] <- T0T1
      Corr_mat_hier[3,1] <- Corr_mat_hier[1,3] <- T0S0
      Corr_mat_hier[4,1] <- Corr_mat_hier[1,4] <- T0S1
      Corr_mat_hier[3,2] <- Corr_mat_hier[2,3] <- T1S0
      Corr_mat_hier[4,2] <- Corr_mat_hier[2,4] <- T1S1
      Corr_mat_hier[4,3] <- Corr_mat_hier[3,4] <- S0S1
      
      stdevs <- c(sqrt(T0T0), sqrt(T1T1), sqrt(S0S0), sqrt(S1S1))  # standard deviations
      b <- stdevs %*% t(stdevs)  
      Sigma_Matrix_hier <- b * Corr_mat_hier
      
      # covariances and variances Delta S and Delta T
      cov_deltaT_deltaS <- Sigma_Matrix_hier[4,2] - Sigma_Matrix_hier[3,2] - Sigma_Matrix_hier[4,1] + Sigma_Matrix_hier[3,1]
      var_delta_S <- Sigma_Matrix_hier[4,4] + Sigma_Matrix_hier[3,3] - (2 * Sigma_Matrix_hier[4,3])
      var_delta_T <- Sigma_Matrix_hier[2,2] + Sigma_Matrix_hier[1,1] - (2 * Sigma_Matrix_hier[2,1])
      
      # beta 0 en 1 of simple regression model regressing Delta T on Delta S
      b1 <- cov_deltaT_deltaS / var_delta_S
      b0 <- Mean_DeltaT - b1 * Mean_DeltaS
      
      # Predicted Delta_T
      Predicted_Delta_T <- b0 + b1 * Delta_S
      Exp_Delta_S_Delta_T_equal_O <- -b0 / b1
      
      # Predict Delta T
      MSE <- ((N-1)/(N-2)) * (var_delta_T - ((cov_deltaT_deltaS**2) / var_delta_S)) #; MSE 
      
      #var_delta_T * (1 - ICA.ContCont.Fitted$ICA[1]**2). p59 K
      t_val <- qt(c(1-zeta.PI/2), df=N-2)
      S_squared_pred <- MSE * (1 + (1/N) + (((Delta_S - Mean_DeltaS)**2)/(var_delta_S*(N-1))))
      
      PI_Interval_Low <- (b0 + b1 * Delta_S) - (t_val * sqrt(S_squared_pred)) #; PI_interval_min
      PI_Interval_Up <- (b0 + b1 * Delta_S) + (t_val * sqrt(S_squared_pred)) #; PI_interval_max
      
      # compute Delta S which has lower and upper PI equal to 0
      A <- b0
      B <- b1
      M <- MSE
      V <- var_delta_S
      G <- Mean_DeltaS
      T2 <- t_val**2 
      
      
      STE_Low_PI_0 <- STE_Up_PI_0 <- NA
      options(warn=-1)
      # Compute Delta S value for which lower bound of PI around predicted Delta T = 0 
      try(STE_Low_PI_0 <- (-sqrt((-2 * A * B * N * V + 2 * A * B * V - 2 * G * M *T2)**2 - 4 * (B**2 * -N * V + B**2 * V + M * T2)*
                               (A**2 * -N * V + A**2 * V + G**2 * M * T2 + M * N * T2 * V - (M * T2*V)/N))+
                         2 * A * B * N * V - 2 * A * B * V + 2 * G * M *T2)/
        (2 * (B**2 * -N * V + B**2 * V + M *T2)), silent = TRUE)
      
      
      # Compute Delta S value for which upper bound of PI around predicted Delta T = 0. Only differce wth above equation = change minus sign in front of sqrt 
      try(STE_Up_PI_0 <- (sqrt((-2 * A * B * N * V + 2 * A * B * V - 2 * G * M *T2)**2 - 4 * (B**2 * -N * V + B**2 * V + M * T2)*
                             (A**2 * -N * V + A**2 * V + G**2 * M * T2 + M * N * T2 * V - (M * T2*V)/N))+
                        2 * A * B * N * V - 2 * A * B * V + 2 * G * M *T2)/
        (2 * (B**2 * -N * V + B**2 * V + M *T2)), silent=TRUE)
      options(warn=0)
      
      # below is code to consider other values than 0 for bounds of PI around delta T
      # Function to compute Delta S value for which lower bound of PI around predicted Delta T = 0 (or other value specified by PI.Bound)
      STE.fun.Low <- function(X,A=b0, B=b1, M=MSE, V=var_delta_S, G=Mean_DeltaS) {
        (A+B*X)-t_val*sqrt(MSE *(1+1/N+((X-G)**2/(V*(N-1)))  
        ))
      }
      # Function to compute  Delta S value for which upper bound of PI around predicted Delta T = 0 (or other value specified by PI.Bound)
      STE.fun.Up <- function(X,A=b0, B=b1, M=MSE, V=var_delta_S, G=Mean_DeltaS) {
        (A+B*X)+t_val*sqrt(MSE *(1+1/N+((X-G)**2/(V*(N-1)))  
        ))
      }
      u1 <- u2 <- NA
      try(u1 <- uniroot(function(x) STE.fun.Low(x)-PI.Bound,c(-100000000,100000000), extendInt = "yes"), silent=TRUE)
      try(u2 <- uniroot(function(x) STE.fun.Up(x)-PI.Bound,c(-100000000,100000000), extendInt = "yes"), silent=TRUE)
      
      # Delta S value for which lower bound of PI around predicted Delta T = 0 (or other value specified by PI.Bound)
      if (is.na(u1[1])==FALSE){STE_Low_PI <- u1$root}
      if (is.na(u1[1])==TRUE){STE_Low_PI <- NA}
      # Delta S value for which upper bound of PI around predicted Delta T = 0 (or other value specified by PI.Bound)
      if (is.na(u2[1])==FALSE){STE_Up_PI <- u2$root}
      if (is.na(u2[1])==TRUE){STE_Up_PI <- NA}
      
      
      
      # plot
      if (is.na(STE_Low_PI[1])==FALSE & is.na(PI_Interval_Up[1])==FALSE & 
          is.na(PI_Interval_Low[1])==FALSE & is.na(STE_Low_PI[1])==FALSE &
        Show.Prediction.Plots==TRUE){
        max_val_y <- max(PI_Interval_Up, STE_Low_PI)
        min_val_y <- min(PI_Interval_Low, STE_Low_PI)
        plot(main=paste("Run", total), x=Delta_S, y=(b0 + b1 * Delta_S), 
             xlab=expression(paste(Delta, S[0])), 
             ylab=expression(paste("E(", Delta, T, "|", Delta, S[0], ")")), type="l", 
             ylim=c(min_val_y, max_val_y), 
         #    xlim=c(STE_Low_PI*-1.5, STE_Low_PI*1.5),
             lwd=2)
        lines(x=Delta_S, y=PI_Interval_Low, col=1, lty=2)
        lines(x=Delta_S, y=PI_Interval_Up, col=1, lty=2)
        abline(h=PI.Bound, col="grey")
        arrows(x0 = STE_Low_PI, y0 = PI.Bound, x1 = STE_Low_PI, y1 = min_val_y, lwd=2)
      #  abline(v=STE_Low_PI, col="red")
      }
      
      if (Save.Plots!="No"){
        options(warn=-1)
        dir.create(Save.Plots)
        options(warn=0)
        
        pdf(paste(Save.Plots, "/Run", total, ".pdf", sep=""),width=6,height=4,paper='special') 
        max_val_y <- max(PI_Interval_Up, STE_Low_PI)
        min_val_y <- min(PI_Interval_Low, STE_Low_PI)
        plot(main=paste("Run", total), x=Delta_S, y=(b0 + b1 * Delta_S), xlab=expression(paste(Delta, S)), 
             ylab=expression(paste("E(", Delta, T, ")")), type="l", ylim=c(min_val_y, max_val_y))
        lines(x=Delta_S, y=PI_Interval_Low, col="grey")
        lines(x=Delta_S, y=PI_Interval_Up, col="grey")
        abline(v=STE_Low_PI, col="red")
        abline(h=PI.Bound, col="red")
        dev.off()
      }
      
      
      # save results
      STE_Low_PI_0_all <- c(STE_Low_PI_0_all, STE_Low_PI_0)
      STE_Up_PI_0_all <- c(STE_Up_PI_0_all, STE_Up_PI_0)
      STE_Low_PI_all <- c(STE_Low_PI_all, STE_Low_PI)
      STE_Up_PI_all <- c(STE_Up_PI_all, STE_Up_PI)
      MSE_all <- c(MSE_all, MSE) 
      b0_all <- c(b0_all, b0)
      b1_all <- c(b1_all, b1)
      Exp_Delta_S_Delta_T_equal_O_all <- c(Exp_Delta_S_Delta_T_equal_O_all, Exp_Delta_S_Delta_T_equal_O)
      S_squared_pred_all <- c(S_squared_pred_all, S_squared_pred)
      Predicted_Delta_T_all <- rbind(Predicted_Delta_T_all, Predicted_Delta_T)
      PI_Interval_Low_all <- rbind(PI_Interval_Low_all, PI_Interval_Low)
      PI_Interval_Up_all <- rbind(PI_Interval_Up_all, PI_Interval_Up)
      T0T0_all <- c(T0T0_all, T0T0) 
      T1T1_all <- c(T1T1_all, T1T1) 
      S0S0_all <- c(S0S0_all, S0S0) 
      S1S1_all <- c(S1S1_all, S1S1) 
      Mean_DeltaT_all <- c(Mean_DeltaT_all, Mean_DeltaT)
      Mean_DeltaS_all <- c(Mean_DeltaS_all, Mean_DeltaS)
        
      ########################  end pos def 
    }
  }
  Results.ICA <- data.frame(Results.ICA)
  rownames(Results.ICA) <- NULL
  Total.Num.Matrices <- dim(Results.ICA)[1]
  
  # ISTE results
  Predicted_Delta_T_all <- data.frame(Predicted_Delta_T_all); names(Predicted_Delta_T_all) <- Delta_S
  PI_Interval_Low_all <- data.frame(PI_Interval_Low_all); names(PI_Interval_Low_all) <- Delta_S
  PI_Interval_Up_all <- data.frame(PI_Interval_Up_all); names(PI_Interval_Up_all) <- Delta_S
  
  
  fit <- 
    list(#ISTE_Low_PI_0=STE_Low_PI_0_all, ISTE_Up_PI_0=STE_Up_PI_0_all,
         ISTE_Low_PI=STE_Low_PI_all, ISTE_Up_PI=STE_Up_PI_all, MSE=MSE_all, gamma0=b0_all, gamma1=b1_all, 
         Delta_S_For_Which_Delta_T_equal_0=Exp_Delta_S_Delta_T_equal_O_all, S_squared_pred=S_squared_pred_all, 
         Predicted_Delta_T=Predicted_Delta_T_all, PI_Interval_Low=PI_Interval_Low_all, PI_Interval_Up=PI_Interval_Up_all,
         T0T0=T0T0_all, T1T1=T1T1_all, S0S0=S0S0_all, S1S1=S1S1_all, 
         Mean_DeltaT=Mean_DeltaT_all, Mean_DeltaS=Mean_DeltaS_all,
         Total.Num.Matrices=Total.Num.Matrices, Pos.Def=Results.ICA[,1:6], ICA=Results.ICA$ICA,
         zeta.PI=0.05, PI.Bound=PI.Bound, PI.Lower=PI.Lower, #GoodSurr=Results.ICA[,7:9], 
         Delta_S=Delta_S,
         Call=match.call())
  
  class(fit) <- "ISTE.ContCont"
  fit
}



# PLOT FUNCTION

plot.ISTE.ContCont <- function(x, Outcome="ISTE", breaks=50, ...){
     Object <- x

     if(Object$PI.Bound != 0){
       cat("\nNote that PI.Bound =", Object$PI.Bound, "was requested in the computation of ISTE.\n\n")
     }
     
     if (Outcome=="ISTE"){
       if (x$PI.Bound==0 & x$PI.Lower==TRUE){
         hist(x$ISTE_Low_PI, xlab="ISTE", main=" ", col="grey", breaks=breaks, ...)
          }
       if (x$PI.Bound==0 & x$PI.Lower==FALSE){
         hist(x$ISTE_Up_PI, xlab="ISTE", main=" ", col="grey",  breaks=breaks, ...)
          }
       if (x$PI.Bound!=0 & x$PI.Lower==TRUE){
         hist(x$ISTE_Low_PI, xlab="ISTE", main=" ", col="grey", breaks=breaks,  ...)
          }
       if (x$PI.Bound!=0 & x$PI.Lower==FALSE){
         hist(x$ISTE_Up_PI, xlab="ISTE", main=" ", col="grey", breaks=breaks,  ...)
          }
       }

     if (Outcome=="MSE"){
       hist(x$MSE, xlab="MSE", main=" ", col="grey",  breaks=breaks, ...)
     }

     if (Outcome=="gamma0"){
       hist(x$gamma0, xlab=expression(gamma[0]), main=" ", col="grey",  breaks=breaks, ...)
     }

     if (Outcome=="gamma1"){
       hist(x$gamma1, xlab=expression(gamma[1]), main=" ", col="grey",  breaks=breaks, ...)
     }

     if (Outcome=="Exp.DeltaT"){
       
       for (i in 1: length(x$Delta_S)){
         Delta_S_hier <- as.numeric(x$Delta_S[i])
         
         hist(x$Predicted_Delta_T[,i], col="grey", xlab=" ", ...,  breaks=breaks, 
              main=eval(substitute(paste("E(Delta T|Delta S = ", Delta_S_hier, ")", sep=""))))
       }
     }
       
      if (Outcome=="Exp.DeltaT.Low.PI"){
         
         for (i in 1: length(x$Delta_S)){
           Delta_S_hier <- as.numeric(x$Delta_S[i])
          hist(x$PI_Interval_Low[,i], col="grey", xlab=" ", breaks=breaks, ...,  
                main=eval(substitute(paste("Lower PI E(Delta T|Delta S = ", Delta_S_hier, ")", sep=""))))
         }
      }
         
        if (Outcome=="Exp.DeltaT.Up.PI"){
          
          for (i in 1: length(x$Delta_S)){
            Delta_S_hier <- as.numeric(x$Delta_S[i])
            hist(x$PI_Interval_Up[,i], col="grey", xlab=" ", breaks=breaks, ...,  
                 main=eval(substitute(paste("Upper PI E(Delta T|Delta S = ", Delta_S_hier, ")", sep=""))))
          }
          
        }
     
     if (Outcome=="Delta_S_For_Which_Delta_T_equal_0"){
       hist(x$Delta_S_For_Which_Delta_T_equal_0, 
      xlab=expression(paste("E(", Delta, "T|", Delta, "S>", omega, ")>0")),
      main=" ", col="grey",  breaks=breaks, ...)
     }
     
 }



summary.ISTE.ContCont <- function(object, ..., Object){
  if (missing(Object)){Object <- object}

  if (Object$PI.Bound==0 & Object$PI.Lower==TRUE){
    Out_here <- Object$ISTE_Low_PI
  }
  if (Object$PI.Bound==0 & Object$PI.Lower==FALSE){
    Out_here <- Object$ISTE_Up_PI
  }
  if (Object$PI.Bound!=0 & Object$PI.Lower==TRUE){
    Out_here <- Object$ISTE_Low_PI
  }
  if (Object$PI.Bound!=0 & Object$PI.Lower==FALSE){
    Out_here <- Object$ISTE_Up_PI
  }
  
  Out_here <- na.exclude(Out_here)
  
  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }

  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Total number of positive definite Sigma matrices and ISTE values")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat(nrow(Object$Pos.Def))
  cat("\n\n\n# ISTE results")
  cat("\n#~~~~~~~~~~~~~\n\n")
  cat("Mean (SD) ICA: ", format(round(mean(Out_here), 4), nsmall = 4), " (", format(round(sd(Out_here), 4), nsmall = 4), ")",
      "  [min: ", format(round(min(Out_here), 4), nsmall = 4), "; max: ",  format(round(max(Out_here), 4), nsmall = 4), "]", sep="")
  cat("\nMode ICA: ", format(round(mode(Out_here)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the ICA distribution: \n\n")
  quant <- quantile(Out_here, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
  
  if(Object$PI.Bound != 0){
    cat("\nNote: PI.Bound =", Object$PI.Bound, "was requested in the computation of ISTE.")
  }
  
}

