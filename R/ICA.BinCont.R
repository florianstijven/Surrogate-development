ICA.BinCont <- function(Dataset, Surr, True, Treat, 
                        G_pi_10 = c(0,1),
                        G_rho_01_00=c(-1,1), 
                        G_rho_01_01=c(-1,1), 
                        G_rho_01_10=c(-1,1), 
                        G_rho_01_11=c(-1,1), 
                        Theta.S_0, 
                        Theta.S_1,
                        M=1000, Seed=123, 
                        Monotonicity=FALSE,
                        Independence=FALSE,
                        Plots=TRUE, Save.Plots="No", 
                        Show.Details=FALSE){          

  Nboot <- 1000
  
  set.seed(Seed)
  
  totaal <- M*length(G_pi_10)*length(G_rho_01_00)*length(G_rho_01_01)*length(G_rho_01_10)*length(G_rho_01_11)
  
  Diff.Sigma = FALSE

  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  
  if (min(na.exclude(Treat))!=c(-1)) {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (max(na.exclude(Treat))!=c(1))  {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (length(unique(na.exclude(True)))>2) {stop("\nThe true endpoint should be binary.")}
  if (min(na.exclude(True))!=c(0)) {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  if (max(na.exclude(True))!=c(1))  {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  
  data_no_miss <- data.frame(na.exclude(cbind(Surr, True, Treat)))
  data_conttreat <- subset(data_no_miss, Treat=="-1")      
  data_exptreat <- subset(data_no_miss, Treat=="1")      
  S_0 <- data_conttreat$Surr
  S_1 <- data_exptreat$Surr
  
  Surr <- na.exclude(Surr)   #LS
  True <- na.exclude(True)   #LS
  Treat <- na.exclude(Treat) #LS

  pi_punt1 <- mean(data_exptreat$True) # E(T|Z=1)
  pi_punt0 <- 1 - pi_punt1 # 1-E(T|Z=1)
  pi_1punt <- mean(data_conttreat$True) # E(T|Z=0)
  pi_0punt <- 1 - pi_1punt # 1 - E(T|Z=0)

  if(Monotonicity==TRUE){
    if (pi_punt0 > pi_0punt | pi_1punt > pi_punt1) {stop("\nThe monotonicity assumption is not compatible with the data")}
  }  

  count <- 0
  R2_H_all <- R2_C_all <- PD_OK_all <- pi_00_all <- pi_01_all <- pi_10_all <- pi_11_all <- NULL
  G_rho_01_00_all <- G_rho_01_01_all <- G_rho_01_10_all <- G_rho_01_11_all <- NULL
  pi_Delta_T_min1_all <- pi_Delta_T_0_all <- pi_Delta_T_1_all <- NULL
  pi_0_00_e_all <- pi_0_01_e_all <- pi_0_10_e_all <- pi_0_11_e_all <- NULL 
  mu_0_00_all <- mu_0_01_all <- mu_0_10_all <- mu_0_11_all <- NULL
  sigma_00_00_all <- sigma_00_01_all <- sigma_00_10_all <- sigma_00_11_all <- NULL
  pi_1_00_e_all <- pi_1_01_e_all <- pi_1_10_e_all <- pi_1_11_e_all <- NULL 
  mu_1_00_all <- mu_1_01_all <- mu_1_10_all <- mu_1_11_all <- NULL
  sigma_11_00_all <- sigma_11_01_all <- sigma_11_10_all <- sigma_11_11_all <- NULL
  mean_Y_S0_all <- mean_Y_S1_all <- NULL
  var_Y_S0_all <- var_Y_S1_all <- NULL
  dev_S0_all <- dev_S1_all <- NULL
  aantal <- 0
  voor_seed <- Seed

  while (aantal < M){
    
    count <- count+1
    voor_seed <- voor_seed + 1

    G_rho_01_00_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
    G_rho_01_01_hier <- runif(1, min=G_rho_01_01[1], max=G_rho_01_01[2])
    G_rho_01_10_hier <- runif(1, min=G_rho_01_10[1], max=G_rho_01_10[2])
    G_rho_01_11_hier <- runif(1, min=G_rho_01_11[1], max=G_rho_01_11[2])

    if(Show.Details==TRUE){
      flush.console()
      cat("\n\nG_rho:", cbind(G_rho_01_00_hier, G_rho_01_01_hier, G_rho_01_10_hier, G_rho_01_11_hier), "\n")
    }

    set.seed(count)
    if(Monotonicity==TRUE){
      G_pi_10 = c(0,0)
    }
    
    if(Independence==FALSE){
      pi_10_hier <- runif(1, min=G_pi_10[1], max=G_pi_10[2])   
      pi_00_hier <- pi_punt0 - pi_10_hier #UD
      pi_01_hier <- pi_0punt -  pi_00_hier #UD
      pi_11_hier <- pi_1punt -  pi_10_hier #UD
    }
    
    if(Independence==TRUE){
      pi_10_hier <- pi_1punt*pi_punt0   
      pi_00_hier <- pi_0punt*pi_punt0
      pi_01_hier <- pi_0punt*pi_punt1
      pi_11_hier <- pi_1punt*pi_punt1
    }
    
    
    if ((pi_00_hier >= 0 & pi_01_hier >= 0 & pi_10_hier >= 0 & pi_11_hier >= 0) & 
        (pi_00_hier <= 1 & pi_01_hier <= 1 & pi_10_hier <= 1 & pi_11_hier <= 1)){
      
      if(Show.Details==TRUE){
        flush.console()
        cat("\n\nPi's used here:", cbind(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), "\n")
      }
      
      pi_Delta_T_min1 <- pi_10_hier
      pi_Delta_T_0 <- pi_00_hier + pi_11_hier
      pi_Delta_T_1 <- pi_01_hier 

      mix1 <- mix2 <- mix.object <- Sum_M_S_0 <- Sum_M_S_1 <- M1 <- M2 <- M3 <- M4 <- NULL

      # Function to minimize deviance of normal mixtures
      ## Unequal sigma's
      if (Diff.Sigma==TRUE){
        mixt.deviance <- function(theta,data){
          pi_00 <- pi_00_hier
          pi_01 <- pi_01_hier
          pi_10 <- pi_10_hier
          pi_11 <- pi_11_hier
  
          mu_1 <- theta[1]
          mu_2 <- theta[2]
          mu_3 <- theta[3]
          mu_4 <- theta[4]
  
          sigma_1 <- theta[5]
          sigma_2 <- theta[6]
          sigma_3 <- theta[7]
          sigma_4 <- theta[8]
  
          pdf <- pi_00*dnorm(data,mu_1,sigma_1) + 
                 pi_01*dnorm(data,mu_2,sigma_2) +
                 pi_10*dnorm(data,mu_3,sigma_3) + 
                 pi_11*dnorm(data,mu_4,sigma_4) 
  
          deviance <- -2*sum(log(pdf))
          return(deviance)
        }
      }

      ## Equal sigma's
      if (Diff.Sigma==FALSE){
        mixt.deviance <- function(theta,data){
          pi_00 <- pi_00_hier
          pi_01 <- pi_01_hier
          pi_10 <- pi_10_hier
          pi_11 <- pi_11_hier
    
          mu_1 <- theta[1]
          mu_2 <- theta[2]
          mu_3 <- theta[3]
          mu_4 <- theta[4]
    
          sigma_all <- theta[5]

          pdf <- pi_00*dnorm(data,mu_1,sigma_all) + 
                 pi_01*dnorm(data,mu_2,sigma_all) +
                 pi_10*dnorm(data,mu_3,sigma_all) + 
                 pi_11*dnorm(data,mu_4,sigma_all) 
    
          deviance <- -2*sum(log(pdf))
          return(deviance)
        }
      }

      # nonlinear minimization for the mixture deviance function using Newton-Raphson type algorithm
      options(warn = -999)
      mix1 <- nlm(f = mixt.deviance, p = Theta.S_0, na.exclude(S_0), iterlim = 100000)
      mix2 <- nlm(f = mixt.deviance, p = Theta.S_1, na.exclude(S_1), iterlim = 100000)
      options(warn = 1)

      if (exists("mix1")==TRUE & exists("mix2")==TRUE & is.null(mix1)==FALSE & is.null(mix2)==FALSE &
          mix1$code<=2 & mix2$code<=2){
        
        # S0
        mix.object <- mix1
        k <- 4
        x_S_0 <- x <- sort(na.exclude(S_0))
        a <- hist(x, plot = FALSE) 

        if (Diff.Sigma==FALSE){
          M1 <- pi_00_hier * dnorm(x, mean = mix.object$estimate[1], sd = mix.object$estimate[5])
          M2 <- pi_01_hier * dnorm(x, mean = mix.object$estimate[2], sd = mix.object$estimate[5])
          M3 <- pi_10_hier * dnorm(x, mean = mix.object$estimate[3], sd = mix.object$estimate[5]) 
          M4 <- pi_11_hier * dnorm(x, mean = mix.object$estimate[4], sd = mix.object$estimate[5])
        }

        if (Diff.Sigma==TRUE){  
          M1 <- pi_00_hier * dnorm(x, mean = mix.object$estimate[1], sd = mix.object$estimate[5])
          M2 <- pi_01_hier * dnorm(x, mean = mix.object$estimate[2], sd = mix.object$estimate[6])
          M3 <- pi_10_hier * dnorm(x, mean = mix.object$estimate[3], sd = mix.object$estimate[7]) 
          M4 <- pi_11_hier * dnorm(x, mean = mix.object$estimate[4], sd = mix.object$estimate[8])
        }  
        
        Sum_M_S_0 <- M1 + M2 + M3 + M4
        S0_M1 <- M1; S0_M2 <- M2; S0_M3 <- M3; S0_M4 <- M4 
        dev_S0 <- -2*sum(log(Sum_M_S_0))

        # S1
        mix.object <- mix2
        k <- 4
        x_S_1 <- x <- sort(na.exclude(S_1))
        a <- hist(x, plot = FALSE) 

        if (Diff.Sigma==FALSE){
          M1 <- pi_00_hier * dnorm(x, mean = mix.object$estimate[1], sd = mix.object$estimate[5])
          M2 <- pi_01_hier * dnorm(x, mean = mix.object$estimate[2], sd = mix.object$estimate[5])
          M3 <- pi_10_hier * dnorm(x, mean = mix.object$estimate[3], sd = mix.object$estimate[5]) 
          M4 <- pi_11_hier * dnorm(x, mean = mix.object$estimate[4], sd = mix.object$estimate[5])
        }

        if (Diff.Sigma==TRUE){  
          M1 <- pi_00_hier * dnorm(x, mean = mix.object$estimate[1], sd = mix.object$estimate[5])
          M2 <- pi_01_hier * dnorm(x, mean = mix.object$estimate[2], sd = mix.object$estimate[6])
          M3 <- pi_10_hier * dnorm(x, mean = mix.object$estimate[3], sd = mix.object$estimate[7]) 
          M4 <- pi_11_hier * dnorm(x, mean = mix.object$estimate[4], sd = mix.object$estimate[8])
        }  

        Sum_M_S_1 <- M1 + M2 + M3 + M4
        S1_M1 <- M1; S1_M2 <- M2; S1_M3 <- M3; S1_M4 <- M4 
        dev_S1 <- -2*sum(log(Sum_M_S_1))

        if (Plots == TRUE){
          if (Save.Plots!="No"){
            options(warn=-1)
            dir.create(Save.Plots)
            options(warn=0)
  
            pdf(paste(Save.Plots, "/Run", (aantal+1), ".pdf", sep=""),width=6,height=4,paper='special') 
            par(mfrow=c(1,2))
            max_val <- max(max(density(S_0)$y), max(Sum_M_S_0))
            hist(S_0, freq = F, ylim=c(0, max_val), xlab=expression(S[0]), 
            main=bquote(paste(S[0], " (run "~.(aantal+1), ")")), col="grey")
            lines(sort(x_S_0), Sum_M_S_0, col="red", lwd=4, lty=1)
            lines(sort(x_S_0), S0_M1, col="blue", lwd=2, lty=3)
            lines(sort(x_S_0), S0_M2, col="blue", lwd=2, lty=3)
            lines(sort(x_S_0), S0_M3, col="blue", lwd=2, lty=3)
            lines(sort(x_S_0), S0_M4, col="blue", lwd=2, lty=3)
  
            max_val <- max(max(density(S_1)$y), max(Sum_M_S_1))
            hist(S_1, freq = F, ylim=c(0, max_val), xlab=expression(S[1]), 
            main=bquote(paste(S[1], " (run "~.(aantal+1), ")")), col="grey")
            lines(sort(x_S_1), Sum_M_S_1, col="red", lwd=4, lty=1)
            lines(sort(x_S_1), S1_M1, col="blue", lwd=2, lty=3)
            lines(sort(x_S_1), S1_M2, col="blue", lwd=2, lty=3)
            lines(sort(x_S_1), S1_M3, col="blue", lwd=2, lty=3)
            lines(sort(x_S_1), S1_M4, col="blue", lwd=2, lty=3)
            dev.off()
          }
 
          # Show plots on screen   
          par(mfrow=c(1,2))
          max_val <- max(max(density(S_0)$y), max(Sum_M_S_0))
          hist(S_0, freq = F, ylim=c(0, max_val), xlab=expression(S[0]), 
          main=bquote(paste(S[0], " (run "~.(aantal+1), ")")), col="grey")
          lines(sort(x_S_0), Sum_M_S_0, col="red", lwd=4, lty=1)
          lines(sort(x_S_0), S0_M1, col="blue", lwd=2, lty=3)
          lines(sort(x_S_0), S0_M2, col="blue", lwd=2, lty=3)
          lines(sort(x_S_0), S0_M3, col="blue", lwd=2, lty=3)
          lines(sort(x_S_0), S0_M4, col="blue", lwd=2, lty=3)
  
          max_val <- max(max(density(S_1)$y), max(Sum_M_S_1))
          hist(S_1, freq = F, ylim=c(0, max_val), xlab=expression(S[1]), 
          main=bquote(paste(S[1], " (run "~.(aantal+1), ")")), col="grey")
          lines(sort(x_S_1), Sum_M_S_1, col="red", lwd=4, lty=1)
          lines(sort(x_S_1), S1_M1, col="blue", lwd=2, lty=3)
          lines(sort(x_S_1), S1_M2, col="blue", lwd=2, lty=3)
          lines(sort(x_S_1), S1_M3, col="blue", lwd=2, lty=3)
          lines(sort(x_S_1), S1_M4, col="blue", lwd=2, lty=3)
          par(mfrow=c(1,1))
        }

        # mixture components f(S_0)
        pi_0_00_e <- pi_00_hier
        pi_0_01_e <- pi_01_hier  
        pi_0_10_e <- pi_10_hier
        pi_0_11_e <- pi_11_hier
        mu_0_00 <- mix1$estimate[1] 
        mu_0_01 <- mix1$estimate[2] 
        mu_0_10 <- mix1$estimate[3]
        mu_0_11 <- mix1$estimate[4]
        if (Diff.Sigma==TRUE){
          sigma_00_00 <- mix1$estimate[5]**2  
          sigma_00_01 <- mix1$estimate[6]**2 
          sigma_00_10 <- mix1$estimate[7]**2 
          sigma_00_11 <- mix1$estimate[8]**2} 
        if (Diff.Sigma==FALSE){sigma_00_00 <- sigma_00_01 <- sigma_00_10 <- sigma_00_11 <- mix1$estimate[5]**2} 

        # mixture components f(S_1)
        pi_1_00_e <- pi_00_hier 
        pi_1_01_e <- pi_01_hier
        pi_1_10_e <- pi_10_hier
        pi_1_11_e <- pi_11_hier
        mu_1_00 <- mix2$estimate[1]
        mu_1_01 <- mix2$estimate[2]
        mu_1_10 <- mix2$estimate[3]
        mu_1_11 <- mix2$estimate[4] 
        if (Diff.Sigma==TRUE){
          sigma_11_00 <- mix2$estimate[5]**2
          sigma_11_01 <- mix2$estimate[6]**2
          sigma_11_10 <- mix2$estimate[7]**2 
          sigma_11_11 <- mix2$estimate[8]**2}
        if (Diff.Sigma==FALSE){sigma_11_00 <- sigma_11_01 <- sigma_11_10 <- sigma_11_11 <- mix2$estimate[5]**2}

        # Check PD covariance matrices 
        mat_a <- matrix(c(sigma_00_00, 
                         (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))), 
                         (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))), 
                          sigma_11_00), nrow = 2, byrow = TRUE)
        eigen_mat_a <- min(eigen((mat_a))$values)

        mat_b <- matrix(c(sigma_00_01, 
                         (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))), 
                         (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))), 
                          sigma_11_01), nrow = 2, byrow = TRUE)
        eigen_mat_b <- min(eigen((mat_b))$values)

        mat_c <- matrix(c(sigma_00_10, 
                         (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))), 
                         (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))), 
                          sigma_11_10), nrow = 2, byrow = TRUE)
        eigen_mat_c <- min(eigen((mat_c))$values)

        mat_d <- matrix(c(sigma_00_11, 
                         (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))), 
                         (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))), 
                          sigma_11_11), nrow = 2, byrow = TRUE)
        eigen_mat_d <- min(eigen((mat_d))$values)

        # check CN 
        singular_mat_a <- svd(mat_a)$d
        CN_mat_a <- max(singular_mat_a)/min(singular_mat_a)
        singular_mat_b <- svd(mat_b)$d
        CN_mat_b <- max(singular_mat_b)/min(singular_mat_b)
        singular_mat_c <- svd(mat_c)$d
        CN_mat_c <- max(singular_mat_c)/min(singular_mat_c)
        singular_mat_d <- svd(mat_d)$d
        CN_mat_d <- max(singular_mat_d)/min(singular_mat_d)

        if (eigen_mat_a > 0 & eigen_mat_b > 0 & eigen_mat_c > 0 & eigen_mat_d > 0 & 
            CN_mat_a < 50 & CN_mat_b < 50 & CN_mat_c < 50 & CN_mat_d < 50) {
          
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nMin Eigenvalues:", cbind(eigen_mat_a, eigen_mat_b, eigen_mat_c, eigen_mat_d), "\n")
          }
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nCondition number:", cbind(CN_mat_a, CN_mat_b, CN_mat_c, CN_mat_d), "\n")
          }
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nChecks passed...")
          }
          
          # Mixture check
          mean_Y_S0 <- (pi_0_00_e*mu_0_00) + (pi_0_01_e*mu_0_01) + (pi_0_10_e*mu_0_10) + (pi_0_11_e*mu_0_11)
          mean_Y_S1 <- (pi_1_00_e*mu_1_00) + (pi_1_01_e*mu_1_01) + (pi_1_10_e*mu_1_10) + (pi_1_11_e*mu_1_11)
          var_Y_S0 <- ((pi_0_00_e*mu_0_00^2)+(pi_0_01_e*mu_0_01^2)+(pi_0_10_e*mu_0_10^2)+(pi_0_11_e*mu_0_11^2)) - mean_Y_S0^2 +
                      ((pi_0_00_e*sigma_00_00)+(pi_0_01_e*sigma_00_01)+(pi_0_10_e*sigma_00_10)+(pi_0_11_e*sigma_00_11))
          var_Y_S1 <- ((pi_1_00_e*mu_1_00^2)+(pi_1_01_e*mu_1_01^2)+(pi_1_10_e*mu_1_10^2)+(pi_1_11_e*mu_1_11^2)) - mean_Y_S1^2 +
                      ((pi_1_00_e*sigma_11_00)+(pi_1_01_e*sigma_11_01)+(pi_1_10_e*sigma_11_10)+(pi_1_11_e*sigma_11_11))

          # Calculation of ICA
          mu_s00 <- mu_1_00 - mu_0_00
          mu_s01 <- mu_1_01 - mu_0_01
          mu_s10 <- mu_1_10 - mu_0_10
          mu_s11 <- mu_1_11 - mu_0_11
          omega_1 <- pi_00_hier / (pi_00_hier + pi_11_hier); omega_2 <- 1 - omega_1

          sigma_s00 <- sigma_00_00 + sigma_11_00 - (2 * (sqrt(sigma_00_00 * sigma_11_00) * G_rho_01_00_hier))
          sigma_s01 <- sigma_00_01 + sigma_11_01 - (2 * (sqrt(sigma_00_01 * sigma_11_01) * G_rho_01_01_hier))
          sigma_s10 <- sigma_00_10 + sigma_11_10 - (2 * (sqrt(sigma_00_10 * sigma_11_10) * G_rho_01_10_hier))
          sigma_s11 <- sigma_00_11 + sigma_11_11 - (2 * (sqrt(sigma_00_11 * sigma_11_11) * G_rho_01_11_hier))

          f_Delta_S <- function(val){
            v <- NA
            v <- (pi_00_hier * dnorm(val, mu_s00, sd = sqrt(sigma_s00))) + 
                 (pi_01_hier * dnorm(val, mu_s01, sd = sqrt(sigma_s01))) +
                 (pi_10_hier * dnorm(val, mu_s10, sd = sqrt(sigma_s10))) +
                 (pi_11_hier * dnorm(val, mu_s11, sd = sqrt(sigma_s11)))
            return(v)
          }

          # I_10
          set.seed(1)
          S1 <- rnorm(n=10000, mean = mu_s10, sd = sqrt(sigma_s10))
          I_10 <- mean(log(dnorm(S1, mu_s10, sd = sqrt(sigma_s10)) / f_Delta_S(S1)))

          # I_01
          set.seed(2)
          S2 <- rnorm(n=10000, mean = mu_s01, sd = sqrt(sigma_s01))
          I_01 <- mean(log(dnorm(S2, mu_s01, sd = sqrt(sigma_s01)) / f_Delta_S(S2)))

          # I_00
          set.seed(3)
          S3 <- rnorm(n=10000, mean = mu_s00, sd = sqrt(sigma_s00))
          I_00 <- mean(log(((omega_1 * dnorm(S3, mu_s00, sd = sqrt(sigma_s00))) +
                          ((omega_2 * dnorm(S3, mu_s11, sd = sqrt(sigma_s11))))) / 
                          f_Delta_S(S3)))

          # I_11
          set.seed(4)
          S4 <- rnorm(n=10000, mean = mu_s11, sd = sqrt(sigma_s11))
          I_11 <-mean(log(((omega_1 * dnorm(S4, mu_s00, sd = sqrt(sigma_s00))) + 
                         ((omega_2 * dnorm(S4, mu_s11, sd = sqrt(sigma_s11))))) / 
                         f_Delta_S(S4)))

          I_Delta_T_Delta_S <- (pi_00_hier * I_00) + (pi_01_hier * I_01) + (pi_10_hier * I_10) + (pi_11_hier * I_11)
          H_Delta_T <- -((pi_Delta_T_min1 * log(pi_Delta_T_min1)) + 
                         (pi_Delta_T_0 * log(pi_Delta_T_0)) +
                         (pi_Delta_T_1 * log(pi_Delta_T_1)))
          if(Monotonicity==TRUE){
            H_Delta_T <- -((pi_Delta_T_0 * log(pi_Delta_T_0)) +
                           (pi_Delta_T_1 * log(pi_Delta_T_1)))
          }
          R2_H <- I_Delta_T_Delta_S / H_Delta_T
          R2_C <- 1 - exp(-2*I_Delta_T_Delta_S)

          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nR2H:", R2_H, "\n")
          }
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nR2C:", R2_C, "\n")
          }

          # aantal
          aantal <- aantal + 1
          cat("\n", (aantal/M)*100, "% done. \n", sep="") 

          # save results for output
          R2_H_all <- c(R2_H_all, R2_H)
          R2_C_all <- c(R2_C_all, R2_C)
          G_rho_01_00_all <- cbind(G_rho_01_00_all, G_rho_01_00_hier)
          G_rho_01_01_all <- cbind(G_rho_01_01_all, G_rho_01_01_hier)
          G_rho_01_10_all <- cbind(G_rho_01_10_all, G_rho_01_10_hier)
          G_rho_01_11_all <- cbind(G_rho_01_11_all, G_rho_01_11_hier)
          pi_00_all <- cbind(pi_00_all, pi_00_hier) 
          pi_01_all <- cbind(pi_01_all, pi_01_hier)
          pi_10_all <- cbind(pi_10_all, pi_10_hier)
          pi_11_all <- cbind(pi_11_all, pi_11_hier) 
          pi_Delta_T_min1_all <- cbind(pi_Delta_T_min1_all, pi_Delta_T_min1)
          pi_Delta_T_0_all <- cbind(pi_Delta_T_0_all, pi_Delta_T_0)
          pi_Delta_T_1_all <- cbind(pi_Delta_T_1_all, pi_Delta_T_1)
          pi_0_00_e_all <- cbind(pi_0_00_e_all, pi_0_00_e)
          pi_0_01_e_all <- cbind(pi_0_01_e_all, pi_0_01_e)
          pi_0_10_e_all <- cbind(pi_0_10_e_all, pi_0_10_e)
          pi_0_11_e_all <- cbind(pi_0_11_e_all, pi_0_11_e)
          mu_0_00_all <- cbind(mu_0_00_all, mu_0_00)
          mu_0_01_all <- cbind(mu_0_01_all, mu_0_01)
          mu_0_10_all <- cbind(mu_0_10_all, mu_0_10)
          mu_0_11_all <- cbind(mu_0_11_all, mu_0_11)
          sigma_00_00_all <- cbind(sigma_00_00_all, sigma_00_00)
          sigma_00_01_all <- cbind(sigma_00_01_all, sigma_00_01)
          sigma_00_10_all <- cbind(sigma_00_10_all, sigma_00_10)
          sigma_00_11_all <- cbind(sigma_00_11_all, sigma_00_11)
          pi_1_00_e_all <- cbind(pi_1_00_e_all, pi_1_00_e)
          pi_1_01_e_all <- cbind(pi_1_01_e_all, pi_1_01_e)
          pi_1_10_e_all <- cbind(pi_1_10_e_all, pi_1_10_e)
          pi_1_11_e_all <- cbind(pi_1_11_e_all, pi_1_11_e)
          mu_1_00_all <- cbind(mu_1_00_all, mu_1_00)
          mu_1_01_all <- cbind(mu_1_01_all, mu_1_01)
          mu_1_10_all <- cbind(mu_1_10_all, mu_1_10)
          mu_1_11_all <- cbind(mu_1_11_all, mu_1_11)
          sigma_11_00_all <- cbind(sigma_11_00_all, sigma_11_00)
          sigma_11_01_all <- cbind(sigma_11_01_all, sigma_11_01)
          sigma_11_10_all <- cbind(sigma_11_10_all, sigma_11_10)
          sigma_11_11_all <- cbind(sigma_11_11_all, sigma_11_11)
          mean_Y_S0_all <- cbind(mean_Y_S0_all, mean_Y_S0)
          mean_Y_S1_all <- cbind(mean_Y_S1_all, mean_Y_S1)
          var_Y_S0_all <- cbind(var_Y_S0_all, var_Y_S0)
          var_Y_S1_all <- cbind(var_Y_S1_all, var_Y_S1)
          dev_S0_all <- cbind(dev_S0_all, dev_S0)
          dev_S1_all <- cbind(dev_S1_all, dev_S1)

        } # min eigen and CN requirements 
      } # mixture distribution available
    } # pi's requirements (>=0 and <=1)
  } # end the while loop
       
  fit <- list(R2_H=R2_H_all, R2_C=R2_C_all, 
              pi_00=as.numeric(pi_00_all), pi_01=as.numeric(pi_01_all), 
              pi_10=as.numeric(pi_10_all), pi_11=as.numeric(pi_11_all), 
              G_rho_01_00=as.numeric(G_rho_01_00_all), 
              G_rho_01_01=as.numeric(G_rho_01_01_all),
              G_rho_01_10=as.numeric(G_rho_01_10_all), 
              G_rho_01_11=as.numeric(G_rho_01_11_all),
              pi_Delta_T_min1=as.numeric(pi_Delta_T_min1_all), 
              pi_Delta_T_0=as.numeric(pi_Delta_T_0_all),
              pi_Delta_T_1=as.numeric(pi_Delta_T_1_all),
              pi_0_00=as.numeric(pi_0_00_e_all),
              pi_0_01=as.numeric(pi_0_01_e_all),
              pi_0_10=as.numeric(pi_0_10_e_all),
              pi_0_11=as.numeric(pi_0_11_e_all), 
              mu_0_00=as.numeric(mu_0_00_all),
              mu_0_01=as.numeric(mu_0_01_all),
              mu_0_10=as.numeric(mu_0_10_all),
              mu_0_11=as.numeric(mu_0_11_all),
              sigma2_00_00=as.numeric(sigma_00_00_all),
              sigma2_00_01=as.numeric(sigma_00_01_all),
              sigma2_00_10=as.numeric(sigma_00_10_all),
              sigma2_00_11=as.numeric(sigma_00_11_all),
              pi_1_00=as.numeric(pi_1_00_e_all),
              pi_1_01=as.numeric(pi_1_01_e_all),
              pi_1_10=as.numeric(pi_1_10_e_all),
              pi_1_11=as.numeric(pi_1_11_e_all), 
              mu_1_00=as.numeric(mu_1_00_all),
              mu_1_01=as.numeric(mu_1_01_all),
              mu_1_10=as.numeric(mu_1_10_all),
              mu_1_11=as.numeric(mu_1_11_all),
              sigma2_11_00=as.numeric(sigma_11_00_all),
              sigma2_11_01=as.numeric(sigma_11_01_all),
              sigma2_11_10=as.numeric(sigma_11_10_all),
              sigma2_11_11=as.numeric(sigma_11_11_all),
              mean_Y_S0=as.numeric(mean_Y_S0_all),
              mean_Y_S1=as.numeric(mean_Y_S1_all),
              var_Y_S0=as.numeric(var_Y_S0_all),
              var_Y_S1=as.numeric(var_Y_S1_all),
              dev_S0=as.numeric(dev_S0_all),
              dev_S1=as.numeric(dev_S1_all),
              mean.S0=mean(S_0, na.rm=TRUE),
              var.S0=var(S_0, na.rm=TRUE),
              mean.S1=mean(S_1, na.rm=TRUE),
              var.S1=var(S_1, na.rm=TRUE),
              Call=match.call())   

  class(fit) <- "ICA.BinCont"
  fit

} # end of function  

