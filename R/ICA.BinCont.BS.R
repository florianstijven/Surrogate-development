ICA.BinCont.BS <- function(Dataset, Surr, True, Treat,
                           BS = TRUE,
                           nb = 300,
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
                           HAA=FALSE,
                           Cond_ind=FALSE,
                           Plots=TRUE, Save.Plots="No", 
                           Show.Details=FALSE){
  
  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  
  data_no_miss.original <- data.frame(na.exclude(cbind(Surr, True, Treat)))
  if (min(na.exclude(Treat))!=c(-1)) {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (max(na.exclude(Treat))!=c(1))  {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (length(unique(na.exclude(True)))>2) {stop("\nThe true endpoint should be binary.")}
  if (min(na.exclude(True))!=c(0)) {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  if (max(na.exclude(True))!=c(1))  {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  
  no_boots_all <- NULL
  R2_H_all <- NULL
  pi_00_all <- pi_01_all <- pi_10_all <- pi_11_all <- NULL
  G_rho_01_00_all <- G_rho_01_01_all <- G_rho_01_10_all <- G_rho_01_11_all <- NULL
  mu_0_00_all <- mu_0_01_all <- mu_0_10_all <- mu_0_11_all <- NULL
  mu_1_00_all <- mu_1_01_all <- mu_1_10_all <- mu_1_11_all <- NULL
  sigma_00_all <- sigma_11_all <- NULL
  voor_seed <- Seed
  
  for (i in 1:nb) {
    data_no_miss <- data_no_miss.original
    
    voor_seed <- voor_seed + 1
    set.seed(voor_seed)
    if (BS==TRUE) {
      #Boot_Sample <- data_no_miss[sample(1:dim(data_no_miss)[1], replace = TRUE),]
      A1 <- subset(data_no_miss,data_no_miss$Treat==-1)  # keep the number of treatment group the same
      B1 <- subset(data_no_miss,data_no_miss$Treat==1)
      A2 <- A1[sample(nrow(A1),replace=T),]
      B2 <- B1[sample(nrow(B1),replace=T),]
      data_no_miss <- as.data.frame(rbind(A2,B2))                        
    }
    
    if (BS==FALSE) {
      data_no_miss <- as.data.frame(data_no_miss)                        
    }
    
    Dataset <- data_no_miss
    Surr <- data_no_miss$Surr
    True <- data_no_miss$True
    Treat <- data_no_miss$Treat
    no_boots <- i
    
    fit.ica.bincont <- ICA.BinCont(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat,
                                   G_pi_10 = G_pi_10,
                                   G_rho_01_00 = G_rho_01_00, 
                                   G_rho_01_01 = G_rho_01_01, 
                                   G_rho_01_10 = G_rho_01_10, 
                                   G_rho_01_11 = G_rho_01_11, 
                                   Theta.S_0 = Theta.S_0, 
                                   Theta.S_1 = Theta.S_1,
                                   M=M, Seed=Seed, 
                                   Monotonicity=Monotonicity,
                                   Independence=Independence,
                                   HAA=HAA,
                                   Cond_ind=Cond_ind,
                                   Plots=Plots, Save.Plots=Save.Plots, 
                                   Show.Details=Show.Details)
    
    no_boots_all <- rbind(no_boots_all, no_boots)
    R2_H_all <- rbind(R2_H_all, fit.ica.bincont$R2_H)
    pi_00_all <- rbind(pi_00_all, fit.ica.bincont$pi_00)
    pi_01_all <- rbind(pi_01_all, fit.ica.bincont$pi_01)
    pi_10_all <- rbind(pi_10_all, fit.ica.bincont$pi_10)
    pi_11_all <- rbind(pi_11_all, fit.ica.bincont$pi_11)
    G_rho_01_00_all <- rbind(G_rho_01_00_all, fit.ica.bincont$G_rho_01_00)
    G_rho_01_01_all <- rbind(G_rho_01_01_all, fit.ica.bincont$G_rho_01_01)
    G_rho_01_10_all <- rbind(G_rho_01_10_all, fit.ica.bincont$G_rho_01_10)
    G_rho_01_11_all <- rbind(G_rho_01_11_all, fit.ica.bincont$G_rho_01_11)
    mu_0_00_all <- rbind(mu_0_00_all, fit.ica.bincont$mu_0_00)
    mu_0_01_all <- rbind(mu_0_01_all, fit.ica.bincont$mu_0_01)
    mu_0_10_all <- rbind(mu_0_10_all, fit.ica.bincont$mu_0_10)
    mu_0_11_all <- rbind(mu_0_11_all, fit.ica.bincont$mu_0_11)
    mu_1_00_all <- rbind(mu_1_00_all, fit.ica.bincont$mu_1_00)
    mu_1_01_all <- rbind(mu_1_01_all, fit.ica.bincont$mu_1_01)
    mu_1_10_all <- rbind(mu_1_10_all, fit.ica.bincont$mu_1_10)
    mu_1_11_all <- rbind(mu_1_11_all, fit.ica.bincont$mu_1_11)
    sigma_00_all <- rbind(sigma_00_all, fit.ica.bincont$sigma2_00_00)
    sigma_11_all <- rbind(sigma_11_all, fit.ica.bincont$sigma2_11_00)

  }
  
  fit <- list(nboots=as.numeric(no_boots_all), 
              R2_H=as.numeric(R2_H_all),
              pi_00=as.numeric(pi_00_all),
              pi_01=as.numeric(pi_01_all),
              pi_10=as.numeric(pi_10_all),
              pi_11=as.numeric(pi_11_all),
              rho_01_00=as.numeric(G_rho_01_00_all),
              rho_01_01=as.numeric(G_rho_01_01_all),
              rho_01_10=as.numeric(G_rho_01_10_all),
              rho_01_11=as.numeric(G_rho_01_11_all),
              mu_0_00=as.numeric(mu_0_00_all),
              mu_0_01=as.numeric(mu_0_01_all),
              mu_0_10=as.numeric(mu_0_10_all),
              mu_0_11=as.numeric(mu_0_11_all),
              mu_1_00=as.numeric(mu_1_00_all),
              mu_1_01=as.numeric(mu_1_01_all),
              mu_1_10=as.numeric(mu_1_10_all),
              mu_1_11=as.numeric(mu_1_11_all),
              sigma_00=as.numeric(sigma_00_all),
              sigma_11=as.numeric(sigma_11_all))   
  class(fit) <- "ICA.BinCont"
  fit
  
}
