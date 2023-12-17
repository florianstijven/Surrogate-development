MufixedContCont.MultS <- function(Dataset, Endpoints=True~Surr.1+Surr.2, 
 Treat="Treat", Trial.ID="Trial.ID", Pat.ID="Pat.ID", 
 Model=c("Full"), Weighted=TRUE, Min.Trial.Size=2, Alpha=.05, 
 Number.Bootstraps=0, Seed=123){          
  
  ci.R2 <- resid.Fitted.Model. <- NULL
  
  if ((Model==c("Full") | Model==c("Reduced"))==FALSE) {stop ("The specification of the Model=c(\"...\") argument of the call is incorrect. Use either Model=c(\"Full\") or Model=c(\"Reduced\")).")}     
  
  Surr.1 <- endpoint <- outcome <- NULL
  Vars <- all.vars(Endpoints)
  True <- Vars[1]
  Surr <- Vars[-1]
  Dataset.Wide <- na.exclude(Dataset[,c(Surr, True, Treat, Trial.ID, Pat.ID)])
  names(Dataset.Wide) <- c(paste("Surr", 1:1000, sep=".")[1:length(Surr)], "True", "Treat", "Trial.ID", "Pat.ID")
  Dataset.Wide$Trial.ID <- as.character(Dataset.Wide$Trial.ID) # LS
      # Remove trials with less patients than specified in Min.Trial.Size
  Select.Trial.IDs <- names(table(Dataset.Wide$Trial.ID))[table(Dataset.Wide$Trial.ID)>=Min.Trial.Size]
  Dataset.Wide <- Dataset.Wide[Dataset.Wide$Trial.ID %in% Select.Trial.IDs,]
  Obs.per.trial <- table(Dataset.Wide$Trial.ID)
  N.total <- dim(Dataset.Wide)[1]
  data <- na.exclude(tidyr::gather(Dataset.Wide, endpoint, outcome, c(Surr.1):True))
  data <- data[order(c(data$Pat.ID)),]
  
  # stage 1
  if (Model==c("Full")){
    Fitted.Model <- nlme::gls(outcome ~ -1 + as.factor(Trial.ID):Treat:as.factor(endpoint) + as.factor(Trial.ID):as.factor(endpoint),
                                   correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), 
                                   weights=nlme::varIdent(form=~1|endpoint), data=data, control = list(msVerbose = FALSE, 
                                   optimizer = "nlm", niterEM = 25, msMaxIter=5000)) 
    summary(Fitted.Model)
    Fixed.eff <- data.frame(matrix(coef(Fitted.Model), nrow = dim(Obs.per.trial)))
    names(Fixed.eff) <- c(paste("Intercept.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Intercept.T", 
                          paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
    Results.Stage.1 <- data.frame(Obs.per.trial, Fixed.eff, stringsAsFactors = TRUE)
    colnames(Results.Stage.1)[1:2] <- c("Trial", "Obs.per.trial")
    rownames(Results.Stage.1) <- NULL 
    D.equiv <- var(Results.Stage.1[,-c(1:2)])
    Residuals.Stage.1 <- as.numeric(residuals(Fitted.Model))
  }
  
  if (Model==c("Reduced")){
    Fitted.Model <- nlme::gls(outcome ~ -1 + as.factor(Trial.ID):Treat:as.factor(endpoint) + as.factor(endpoint), 
              correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), 
              weights=nlme::varIdent(form=~1|endpoint), data=data, control = list(msVerbose = FALSE, 
                  optimizer = "nlm", niterEM = 25, msMaxIter=5000))
    summary(Fitted.Model)
    Fixed.eff <- data.frame(matrix(coef(Fitted.Model)[-c(1:length(Vars))],  # coef 1 to length(vars) are the overall (i.e., not trial-specific) intercepts for S1, S2, ..., SK, T. These are not needed further on
                                   nrow = dim(Obs.per.trial)))
    names(Fixed.eff) <- c(paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
    Results.Stage.1 <- data.frame(Obs.per.trial, Fixed.eff, stringsAsFactors = TRUE)
    colnames(Results.Stage.1)[1:2] <- c("Trial", "Obs.per.trial")
    rownames(Results.Stage.1) <- NULL 
    D.equiv <- var(Results.Stage.1[,-c(1:2)])
  }
  
  # Stage 2
  Results.Stage.1.For.Stage.2.Analysis <- Results.Stage.1[,-c(1:2)]
  Select.Cols <- names(Results.Stage.1.For.Stage.2.Analysis)!="Intercept.T"   # Exclude intercept T, not needed for prediction Treatment effect T further on
  Results.Stage.1.For.Stage.2.Analysis <- Results.Stage.1.For.Stage.2.Analysis[,Select.Cols]
  
  if (Model==c("Full")){
    if (Weighted == FALSE){
      Results.Stage.2 <- lm(Treatment.T ~ ., data=Results.Stage.1.For.Stage.2.Analysis)
      }
    if (Weighted == TRUE){
      Results.Stage.2 <- lm(Treatment.T ~ ., data=Results.Stage.1.For.Stage.2.Analysis, weights=Results.Stage.1$Obs.per.trial)
      }
    }
  if (Model==c("Reduced")){
    if (Weighted == FALSE){
      Results.Stage.2 <- lm(Treatment.T ~ ., data=Results.Stage.1.For.Stage.2.Analysis)
      }
    if (Weighted == TRUE){
      Results.Stage.2 <- lm(Treatment.T ~ ., data=Results.Stage.1.For.Stage.2.Analysis, weights=Results.Stage.1$Obs.per.trial)
      }
    }

  # R2 trial
  Trial.R2.value <- as.numeric(summary(Results.Stage.2)[c("r.squared")])
  # R2 trial
  Trial.R2.Adj.value <- as.numeric(summary(Results.Stage.2)[c("adj.r.squared")])
  
    
  Boot.R2.Trial <- Boot.R2.Trial.Adj <- Boot.R2.Indiv <- NULL
  
  # Bootstrap CI for R2 trial and R2 indiv
  if (Number.Bootstraps  > 0){
    
    set.seed(Seed)
    for (z in 1: Number.Bootstraps){
      Pat.IDs <- unique(data$Pat.ID)
      Select.Pat.Bootstrap <- sample(Pat.IDs, size = N.total, replace = TRUE)
      data.boot <- data[data$Pat.ID %in% Select.Pat.Bootstrap,]
      
      if (exists("Fitted.Model.Boot")){rm(Fitted.Model.Boot)} 
      # stage 1 bootstrap
      if (Model==c("Full")){
        try(Fitted.Model.Boot <- nlme::gls(outcome ~ -1 + as.factor(Trial.ID):Treat:as.factor(endpoint) + as.factor(Trial.ID):as.factor(endpoint),
                                       correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), 
                                       weights=nlme::varIdent(form=~1|endpoint), data=data.boot, control = list(msVerbose = FALSE, 
                                                                                                                optimizer = "nlm", niterEM = 25, msMaxIter=5000)), silent = TRUE) 
        if (exists("Fitted.Model.Boot")){
        Fixed.Eff.Boot <- data.frame(matrix(coef(Fitted.Model.Boot), nrow = dim(Obs.per.trial)))
        names(Fixed.Eff.Boot) <- c(paste("Intercept.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Intercept.T", 
                                   paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
        Results.Stage.1.Boot <- data.frame(Obs.per.trial, Fixed.Eff.Boot, stringsAsFactors = TRUE)
        colnames(Results.Stage.1.Boot)[1:2] <- c("Trial", "Obs.per.trial")
        rownames(Results.Stage.1.Boot) <- NULL 
        D.equiv.Boot <- var(Results.Stage.1.Boot[,-c(1:2)])
        Residuals.Stage.1.Boot <- as.numeric(residuals(Fitted.Model.Boot))
        }
      }
      
      if (Model==c("Reduced")){
        try(Fitted.Model.Boot <- nlme::gls(outcome ~ -1 + as.factor(Trial.ID):Treat:as.factor(endpoint) + as.factor(endpoint), 
                                       correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), 
                                       weights=nlme::varIdent(form=~1|endpoint), data=data.boot, control = list(msVerbose = FALSE, 
                                                                                                                optimizer = "nlm", niterEM = 25, msMaxIter=5000)), silent=TRUE)
        if (exists("Fitted.Model.Boot")){
        Fixed.Eff.Boot <- data.frame(matrix(coef(Fitted.Model.Boot)[-c(1:length(Vars))],  # coef 1 to length(vars) are the overall (i.e., not trial-specific) intercepts for S1, S2, ..., SK, T. These are not needed further on
                                            nrow = dim(Obs.per.trial)))
        names(Fixed.Eff.Boot) <- c(paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
        Results.Stage.1.Boot <- data.frame(Obs.per.trial, Fixed.Eff.Boot, stringsAsFactors = TRUE)
        colnames(Results.Stage.1.Boot)[1:2] <- c("Trial", "Obs.per.trial")
        rownames(Results.Stage.1.Boot) <- NULL 
        D.equiv.Boot <- var(Results.Stage.1.Boot[,-c(1:2)])
        }
      }
      
      if (exists("Fitted.Model.Boot")){
      # Stage 2
      Results.Stage.1.Boot.For.Stage.2.Analysis.Boot <- Results.Stage.1.Boot[,-c(1:2)]
      Select.Cols <- names(Results.Stage.1.Boot.For.Stage.2.Analysis.Boot)!="Intercept.T"   # Exclude intercept T, not needed for prediction Treatment effect T further on
      Results.Stage.1.Boot.For.Stage.2.Analysis.Boot <- Results.Stage.1.Boot.For.Stage.2.Analysis.Boot[,Select.Cols]
      
      if (Model==c("Full")){
        if (Weighted == FALSE){
          Results.Stage.2.Boot <- lm(Treatment.T ~ ., data=Results.Stage.1.Boot.For.Stage.2.Analysis.Boot)
        }
        if (Weighted == TRUE){
          Results.Stage.2.Boot <- lm(Treatment.T ~ ., data=Results.Stage.1.Boot.For.Stage.2.Analysis.Boot, weights=Results.Stage.1.Boot$Obs.per.trial)
        }
      }
      if (Model==c("Reduced")){
        if (Weighted == FALSE){
          Results.Stage.2.Boot <- lm(Treatment.T ~ ., data=Results.Stage.1.Boot.For.Stage.2.Analysis.Boot)
        }
        if (Weighted == TRUE){
          Results.Stage.2.Boot <- lm(Treatment.T ~ ., data=Results.Stage.1.Boot.For.Stage.2.Analysis.Boot, weights=Results.Stage.1.Boot$Obs.per.trial)
        }
      }
      
      # R2 trial
      Trial.R2.value.Boot <- as.numeric(summary(Results.Stage.2.Boot)[c("r.squared")])
      # R2 trial
      Trial.R2.Adj.value.Boot <- as.numeric(summary(Results.Stage.2.Boot)[c("adj.r.squared")])
      
      # R2 indiv
      cors.Boot <- nlme::corMatrix(Fitted.Model.Boot$modelStruct$corStruct)[[1]]
      rho.sigma.ST.Boot <- matrix(cors.Boot[,c(length(Vars))][-length(cors.Boot[,c(length(Vars))])], ncol = 1)#; rho.sigma.ST
      rho.sigma.SS.Boot <- cors.Boot[c(1:c(length(Vars)-1)), c(1:c(length(Vars)-1))]
      R2ind.Boot <- t(rho.sigma.ST.Boot) %*% solve(rho.sigma.SS.Boot) %*% rho.sigma.ST.Boot
      } 
      
      if (exists("Fitted.Model.Boot")==FALSE){
        Trial.R2.value.Boot <- NA
        Trial.R2.Adj.value.Boot <- NA
        R2ind.Boot <- NA
      }
      
      Boot.R2.Trial <- na.exclude(c(Boot.R2.Trial, Trial.R2.value.Boot))
      Boot.R2.Trial.Adj <- na.exclude(c(Boot.R2.Trial.Adj, Trial.R2.Adj.value.Boot))
      Boot.R2.Indiv <- na.exclude(c(Boot.R2.Indiv, R2ind.Boot))
      
    }   
  }# END bootstrap
  
  
  # CIs for R2 trial 
  
  # CI Lee
  N.trial <- length(Obs.per.trial)
  Num.Preds.Stage.2 <- c(dim(Results.Stage.1.For.Stage.2.Analysis)[2]-1)
  Trial.R2.lb.Lee <- MBESS::ci.R2(R2 = Trial.R2.value, conf.level = 1-Alpha, N=N.trial, p = Num.Preds.Stage.2)$Lower.Conf.Limit.R2
  Trial.R2.ub.Lee <- MBESS::ci.R2(R2 = Trial.R2.value, conf.level = 1-Alpha, N=N.trial, p = Num.Preds.Stage.2)$Upper.Conf.Limit.R2
  Trial.R2 <- data.frame(cbind(Trial.R2.value, Trial.R2.lb.Lee, Trial.R2.ub.Lee), stringsAsFactors = TRUE)
  colnames(Trial.R2) <- c("R2 Trial", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ")
  Trial.R2.Lee.Method <- Trial.R2
  
  # CI bootstrap
  Trial.R2 <- data.frame(cbind(Trial.R2.value, quantile(Boot.R2.Trial, probs = c(Alpha/2)), 
                               quantile(Boot.R2.Trial, probs = c(1-Alpha/2))), stringsAsFactors = TRUE)
  colnames(Trial.R2) <- c("R2 Trial", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ")
  Trial.R2.Boot.Method <- Trial.R2
  
  
  # CIs for adjusted R2 trial 
  
  # CI Lee
  Trial.R2.lb.Lee <- MBESS::ci.R2(R2 = Trial.R2.Adj.value, conf.level = 1-Alpha, N=N.trial, p = Num.Preds.Stage.2)$Lower.Conf.Limit.R2
  Trial.R2.ub.Lee <- MBESS::ci.R2(R2 = Trial.R2.Adj.value, conf.level = 1-Alpha, N=N.trial, p = Num.Preds.Stage.2)$Upper.Conf.Limit.R2
  Trial.R2 <- data.frame(cbind(Trial.R2.Adj.value, Trial.R2.lb.Lee, Trial.R2.ub.Lee), stringsAsFactors = TRUE)
  colnames(Trial.R2) <- c("R2 Trial", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ")
  Trial.R2.Adj.Lee.Method <- Trial.R2
  
  # CI bootstrap
  Trial.R2 <- data.frame(cbind(Trial.R2.Adj.value, quantile(Boot.R2.Trial.Adj, probs = c(Alpha/2)), 
                               quantile(Boot.R2.Trial.Adj, probs = c(1-Alpha/2))), stringsAsFactors = TRUE)
  colnames(Trial.R2) <- c("R2 Trial", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ")
  Trial.R2.Adj.Boot.Method <- Trial.R2
  
  
  

  # Individual-Level Surrogacy and CI
  cors <- nlme::corMatrix(Fitted.Model$modelStruct$corStruct)[[1]]  
  # CHECK
  Resid.Stage.1 <- data.frame(resid(Fitted.Model), data)
  Resid.Stage.1$Patient.ID <- paste(Resid.Stage.1$Trial.ID, Resid.Stage.1$Pat.ID)
  For.Sigma <- tidyr::spread(Resid.Stage.1[c("endpoint", "resid.Fitted.Model.", "Patient.ID")], key=endpoint, value=resid.Fitted.Model.)
  #cors <- cor(For.Sigma[,-1])
    # Fit lm on residuals. Useful for LR test to compare R2_indiv different models
  For.Sigma <- For.Sigma[,-1]
  Model.R2.Indiv <- lm(data = For.Sigma, formula = True~.)

  # CI based on Lee
  rho.sigma.ST <- matrix(cors[,c(length(Vars))][-length(cors[,c(length(Vars))])], ncol = 1)#; rho.sigma.ST
  rho.sigma.SS <- cors[c(1:c(length(Vars)-1)), c(1:c(length(Vars)-1))]
  R2ind <- t(rho.sigma.ST) %*% solve(rho.sigma.SS) %*% rho.sigma.ST
  Num.Preds <- (length(Vars)-1)
  Indiv.R2.lb.Lee <- MBESS::ci.R2(R2 = R2ind, conf.level = 1-Alpha, N=N.total, p = Num.Preds)$Lower.Conf.Limit.R2
  Indiv.R2.ub.Lee <- MBESS::ci.R2(R2 = R2ind, conf.level = 1-Alpha, N=N.total, p = Num.Preds)$Upper.Conf.Limit.R2
  Indiv.R2 <- data.frame(cbind(R2ind, Indiv.R2.lb.Lee, Indiv.R2.ub.Lee), stringsAsFactors = TRUE)
  colnames(Indiv.R2) <- c("R2 Indiv", "CI lower limit", "CI upper limit")
  rownames(Indiv.R2) <- c(" ")
  Indiv.R2.Lee <- Indiv.R2
  
  # CI based on Bootstrap
  Indiv.R2 <- data.frame(cbind(R2ind, quantile(Boot.R2.Indiv, probs = c(Alpha/2)), 
                               quantile(Boot.R2.Indiv, probs = c(1-Alpha/2))), stringsAsFactors = TRUE)
  colnames(Indiv.R2) <- c("R2 Indiv", "CI lower limit", "CI upper limit")
  rownames(Indiv.R2) <- c(" ")
  Indiv.R2.Boot <- Indiv.R2
  
  Min.Eigen.CorResid <- min(eigen(cors)$values)    # lowest eigenvalue
  
  if (Min.Eigen.CorResid <= 0) warning(paste("The R-square Individual estimate may be invalid, because its calculation is based on a non-positive definite correlation matrix. "))
  
  fit <- 
    list(Data.Analyze=Dataset.Wide, Obs.Per.Trial=Obs.per.trial, Results.Stage.1=Results.Stage.1, 
         Results.Stage.2=Results.Stage.2,  
        Trial.R2.Lee=Trial.R2.Lee.Method, Trial.R2.Boot=Trial.R2.Boot.Method,
        Trial.R2.Adj.Lee=Trial.R2.Adj.Lee.Method, Trial.R2.Adj.Boot=Trial.R2.Adj.Boot.Method,
        Indiv.R2.Lee=Indiv.R2.Lee, Indiv.R2.Boot=Indiv.R2.Boot, Fitted.Model.Stage.1=Fitted.Model, 
        Model.R2.Indiv=Model.R2.Indiv,D.Equiv=D.equiv, 
        Call=match.call())   
  
  class(fit) <- "MufixedContCont.MultS"
  fit
}
