MumixedContCont.MultS <- function(Dataset, Endpoints=True~Surr.1+Surr.2,  
    Treat="Treat", Trial.ID="Trial.ID", Pat.ID="Pat.ID", 
    Model=c("Full"), Min.Trial.Size=2, Alpha=.05, Opt="nlminb"){
  
  if ((Model==c("Full") | Model==c("Reduced"))==FALSE){stop ("The specification of the Model=c(\"...\") argument of the call is incorrect. Use either Model=c(\"Full\") or Model=c(\"Reduced\").")}     
  
  Surr.1 <- endpoint <- outcome <- lmeControl <- ci.R2 <- R2ind <- NULL
  Vars <- all.vars(Endpoints)
  True <- Vars[1]
  Surr <- Vars[-1]
  Dataset.Wide <- na.exclude(Dataset[,c(Surr, True, Treat, Trial.ID, Pat.ID)])
  names(Dataset.Wide) <- c(paste("Surr", 1:1000, sep=".")[1:length(Surr)], "True", "Treat", "Trial.ID", "Pat.ID")
  Dataset.Wide$Trial.ID <- as.character(Dataset.Wide$Trial.ID) # LS
  # Remove trials with less patients than specified in Min.Trial.Size
  Select.Trial.IDs <- names(table(Dataset.Wide$Trial.ID))[table(Dataset.Wide$Trial.ID)>=Min.Trial.Size]
  Dataset.Wide <- Dataset.Wide[Dataset.Wide$Trial.ID %in% Select.Trial.IDs,]
  
  # Check coding
  if (length(unique(Dataset.Wide$Treat))!=2) stop("Please make sure that the treatment variable has only 2 levels.")
  if ((sort(unique(Dataset.Wide$Treat))[1]==c(-0.5)) & (sort(unique(Dataset.Wide$Treat))[2]==c(0.5))){
    Dataset.Wide$Treat <- Dataset.Wide$Treat+.5}
  if (((sort(unique(Dataset.Wide$Treat))[1]==c(0) & sort(unique(Dataset.Wide$Treat))[2]==c(1))==FALSE) & 
      ((sort(unique(Dataset.Wide$Treat))[1]==c(-1) & sort(unique(Dataset.Wide$Treat))[2]==c(1))==FALSE))
    stop("Please make sure that the treatment is either coded as control = -1 and experimental = 1, or as control = 0 and experimental = 1.")
  
  
  Obs.per.trial <- table(Dataset.Wide$Trial.ID)
  N.total <- dim(Dataset.Wide)[1]
  data <- na.exclude(tidyr::gather(Dataset.Wide, endpoint, outcome, c(Surr.1):True))
  
  if (Model==c("Full")){
     ctrl <- nlme::lmeControl(opt=Opt)
     Fitted.Model <- nlme::lme(outcome~ -1 + as.factor(endpoint):Treat + as.factor(endpoint), 
                          random=~ -1 + as.factor(endpoint) + as.factor(endpoint):Treat|as.factor(Trial.ID),
                          correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), data=data,
                          weights=nlme::varIdent(form=~1|endpoint), control=ctrl) 
     Fixed.Effect.Pars <- data.frame(summary(Fitted.Model)$coef$fixed)
     row.names(Fixed.Effect.Pars) <- c(paste("Intercept.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Intercept.T", 
                  paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
     names(Fixed.Effect.Pars) <- "Estimate"
     
     Random.effect.pars <- data.frame(summary(Fitted.Model)$coef$random)
     names(Random.effect.pars) <- c(paste("Intercept.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Intercept.T", 
                                       paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
     Residuals <- residuals(Fitted.Model, type='response')
     }
  
  if (Model==c("Reduced")){
    
    ctrl <- nlme::lmeControl(opt=Opt)
    Fitted.Model <- nlme::lme(outcome~ -1 + as.factor(endpoint):Treat + as.factor(endpoint), 
                         random=~ -1 + as.factor(endpoint):Treat|as.factor(Trial.ID),
                         correlation = nlme::corSymm(form=~1|as.factor(Trial.ID)/as.factor(Pat.ID)), data=data,
                         weights=nlme::varIdent(form=~1|endpoint), control = ctrl) 
    
    Fixed.Effect.Pars <- data.frame(summary(Fitted.Model)$coef$fixed)
    row.names(Fixed.Effect.Pars) <- c(paste("Intercept.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Intercept.T", 
                                      paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
    names(Fixed.Effect.Pars) <- "Estimate"
    
    Random.effect.pars <- data.frame(summary(Fitted.Model)$coef$random)
    names(Random.effect.pars) <- c(paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
    Residuals <- residuals(Fitted.Model, type='response')
   }
  
  
  # Trial-level surrogacy estimates
  if (Model==c("Full")){
    
    D <- matrix(as.matrix(Fitted.Model$modelStruct$reStruct$`as.factor(Trial.ID)`), ncol=length(Vars)*2)
    rownames(D) <- colnames(D) <- c(paste("Intercept.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Intercept.T", 
                                    paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
    Min.Eigen.D <- min(eigen(D)$values)
    if (Min.Eigen.D <= 0) warning(paste("The R-square Trial estimate may be invalid, because its calculation is based on a non-positive definite covariance matrix"))
    singular  <- svd(D)$d
    Cond.Number.D.Matrix <- max(singular)/min(singular)
    if (Cond.Number.D.Matrix > 100) warning(paste("The R-square Trial estimate may not be thrustworthy, as the conditioning number of the D matrix is high and equals",
                                                  deparse(Cond.Number.D.Matrix)))  

    Select.Rows <- grepl(x = row.names(D), ".T")==FALSE
    Select.Col <- grepl(x = row.names(D), "Treatment.T")==TRUE # LS row.names
    Treatment.T.Indicator <- grepl(x = row.names(D), "Treatment.T")==TRUE # LS row.names
    D.ST <- matrix(D[Select.Rows, Select.Col], ncol = 1); D.ST
    Select.Rows <- Select.Cols <- grepl(x = row.names(D), ".S")==TRUE
    D.SS <- D[Select.Rows, Select.Cols]
    D.TT <- D[Treatment.T.Indicator, Treatment.T.Indicator]
    Trial.R2.value <- c(t(D.ST) %*% solve(D.SS) %*% D.ST)/D.TT
    }
  
  if (Model==c("Reduced")){
    
    D <- matrix(as.matrix(Fitted.Model$modelStruct$reStruct$`as.factor(Trial.ID)`), ncol=length(Vars))
    rownames(D) <- colnames(D) <- c(paste("Treatment.S", 1:1000, sep="")[c(1:length(Vars)-1)], "Treatment.T")
    Min.Eigen.D <- min(eigen(D)$values)
    if (Min.Eigen.D <= 0) warning(paste("The R-square Trial estimate may be invalid, because its calculation is based on a non-positive definite covariance matrix"))
    singular  <- svd(D)$d
    Cond.Number.D.Matrix <- max(singular)/min(singular)
    if (Cond.Number.D.Matrix > 100) warning(paste("The R-square Trial estimate may not be thrustworthy, as the conditioning number of the D matrix is high and equals",
                                                  deparse(Cond.Number.D.Matrix)))  
    Select.Rows <- grepl(x = row.names(D), ".T")==FALSE
    Select.Col <- grepl(x = row.names(D), "Treatment.T")==TRUE # LS row.names
    Treatment.T.Indicator <- grepl(x = row.names(D), "Treatment.T")==TRUE # LS row.names
    D.ST <- matrix(D[Select.Rows, Select.Col], ncol = 1); D.ST
    Select.Rows <- Select.Cols <- grepl(x = row.names(D), ".S")==TRUE
    D.SS <- D[Select.Rows, Select.Cols]
    D.TT <- D[Treatment.T.Indicator, Treatment.T.Indicator]
    Trial.R2.value <- c(t(D.ST) %*% solve(D.SS) %*% D.ST)/D.TT
      }
  
  N.trial <- length(Obs.per.trial)
  Num.Preds <- length(D.ST)
  
  # CI Lee
  Trial.R2.lb.Lee <- MBESS::ci.R2(R2 = Trial.R2.value, conf.level = 1-Alpha, N=N.trial, p = Num.Preds)$Lower.Conf.Limit.R2
  Trial.R2.ub.Lee <- MBESS::ci.R2(R2 = Trial.R2.value, conf.level = 1-Alpha, N=N.trial, p = Num.Preds)$Upper.Conf.Limit.R2
  Trial.R2 <- data.frame(cbind(Trial.R2.value, Trial.R2.lb.Lee, Trial.R2.ub.Lee), stringsAsFactors = TRUE)
  colnames(Trial.R2) <- c("R2 Trial", "CI lower limit", "CI upper limit")
  rownames(Trial.R2) <- c(" ")
  Trial.R2.Lee.Method <- Trial.R2
  
  # R2 indiv
  cors <- nlme::corMatrix(Fitted.Model$modelStruct$corStruct)[[1]]
  
  # CI based on Lee
  rho.sigma.ST <- matrix(cors[,c(length(Vars))][-length(cors[,c(length(Vars))])], ncol = 1)#; rho.sigma.ST
  rho.sigma.SS <- cors[c(1:c(length(Vars)-1)), c(1:c(length(Vars)-1))]
  R2ind <- t(rho.sigma.ST) %*% solve(rho.sigma.SS) %*% rho.sigma.ST
  Indiv.R2.lb.Lee <- MBESS::ci.R2(R2 = R2ind, conf.level = 1-Alpha, N=N.total, p = Num.Preds)$Lower.Conf.Limit.R2
  Indiv.R2.ub.Lee <- MBESS::ci.R2(R2 = R2ind, conf.level = 1-Alpha, N=N.total, p = Num.Preds)$Upper.Conf.Limit.R2
  Indiv.R2 <- data.frame(cbind(R2ind, Indiv.R2.lb.Lee, Indiv.R2.ub.Lee), stringsAsFactors = TRUE)
  colnames(Indiv.R2) <- c("R2 Indiv", "CI lower limit", "CI upper limit")
  rownames(Indiv.R2) <- c(" ")
  Indiv.R2.Lee <- Indiv.R2

  
  Min.Eigen.CorResid <- min(eigen(cors)$values)    # lowest eigenvalue
  if (Min.Eigen.CorResid <= 0) warning(paste("The R-square Individual estimate may be invalid, because its calculation is based on a non-positive definite correlation matrix. "))
  singular  <- svd(cors)$d
  Cond.Number.Sigma.Matrix <- max(singular)/min(singular)
  
  fit <- 
    list(Data.Analyze=Dataset.Wide, Obs.Per.Trial=Obs.per.trial, 
        Fixed.Effect=Fixed.Effect.Pars, Random.Effect=Random.effect.pars, 
        Trial.R2.Lee=Trial.R2.Lee.Method, Indiv.R2.Lee=Indiv.R2.Lee, 
        D=D, Cond.Number.D.Matrix=Cond.Number.D.Matrix, Cond.Number.Sigma.Matrix=Cond.Number.Sigma.Matrix, 
        Fitted.Model=Fitted.Model,
        Call=match.call())
  
  class(fit) <- "MumixedContCont.MultS"
  
  fit
  
}
