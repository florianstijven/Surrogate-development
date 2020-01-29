BimixedCbCContCont = function(Dataset, Surr, True, Treat, Trial.ID,Min.Treat.Size=2,Alpha=0.05){
  
  R2TrialFun <- function(D){
    A <- matrix(c(D[1,4], D[2,4]),2,1)
    B <- matrix(c(D[1,1], D[1,2], D[1,2], D[2,2]),2,2)
    C <- D[4,4]
    R2.trial <- crossprod(A,solve(B,A))/C
    return(R2.trial)
  }
  R2IndFun <- function(Sigma){
    R2.ind <- Sigma[1,2]^2/(Sigma[1,1]*Sigma[2,2])
    return(R2.ind)
  }
  pdDajustment = function(D){
    EigenD <- eigen(D)
    E <- diag(EigenD$values)
    diag(E)[diag(E)<=0] <- 0.0001
    L <- EigenD$vector
    DH <- L%*%tcrossprod(E,L)
    return(DH)
  }
  Dataset = Dataset[,c(paste(substitute(Trial.ID)),paste(substitute(Treat)),paste(substitute(Surr)),
                       paste(substitute(True)))]
  Dataset = Dataset[!apply(Dataset,1,anyNA),]

  Trial.ID <- Dataset[,1]  
  Treat <- Dataset[,2]
  Surr <- Dataset[,3]
  True <- Dataset[,4]
    
  Cluster.ID = Trial.ID[order(Trial.ID)]
  n = c(table(Cluster.ID))
  N = length(n)
  Y_i = lapply(split(cbind(Surr,True)[order(Trial.ID),],Cluster.ID),function(x){matrix(x,length(x)/2,2)})
  X_i = lapply(split(cbind(1,Treat[order(Trial.ID)]),Cluster.ID),function(x){matrix(x,length(x)/2,2)})
  
  Est_i = mapply(function(x){
    Y = Y_i[[x]]
    Z = X_i[[x]]
    n.trial = min(table(Z[,2]))
    n = nrow(Y)
    if(n.trial < Min.Treat.Size){return(rep(NA,7))}
    BetaH <- tryCatch(solve(crossprod(Z),crossprod(Z,Y)),error=function(e){NULL})
    if(is.null(BetaH)){return(rep(NA,7))}
    e <- Y-Z%*%BetaH  
    SigmaH <- crossprod(e)/(n-2)
    Output <- c(c(BetaH),ks::vech(SigmaH))
    return(Output)
  },x=1:N)
  valid.outcome = !apply(Est_i,2,anyNA)
  Trial.removed = N - sum(valid.outcome)
  if(Trial.removed > 0){
    warning(paste(Trial.removed,'were removed'))
  }
  n = n[valid.outcome]
  N = length(n)
  X_i = X_i[valid.outcome]  
  BetaH_i = Est_i[1:4,valid.outcome]
  SigmaH_i = Est_i[-c(1:4),valid.outcome]
  a_i = n/sum(n)
  a_i.2 = (n-2)/sum(n-2)

    
  BetaH = apply(BetaH_i,1,weighted.mean,w=a_i)
  SigmaH = ks::invvech(apply(SigmaH_i,1,weighted.mean,w=a_i.2))
  
  R_i = mapply(function(x){
    kronecker(SigmaH,solve(crossprod(X_i[[x]])))
  },x=1:N,SIMPLIFY = F)

  N = length(a_i)
  b_i <- BetaH_i- tcrossprod(BetaH,matrix(1,N,1))  
  Sb = tcrossprod(b_i)
  
  num = 1 - 2*a_i + sum(a_i^2)
  denom = sum(num)
  DH = (Sb - Reduce('+',mapply(function(X,x){X*x},R_i,num,SIMPLIFY = F)))/denom
  DH.pd = min(eigen(DH,only.values = T)$values) > 0
  if(!DH.pd){
    DH = pdDajustment(DH)
    warning(paste("The estimate of D is non-positive definite. Adjustment for non-positive definiteness was required"))
  }
  
  Var.BetaH_i = mapply('+',R_i, MoreArgs = list(DH=DH),SIMPLIFY = F)
  VarBetaH_i.inv = mapply(solve,Var.BetaH_i,SIMPLIFY = FALSE)
  A_i.denom = solve(Reduce('+',VarBetaH_i.inv))
  A_i = mapply('%*%',list(A_i.denom),VarBetaH_i.inv,SIMPLIFY = FALSE)
  
  BetaH = apply(mapply(function(x){A_i[[x]]%*%BetaH_i[,x]},x=1:N),1,sum)
  
  VarBetaH = Reduce('+',mapply(tcrossprod,mapply('%*%',A_i,Var.BetaH_i,SIMPLIFY=FALSE),A_i,SIMPLIFY = FALSE))
  
  R2trial = R2TrialFun(DH)
  R2ind = R2IndFun(SigmaH)

  R2ind.sd <- sqrt((4 * R2ind * ((1 - R2ind)^2))/(N - 3))
  R2ind.lb <- max(0, R2ind + qnorm(Alpha/2) * R2ind.sd)
  R2ind.ub <- min(1, R2ind + qnorm(1 - Alpha/2) * R2ind.sd)
  Indiv.R2 <- data.frame(cbind(R2ind, R2ind.sd, R2ind.lb, R2ind.ub))
  colnames(Indiv.R2) <- c("R2 Indiv", "Standard Error", 
                          "CI lower limit", "CI upper limit")
  
  R2trial.sd <- sqrt((4 * R2trial * (1 - R2trial)^2)/(N - 3))
  R2trial.lb <- max(0, R2trial + qnorm(Alpha/2) * (R2trial.sd))
  R2trial.ub <- min(1, R2trial + qnorm(1 - Alpha/2) * (R2trial.sd))
  Trial.R2 <- data.frame(cbind(R2trial, R2trial.sd,R2trial.lb, R2trial.ub))
  colnames(Trial.R2) <- c("R2 Trial", "Standard Error", 
                          "CI lower limit", "CI upper limit")
  rownames(Trial.R2) = rownames(Indiv.R2) <- c(" ")
  
  
  Obs.Per.Trial = cbind(unique(Cluster.ID)[valid.outcome],t(mapply(function(X){table(X[,2])},X_i)),n)
  colnames(Obs.Per.Trial) = c('Trial', 'Number.cont.Treat', 'Number.exp.Treat', 'Obs.per.trial')
  Obs.Per.Trial = as.data.frame(Obs.Per.Trial)
  
  Fixed.Effects = cbind(BetaH,sqrt(diag(VarBetaH)))[c(3,1,4,2),]
  rownames(Fixed.Effects) = c('mu_S','mu_T','alpha','beta')
  colnames(Fixed.Effects) = c('estimate','standard error')
  Fixed.Effects = as.data.frame(Fixed.Effects)
  
  DH = DH[c(3,1,4,2),c(3,1,4,2)]
  colnames(DH) = rownames(DH) = c('mu_S_i','mu_T_i','a_i','b_i')
  SigmaH = SigmaH
  colnames(SigmaH) = rownames(SigmaH) = c('S_ij','T_ij')
  
  Output = list(Obs.Per.Trial=Obs.Per.Trial, Trials.removed=Trial.removed, Fixed.Effects=Fixed.Effects,
                Trial.R2=Trial.R2,Indiv.R2=Indiv.R2,D=DH,D.pd=DH.pd,Sigma=SigmaH,Call = match.call())
  class(Output) <- "BimixedCbCContCont"
  return(Output)
}

summary.BimixedCbCContCont = function (object, ..., Object){
  if (missing(Object)) {
    Object <- object
  }
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  if(Object$Trials.removed > 0){
    cat("\nTotal number of trials removed: ", Object$Trials.removed)
  }
  
  cat("\nTotal number of patients: ", sum(Object$Obs.Per.Trial$Obs.per.trial))
  cat("\nM(SD) patients per trial: ", format(round(mean((Object$Obs.Per.Trial$Obs.per.trial)), 
                                                   4), nsmall = 4), " (", format(round(sd((Object$Obs.Per.Trial$Obs.per.trial)), 
                                                                                       4), nsmall = 4), ")", "  [min: ", min((Object$Obs.Per.Trial$Obs.per.trial)), 
      "; max: ", max((Object$Obs.Per.Trial$Obs.per.trial)), 
      "]", sep = "")
  cat("\nTotal number of patients in experimental treatment group: ", 
      sum(Object$Obs.Per.Trial$Number.exp.Treat), "\nTotal number of patients in control treatment group: ", 
      sum(Object$Obs.Per.Trial$Number.cont.Treat))
  
  cat("\n\n\n# Meta-analytic results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  cat('Fixed effects summary')
  cat("\n")
  print(format(round(Object$Fixed.Effects,4),nsmall=4))
  cat("\n")
  cat('Surrogacy at the trial level')
  cat("\n")
  print(format(round(Object$Trial.R2, 4), nsmall = 4))
  cat("\n")
  cat('Surrogacy at the individual level')
  cat("\n")
  print(format(round(Object$Indiv.R2, 4), nsmall = 4))
  cat("\n")
  if(Object$D.pd != TRUE){
    cat('Adjustment for non-positive definiteness was needed to estimate D')
    cat("\n")
  }
}
