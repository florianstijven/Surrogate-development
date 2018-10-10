#require(ks)
#require(rootSolve)


ICA.ContCont.MultS.PC = function(M=1000,N,Sigma,Seed=123,Show.Progress=FALSE){
 
  # LIST OF SUPPORT FUNCTIONS
  
  # function to evaluate whether a correlation matrix is PD based on the eigenvalues
  is.PD = function(X,tol=1e-8){
    X[is.na(X)] = 0
    min.lambda = min(eigen(X,only.values=T,symmetric=T)$values)
    return(min.lambda > tol)
  }
  # function to estimate a single correlation r_[j,j+k] based on partial correlations
  r.random = function(j,k,R){
    d = ncol(R)
    r1 = R[j,(j+1):(j+k-1)]
    r3 = R[(j+1):(j+k-1),j+k]
    R2.inv = R[(j+1):(j+k-1),(j+1):(j+k-1)]
    
    R2 = solve(R2.inv)
    r.c = extraDistr::rnsbeta(1,1+0.5*(d-1-k),1+0.5*(d-1-k),-1,1)
    D2 = (1 - tcrossprod(crossprod(r1,R2),r1) )*(1 - tcrossprod(crossprod(r3,R2),r3))
    if(D2 < 0 & D2 > -1e-8){D2 = 0} 
    r = tcrossprod(crossprod(r1,R2),r3) + r.c*sqrt(D2)
    return(r)
  }
  # function to generate a random correlation matrix based on partial correlations
  Correlation.matrix.PC = function(R,Range=c(-1,1)){
   # require(extraDistr)
    Rf = R
    Rf[is.na(Rf)] = 0
    d = ncol(R)
    j.ind = do.call('c',lapply(1:(d-1),function(x){x:1}))
    k.ind = do.call('c',lapply(1:(d-1),function(x){1:x}))
    for(i in 1:(length(j.ind))){
      j = j.ind[i]
      k = k.ind[i]
      if(is.na(R[j,j+k])){
        if(k==1){
          R[j,j+k] = R[j+k,j] = extraDistr::rnsbeta(1,d/2,d/2,Range[1],Range[2])
        }else{
          R[j,j+k] = R[j+k,j]  = r.random(j,k,R)
          
        }
        if(is.nan(R[j,j+k])){stop('error')}
      }
    }
    return(R)
  }
  # function to compute the ICA
  MutivarICA.fun = function(R,Sigma,N){
    d = nrow(R)
    sdMat = diag(sqrt(Sigma))
    rtn = sdMat %*% R %*% t(sdMat)
    var_diff <- function(cov_mat) {
      cov_val <- cov_mat[1, 1] + cov_mat[2, 2] - (2 * 
                                                    cov_mat[1, 2])
      fit <- c(cov_val)
      fit
    }
    cov_2_diffs <- function(cov_mat) {
      cov_val <- (cov_mat[2, 2] - cov_mat[1, 2]) - 
        (cov_mat[2, 1] - cov_mat[1, 1])
      fit <- c(cov_val)
      fit
    }
    A <- matrix(var_diff(cov_mat = rtn[1:2, 1:2]), nrow = 1)
    B <- NULL
    Aantal <- (dim(R)[1] - 2)/2
    rtn_part <- rtn[c(3:dim(rtn)[1]), c(1, 2)]
    for (z in 1:Aantal) {
      cov_mat_here <- rtn_part[c((z * 2) - 1, z * 2), 
                               c(1:2)]
      B <- rbind(B, cov_2_diffs(cov_mat_here))
    }
    Dim <- dim(R)[1]
    Sub_mat_var <- rtn[c(3:Dim), c(3:Dim)]
    C <- matrix(NA, nrow = Aantal, ncol = Aantal)
    for (l in 1:Aantal) {
      for (k in 1:Aantal) {
        Sub_mat_var_hier <- Sub_mat_var[c((k * 2) - 
                                            1, k * 2), c((l * 2) - 1, l * 2)]
        C[k, l] <- cov_2_diffs(cov_mat = Sub_mat_var_hier)
      }
    }
    Delta <- cbind(rbind(A, B), rbind(t(B), C))
    ICA <- (t(B) %*% solve(C) %*% B)/A
    Adj.ICA <- 1 - (1 - ICA) * ((N - 1)/(N - Aantal - 
                                           1))
    return(c(ICA, Adj.ICA))
  }
  
  
  # Start of function ICA.ContCont.MultS.PC
  
  set.seed(Seed)
  d = nrow(Sigma)[1]
  Vars = diag(Sigma)
  IND =  ks::vec(matrix(1:d,ncol=2,byrow = T),byrow = F) 
  Sigma = Sigma[IND,IND]
  R = cov2cor(Sigma)
  if(!is.PD(R)){
    alpha  = uniroot(function(alpha,R,Rfixed,tol=0.0001){
      if(anyNA(R)){
        R[is.na(R)] = 0
      }
      f = alpha*R + (1-alpha)*Rfixed
      min(eigen(f)$values) - tol
    },c(0,1),R=R,Rfixed=diag(d),tol = 1e-8)$root
    R = R*alpha + (1-alpha)*diag(d)
    warning(paste('The initial correlation matrix is not PD. TThe matrix was shrunk by a factor alpha=',alpha,' for correction', sep=''))
  }
  
  Results = mapply(function(x){
    R.random = Correlation.matrix.PC(R)
    IND =  ks::vec(matrix(1:d,ncol=2),byrow = TRUE) 
    R.random = R.random[IND,IND]
    ICA =MutivarICA.fun(R.random,Vars,N)
    if(Show.Progress==TRUE){cat((x/M) * 100, "% done... ", sep = "")}
    return(c(ICA,R.random[lower.tri(R.random)]))
  },x=1:M)
  R2_H = Results[1,]
  Corr.R2_H = Results[2,]
  Outcome = list(R2_H=R2_H,Corr.R2_H=Corr.R2_H,Lower.Dig.Corrs.All=t(Results[-c(1:2),]))
  class(Outcome) <- "ICA.ContCont.MultS"
  return(Outcome)
}
