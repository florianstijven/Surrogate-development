library(testthat)
library(pracma)

ICA_t <- function(df, T0S0, T1S1, T0T0=1, T1T1=1, S0S0=1, S1S1=1,
                  T0T1=seq(-1, 1, by=.1), T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1),
                  S0S1=seq(-1, 1, by=.1)) {

  T0S0_hier <- T0S0[1]
  T1S1_hier <- T1S1[1]

  Results <- na.exclude(matrix(NA, 1, 11))
  colnames(Results) <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "rho", "ICA", "ICA_t", "sigma.delta.T", "sigma.delta.S")
  combins <- expand.grid(T0T1, T0S0_hier, T0S1, T1S0, T1S1_hier, S0S1)
  lengte <- dim(combins)[1]
  zeta_fun <- function(nu) 2*log( (gamma(nu/2)*gamma(1/2))/(sqrt(pi)*gamma( (1+nu)/2 ))*sqrt(nu/2))- (nu/2)*psi(nu/2)+(nu+1)*psi((1+nu)/2)- ((nu+2)/2)*psi((nu+2)/2)

  if (length(T0S0)>1){
    if (length(T0S0)<lengte){stop("The specified vector for T0S0 should be larger than ", lengte) }
    T0S0_vector <- T0S0[1:lengte] # sample
    combins[,2] <- T0S0_vector
  }
  if (length(T1S1)>1){
    if (length(T1S1)<lengte){stop("The specified vector for T1S1 should be larger than ", lengte) }
    T1S1_vector <- T1S1[1:lengte] # sample
    combins[,5] <- T1S1_vector
  }


  for (i in 1: nrow(combins)) {
    T0T1 <- combins[i, 1]
    T0S0 <- combins[i, 2]
    T0S1 <- combins[i, 3]
    T1S0 <- combins[i, 4]
    T1S1 <- combins[i, 5]
    S0S1 <- combins[i, 6]
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

    if (Min.Eigen.Cor > 0) {
      rho <- ((sqrt(S0S0*T0T0)*Cor_c[3,1])+(sqrt(S1S1*T1T1)*Cor_c[4,2])-(sqrt(S0S0*T1T1)*Cor_c[3,2])-(sqrt(S1S1*T0T0)*Cor_c[4,1]))/(sqrt((T0T0+T1T1-(2*sqrt(T0T0*T1T1)*Cor_c[2,1]))*(S0S0+S1S1-(2*sqrt(S0S0*S1S1)*Cor_c[4,3]))))
      ICA<- rho^2
      if ((is.finite(ICA))==TRUE){
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * Cor_c[2,1])
        sigma.delta.S <- S0S0 + S1S1 - (2 * sqrt(S0S0*S1S1) * Cor_c[3,4])
        ICA_t<- 1-(1-rho^2)*exp(-2*zeta_fun(df))
        results.part <- as.vector(cbind(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1, rho, ICA, ICA_t, sigma.delta.T, sigma.delta.S))
        Results <- rbind(Results, results.part)
        if( (df > 343)==TRUE)  { print(" The maximum value of df is 342. When it is greater than 342, it is the same as ICA in the normal causal model. ")  }
        rownames(Results) <- NULL}
    }
  }
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  Total.Num.Matrices <- nrow(combins)

  fit <-
    list(Total.Num.Matrices=Total.Num.Matrices, Pos.Def=Results[,1:6], rho=Results$rho, ICA=Results$ICA, ICA_t=Results$ICA_t, Sigmas=Results[,10:11], Call=match.call())

  class(fit) <- "ICA_t"
  fit

}



##test1
test_that("ICA_t", {
  df = 3
  T0S0 = 0.9597334
  T1S1=0.9644139
  T0T0=544.3285
  T1T1=550.6597
  S0S0=180.6831
  S1S1=180.9433
  T0T1=-0.9
  T0S1=-0.9
  T1S0=-0.9
  S0S1=-0.9
  log_lik = ICA_t(df=3, T0S0 = 0.9597334, T1S1=0.9644139, T0T0=544.3285, T1T1=550.6597, S0S0=180.6831, S1S1=180.9433, T0T1=-0.9,
                  T0S1=-0.9 , T1S0=-0.9 , S0S1=-0.9)
  expect_equal(round(log_lik$ICA_t, 3), 0.964)
})


##test2
test_that("ICA_t", {
  df = 500
  T0S0 = 0.9597334
  T1S1=0.9644139
  T0T0=544.3285
  T1T1=550.6597
  S0S0=180.6831
  S1S1=180.9433
  T0T1=-0.9
  T0S1=-0.9
  T1S0=-0.9
  S0S1=-0.9
  log_lik = ICA_t(df=500, T0S0 = 0.9597334, T1S1=0.9644139, T0T0=544.3285, T1T1=550.6597, S0S0=180.6831, S1S1=180.9433, T0T1=-0.9,
                  T0S1=-0.9 , T1S0=-0.9 , S0S1=-0.9)
  expect_equal(log_lik$ICA_t, NaN)
})

