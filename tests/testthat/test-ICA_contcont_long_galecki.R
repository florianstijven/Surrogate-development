library(testthat)
library(pracma)

ICA_contcont_long_galecki <- function(p=numeric(), T0S0, T1S1, T0T0, T1T1, S0S0, S1S1,
                                      T0T1=seq(-1, 1, by=.1), T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1),
                                      S0S1=seq(-1, 1, by=.1)) {


  Sigma_c <- matrix(NA,nrow=4, ncol=4)

  # Fill covariances
  Sigma_c[4,2] <- Sigma_c[2,4] <- T1S1
  Sigma_c[3,1] <- Sigma_c[1,3] <- T0S0

  Sigma_c[1,1] <-  T0T0
  Sigma_c[2,2] <-  T1T1
  Sigma_c[3,3] <-  S0S0
  Sigma_c[4,4] <-  S1S1

  Cor_c <- cov2cor(Sigma_c)

  T0S0 <- Cor_c[1,3]
  T1S1 <- Cor_c[2,4]
  T0T0 <- 1
  T1T1 <- 1
  S0S0 <- 1
  S1S1 <- 1

  T0S0_hier <- T0S0[1]
  T1S1_hier <- T1S1[1]
  p <- as.numeric(p)
  Results <- na.exclude(matrix(NA, 1, 11))
  colnames_result <- c("T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1", "sigma.delta.T", "sigma.delta.S", "rho_u", "rho_u2", "R2_Lambda")
  colnames(Results) <- c(colnames_result)

  combins <- expand.grid(T0T1, T0S0_hier, T0S1, T1S0, T1S1_hier, S0S1)
  lengte <- dim(combins)[1]

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
      rho_u <- ((sqrt(S0S0*T0T0)*Cor_c[3,1])+(sqrt(S1S1*T1T1)*Cor_c[4,2])-(sqrt(S0S0*T1T1)*Cor_c[3,2])-(sqrt(S1S1*T0T0)*Cor_c[4,1]))/(sqrt((T0T0+T1T1-(2*sqrt(T0T0*T1T1)*Cor_c[2,1]))*(S0S0+S1S1-(2*sqrt(S0S0*S1S1)*Cor_c[4,3]))))
      rho_u2<- rho_u^2
      if ((is.finite(rho_u2))==TRUE){
        sigma.delta.T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * Cor_c[2,1])
        sigma.delta.S <- S0S0 + S1S1 - (2 * sqrt(S0S0*S1S1) * Cor_c[3,4])

        R2_Lambda<-1 - (1-rho_u2)^p
        R2_Lambda<- as.numeric(R2_Lambda)
        results.part_1 <- as.vector(cbind(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1,  sigma.delta.T, sigma.delta.S, rho_u, rho_u2, R2_Lambda))

        if (!exists("Results")) {
          Results <- data.frame(matrix(ncol = length(results.part_1), nrow = 0))
        }

        Results <- rbind(Results, results.part_1)
      }

      # Results <- rbind(Results, results.part)
      # rownames(Results) <- NULL}
    }
  }
  Results <- data.frame(Results)
  rownames(Results) <- NULL
  Total.Num.Matrices <- nrow(combins)

  fit <-
    list(Total.Num.Matrices=Total.Num.Matrices, Pos.Def=Results[,1:6], rho_u=Results$rho_u, Sigmas=Results[,7:8], rho_u2=Results$rho_u2, R2_Lambda=Results$R2_Lambda, Call=match.call())

  class(fit) <- "ICA_contcont_long_galecki"
  fit

}

test_that("ICA_contcont_long_galecki",
         {p=4
         T0S0=266.56
         T1S1=205.17
         T0T0=544.3285
         T1T1=550.6597
         S0S0=180.6831
         S1S1=180.9433
         T0T1=0.9
         T0S1=0.9
         T1S0=0.9
         S0S1=0.9
         loglik= ICA_contcont_long_galecki(p=4,  T0S0 = 266.56, T1S1=205.17,
                                                 T0T0=544.3285, T1T1=550.6597,
                                                 S0S0=180.6831, S1S1=180.9433,
                                                 T0T1= -0.8,  T0S1= -0.8,
                                                 T1S0= -0.8 , S0S1= -0.8)
         expect_equal(loglik$R2_Lambda, 0.9955342, tolerance = 1e-6)})
