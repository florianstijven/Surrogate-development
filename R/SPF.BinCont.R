SPF.BinCont <- function(x, a, b){ 
  
  R2H <- x$R2_H

  mu_00 <- (x$mu_1_00 - x$mu_0_00)
  mu_01 <- (x$mu_1_01 - x$mu_0_01)
  mu_10 <- (x$mu_1_10 - x$mu_0_10)
  mu_11 <- (x$mu_1_11 - x$mu_0_11)

  sigma_00 <- (x$sigma_00) + (x$sigma_11) - (2 * sqrt((x$sigma_00) * (x$sigma_11)) * x$rho_01_00)
  sigma_01 <- (x$sigma_00) + (x$sigma_11) - (2 * sqrt((x$sigma_00) * (x$sigma_11)) * x$rho_01_01)
  sigma_10 <- (x$sigma_00) + (x$sigma_11) - (2 * sqrt((x$sigma_00) * (x$sigma_11)) * x$rho_01_10)
  sigma_11 <- (x$sigma_00) + (x$sigma_11) - (2 * sqrt((x$sigma_00) * (x$sigma_11)) * x$rho_01_11)

  pi_00 <- x$pi_00
  pi_01 <- x$pi_01
  pi_10 <- x$pi_10
  pi_11 <- x$pi_11
  
  myfunc <- function(a,b,c){
    bp <- max.col(data.frame(a,b,c))
    bp.transform <- ifelse(bp == 1, -1, ifelse(bp == 2, 0, 1))
    return(bp.transform)
  }
  
  # delta S in [-inf, a] = region [-1]
  f_DS_min1_DT_00 <- pi_00*(pnorm(c(a), mean = mu_00, sd = sqrt(sigma_00)))            # f(DS|DT)xf(DT)
  f_DS_min1_DT_01 <- pi_01*(pnorm(c(a), mean = mu_01, sd = sqrt(sigma_01)))
  f_DS_min1_DT_10 <- pi_10*(pnorm(c(a), mean = mu_10, sd = sqrt(sigma_10)))
  f_DS_min1_DT_11 <- pi_11*(pnorm(c(a), mean = mu_11, sd = sqrt(sigma_11)))
  f_DS_min1 <- f_DS_min1_DT_00 + f_DS_min1_DT_01 + f_DS_min1_DT_10 + f_DS_min1_DT_11   # F(DS)
  
  r_min1_min1 <- f_DS_min1_DT_10/f_DS_min1                    # P[DT=-1|DS in (-inf,a)]
  r_0_min1 <- (f_DS_min1_DT_00+f_DS_min1_DT_11)/f_DS_min1     # P[DT=0|DS in (-inf,a)]
  r_1_min1 <- f_DS_min1_DT_01/f_DS_min1                       # P[DT=1|DS in (-inf,a)]
  
  ## DS_min1 <- c(as.numeric(f_DS_min1_DT_10), as.numeric((f_DS_min1_DT_00+f_DS_min1_DT_11)), as.numeric(f_DS_min1_DT_01))
  ## DS_min1_max <- DS_min1[which.max(DS_min1)]
  
  DS_min1_max <- pmax(f_DS_min1_DT_10, (f_DS_min1_DT_00+f_DS_min1_DT_11), f_DS_min1_DT_01)
  best.pred.min1 <- myfunc(r_min1_min1, r_0_min1, r_1_min1)
  
  # delta S in [a,b] = region [0]
  f_DS_0_DT_00 <- pi_00*(pnorm(c(b), mean = mu_00, sd = sqrt(sigma_00)) - 
                             pnorm(c(a), mean = mu_00, sd = sqrt(sigma_00)))
  f_DS_0_DT_01 <- pi_01*(pnorm(c(b), mean = mu_01, sd = sqrt(sigma_01)) - 
                             pnorm(c(a), mean = mu_01, sd = sqrt(sigma_01)))
  f_DS_0_DT_10 <- pi_10*(pnorm(c(b), mean = mu_10, sd = sqrt(sigma_10)) - 
                             pnorm(c(a), mean = mu_10, sd = sqrt(sigma_10)))
  f_DS_0_DT_11 <- pi_11*(pnorm(c(b), mean = mu_11, sd = sqrt(sigma_11)) - 
                             pnorm(c(a), mean = mu_11, sd = sqrt(sigma_11)))
  f_DS_0 <- f_DS_0_DT_00 + f_DS_0_DT_01 + f_DS_0_DT_10 + f_DS_0_DT_11
  
  r_min1_0 <- f_DS_0_DT_10/f_DS_0               # P[DT=-1|DS in (a,b)]
  r_0_0 <- (f_DS_0_DT_00+f_DS_0_DT_11)/f_DS_0   # P[DT=0|DS in (a,b)]
  r_1_0 <- f_DS_0_DT_01/f_DS_0                  # P[DT=1|DS in (a,b)]
  
  ## DS_0 <- c(as.numeric(f_DS_0_DT_10), as.numeric((f_DS_0_DT_00+f_DS_0_DT_11)), as.numeric(f_DS_0_DT_01))
  ## DS_0_max <- DS_0[which.max(DS_0)]
  
  DS_0_max <- pmax(f_DS_0_DT_10, (f_DS_0_DT_00+f_DS_0_DT_11), f_DS_0_DT_01)
  best.pred.0 <- myfunc(r_min1_0, r_0_0, r_1_0)
  
  # delta S in [b, inf] = region [1]
  f_DS_1_DT_00 <- pi_00*(1-pnorm(c(b), mean = mu_00, sd = sqrt(sigma_00)))
  f_DS_1_DT_01 <- pi_01*(1-pnorm(c(b), mean = mu_01, sd = sqrt(sigma_01)))
  f_DS_1_DT_10 <- pi_10*(1-pnorm(c(b), mean = mu_10, sd = sqrt(sigma_10)))
  f_DS_1_DT_11 <- pi_11*(1-pnorm(c(b), mean = mu_11, sd = sqrt(sigma_11)))
  f_DS_1 <- f_DS_1_DT_00 + f_DS_1_DT_01 + f_DS_1_DT_10 + f_DS_1_DT_11
  
  r_min1_1 <- f_DS_1_DT_10/f_DS_1               # P[DT=-1|DS in (b, inf)]
  r_0_1 <- (f_DS_1_DT_00+f_DS_1_DT_11)/f_DS_1   # P[DT=0|DS in (b, inf)]
  r_1_1 <- f_DS_1_DT_01/f_DS_1                  # P[DT=1|DS in (b, inf)]
  
  ## DS_1 <- c(as.numeric(f_DS_1_DT_10), as.numeric((f_DS_1_DT_00+f_DS_1_DT_11)), as.numeric(f_DS_1_DT_01))
  ## DS_1_max <- DS_1[which.max(DS_1)]
  
  DS_1_max <- pmax(f_DS_1_DT_10, (f_DS_1_DT_00+f_DS_1_DT_11), f_DS_1_DT_01)
  best.pred.1 <- myfunc(r_min1_1, r_0_1, r_1_1)

  # P[DT=psi(DS)]
  P_DT_psi_DS_max <- DS_min1_max + DS_0_max + DS_1_max
  
  # Special case
  
  # P[DT=i|DS in (0, inf)] 
  f_DS_0_inf_DT_00 <- pi_00*(1-pnorm(0, mean = mu_00, sd = sqrt(sigma_00)))
  f_DS_0_inf_DT_01 <- pi_01*(1-pnorm(0, mean = mu_01, sd = sqrt(sigma_01)))
  f_DS_0_inf_DT_10 <- pi_10*(1-pnorm(0, mean = mu_10, sd = sqrt(sigma_10)))
  f_DS_0_inf_DT_11 <- pi_11*(1-pnorm(0, mean = mu_11, sd = sqrt(sigma_11)))
  f_DS_0_inf <- f_DS_0_inf_DT_00 + f_DS_0_inf_DT_01 + f_DS_0_inf_DT_10 + f_DS_0_inf_DT_11
  
  P_DT_min1_DS_0_inf <- f_DS_0_inf_DT_10/f_DS_0_inf   # false positive if smaller value of S is harmful
  P_DT_0_DS_0_inf <- (f_DS_0_inf_DT_00+f_DS_0_inf_DT_11)/f_DS_0_inf
  P_DT_1_DS_0_inf <- f_DS_0_inf_DT_01/f_DS_0_inf      # false negative if larger value of S is harmful
  
  # P[DT=i|DS in (-inf, 0)] 
  f_DS_mininf_0_DT_00 <- pi_00*(pnorm(0, mean = mu_00, sd = sqrt(sigma_00)))
  f_DS_mininf_0_DT_01 <- pi_01*(pnorm(0, mean = mu_01, sd = sqrt(sigma_01)))
  f_DS_mininf_0_DT_10 <- pi_10*(pnorm(0, mean = mu_10, sd = sqrt(sigma_10)))
  f_DS_mininf_0_DT_11 <- pi_11*(pnorm(0, mean = mu_11, sd = sqrt(sigma_11)))
  f_DS_mininf_0 <- f_DS_mininf_0_DT_00 + f_DS_mininf_0_DT_01 + f_DS_mininf_0_DT_10 + f_DS_mininf_0_DT_11
  
  P_DT_min1_DS_mininf_0 <- f_DS_mininf_0_DT_10/f_DS_mininf_0   # false positive if larger value of S is harmful
  P_DT_0_DS_mininf_0 <- (f_DS_mininf_0_DT_00+f_DS_mininf_0_DT_11)/f_DS_mininf_0
  P_DT_1_DS_mininf_0 <- f_DS_mininf_0_DT_01/f_DS_mininf_0      # false negative if smaller value of S is harmful
  
  # P[DT=0|DS=0]
  comp_00 <- (pi_00/sqrt(sigma_00))*dnorm((mu_00/sqrt(sigma_00)))
  comp_01 <- (pi_01/sqrt(sigma_01))*dnorm((mu_01/sqrt(sigma_01)))
  comp_10 <- (pi_10/sqrt(sigma_10))*dnorm((mu_10/sqrt(sigma_10)))
  comp_11 <- (pi_11/sqrt(sigma_11))*dnorm((mu_11/sqrt(sigma_11)))
  all_comp <- comp_00 + comp_01 + comp_10 + comp_11
  
  P_DT_0_DS_0 <- (comp_00+comp_11)/all_comp
  
  fit <- list(R2H=as.numeric(R2H), 
              a=a, b=b,
              pi_00=as.numeric(pi_00), 
              pi_01=as.numeric(pi_01), 
              pi_10=as.numeric(pi_10), 
              pi_11=as.numeric(pi_11),
              mu_.00=mu_00, 
              mu_.01=mu_01, 
              mu_.10=mu_10, 
              mu_.11=mu_11,
              sigma_.00=sigma_00, 
              sigma_.01=sigma_01, 
              sigma_.10=sigma_10, 
              sigma_.11=sigma_11,
              f_DS_min1_DT_00=f_DS_min1_DT_00,
              f_DS_min1_DT_01=f_DS_min1_DT_01,
              f_DS_min1_DT_10=f_DS_min1_DT_10,
              f_DS_min1_DT_11=f_DS_min1_DT_11,
              f_DS_min1=f_DS_min1,
              DS_min1_max=DS_min1_max,
              r_min1_min1=r_min1_min1,
              r_0_min1=r_0_min1,
              r_1_min1=r_1_min1,
              best.pred.min1=best.pred.min1,
              f_DS_0_DT_00=f_DS_0_DT_00,
              f_DS_0_DT_01=f_DS_0_DT_01,
              f_DS_0_DT_10=f_DS_0_DT_10,
              f_DS_0_DT_11=f_DS_0_DT_11,
              f_DS_0=f_DS_0,
              DS_0_max=DS_0_max,
              r_min1_0=r_min1_0,
              r_0_0=r_0_0,
              r_1_0=r_1_0,
              best.pred.0=best.pred.0,
              f_DS_1_DT_00=f_DS_1_DT_00,
              f_DS_1_DT_01=f_DS_1_DT_01,
              f_DS_1_DT_10=f_DS_1_DT_10,
              f_DS_1_DT_11=f_DS_1_DT_11,
              f_DS_1=f_DS_1,
              DS_1_max=DS_1_max,
              r_min1_1=r_min1_1,
              r_0_1=r_0_1,
              r_1_1=r_1_1,
              best.pred.1=best.pred.1,
              P_DT_psi_DS_max=P_DT_psi_DS_max,
              f_DS_0_inf=f_DS_0_inf,
              P_DT_min1_DS_0_inf=P_DT_min1_DS_0_inf,
              P_DT_0_DS_0_inf=P_DT_0_DS_0_inf,
              P_DT_1_DS_0_inf=P_DT_1_DS_0_inf,
              f_DS_mininf_0=f_DS_mininf_0,
              P_DT_min1_DS_mininf_0=P_DT_min1_DS_mininf_0,
              P_DT_0_DS_mininf_0=P_DT_0_DS_mininf_0,
              P_DT_1_DS_mininf_0=P_DT_1_DS_mininf_0,
              P_DT_0_DS_0=P_DT_0_DS_0,
              Call=match.call())      

  class(fit) <- "SPF.BinCont"
  fit
}



