SPF.BinCont <- function(x, a, b){ 

mu_star_00 <- (x$mu_1_00 - x$mu_0_00)
mu_star_01 <- (x$mu_1_01 - x$mu_0_01)
mu_star_10 <- (x$mu_1_10 - x$mu_0_10)
mu_star_11 <- (x$mu_1_11 - x$mu_0_11)

sigma_star_00 <- (x$sigma2_00_00) + (x$sigma2_11_00) - 
  (2 * sqrt((x$sigma2_00_00) * (x$sigma2_11_00)) * x$G_rho_01_00)
sigma_star_01 <- (x$sigma2_00_01) + (x$sigma2_11_01) - 
  (2 * sqrt((x$sigma2_00_01) * (x$sigma2_11_01)) * x$G_rho_01_01)
sigma_star_10 <- (x$sigma2_00_10) + (x$sigma2_11_10) - 
  (2 * sqrt((x$sigma2_00_10) * (x$sigma2_11_10)) * x$G_rho_01_10)
sigma_star_11 <- (x$sigma2_00_11) + (x$sigma2_11_11) - 
  (2 * sqrt((x$sigma2_00_11) * (x$sigma2_11_11)) * x$G_rho_01_11)

# N equations
comp1 <- x$pi_00 * (pnorm(c(b), mean=mu_star_00, sd=sqrt(sigma_star_00), lower.tail=TRUE) - 
                      pnorm(c(a), mean=mu_star_00, sd=sqrt(sigma_star_00), lower.tail=TRUE))
comp2 <- x$pi_01 * (pnorm(c(b), mean=mu_star_01, sd=sqrt(sigma_star_01), lower.tail=TRUE) - 
                      pnorm(c(a), mean=mu_star_01, sd=sqrt(sigma_star_01), lower.tail=TRUE))
comp3 <- x$pi_10 * (pnorm(c(b), mean=mu_star_10, sd=sqrt(sigma_star_10), lower.tail=TRUE) - 
                      pnorm(c(a), mean=mu_star_10, sd=sqrt(sigma_star_10), lower.tail=TRUE))
comp4 <- x$pi_11 * (pnorm(c(b), mean=mu_star_11, sd=sqrt(sigma_star_11), lower.tail=TRUE) - 
                      pnorm(c(a), mean=mu_star_11, sd=sqrt(sigma_star_11), lower.tail=TRUE))
N <- comp1 + comp2 + comp3 + comp4

# P(delta T = -1 | ...)
Tel <- x$pi_10 * ((pnorm(c(b), mean=mu_star_10, sd=sqrt(sigma_star_10), lower.tail=TRUE) - 
                    pnorm(c(a), mean=mu_star_10, sd=sqrt(sigma_star_10), lower.tail=TRUE)))
P_Delta_T_min1 <- Tel / N

# P(delta T = 0 | ...)
Tel1 <- x$pi_00 * ((pnorm(c(b), mean=mu_star_00, sd=sqrt(sigma_star_00), lower.tail=TRUE) - 
                  pnorm(c(a), mean=mu_star_00, sd=sqrt(sigma_star_00), lower.tail=TRUE)))
Tel2 <- x$pi_11 * ((pnorm(c(b), mean=mu_star_11, sd=sqrt(sigma_star_11), lower.tail=TRUE) - 
                     pnorm(c(a), mean=mu_star_11, sd=sqrt(sigma_star_11), lower.tail=TRUE)))
P_Delta_T_0 <- (Tel1 + Tel2) / N

# P(delta T = 1 | ...)
Tel <- x$pi_01 * ((pnorm(c(b), mean=mu_star_01, sd=sqrt(sigma_star_01), lower.tail=TRUE) - 
                    pnorm(c(a), mean=mu_star_01, sd=sqrt(sigma_star_01), lower.tail=TRUE)))
P_Delta_T_1 <- Tel / N

# Check
#P_Delta_T_min1 + P_Delta_T_0 + P_Delta_T_1

fit <- 
  list(a=a, b=b, P_Delta_T_min1=P_Delta_T_min1, P_Delta_T_0=P_Delta_T_0, P_Delta_T_1=P_Delta_T_1, Call=match.call())      

class(fit) <- "SPF.BinCont"
fit
}




