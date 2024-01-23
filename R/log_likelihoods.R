clayton_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy, sum_observations){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  #last value in para is the association parameter
  theta <- para[2*(k + 2) + 1]

  #survival probabilities
  u = flexsurv::psurvspline(q = X, gamma = gammax, knots = knotsx, lower.tail = FALSE)
  v = flexsurv::psurvspline(q = Y, gamma = gammay, knots = knotsy, lower.tail = FALSE)
  #densities
  du = flexsurv::dsurvspline(x = X, gamma = gammax, knots = knotsx)
  dv = flexsurv::dsurvspline(x = Y, gamma = gammay, knots = knotsy)

  C <- (u^(-theta) + v^(-theta) - 1)^(-1/theta) #copula

  part1 <- ifelse(d1*d2==1,log(1+theta)+(1+2*theta)*log(C) - (theta+1)*log(u)-(theta+1)*log(v) + log(du) + log(dv),0) #both events
  part2 <- ifelse(d1*(1-d2)==1, (theta+1)*log(C) - (theta+1)*log(u) + log(du),0) #non-terminal event only
  part3 <- ifelse((1-d1)*d2==1, (theta+1)*log(C) - (theta+1)*log(v) + log(dv),0) #terminal event only
  part4 <- ifelse((1-d1)*(1-d2)==1,log(C),0) #both events censored

  loglik = part1+part2+part3+part4
  if (sum_observations){
    loglik <- sum(loglik)
  }

  return(loglik)
}

frank_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy, sum_observations){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  #last value in para is the association parameter
  theta <- para[2*(k + 2) + 1]

  #survival probabilities
  u = flexsurv::psurvspline(q = X, gamma = gammax, knots = knotsx, lower.tail = FALSE)
  v = flexsurv::psurvspline(q = Y, gamma = gammay, knots = knotsy, lower.tail = FALSE)
  #densities
  du = flexsurv::dsurvspline(x = X, gamma = gammax, knots = knotsx)
  dv = flexsurv::dsurvspline(x = Y, gamma = gammay, knots = knotsy)

  #copula
  C <- (-1/theta)*log(((1-exp(-theta)-(1-exp(-theta*u))*(1-exp(-theta*v))))/(1-exp(-theta)))

  part1 <- ifelse(d1*d2==1,(log(theta)+theta*C+log(exp(theta*C)-1)-log(exp(theta*u)-1)-log(exp(theta*v)-1)+log(du)+log(dv)),0)
  part2 <- ifelse(d1*(1-d2)==1,log((1-exp(theta*C))/(1-exp(theta*u)))+log(du),0)
  part3 <- ifelse(((1-d1)*(d2))==1,(log((1-exp(theta*C))/(1-exp(theta*v)))+log(dv)),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)

  loglik = part1+part2+part3+part4
  if (sum_observations){
    loglik <- sum(loglik)
  }

  return(loglik)
}

gumbel_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy, sum_observations){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  #last value in para is the association parameter
  theta <- para[2*(k + 2) + 1]

  #survival probabilities
  u = flexsurv::psurvspline(q = X, gamma = gammax, knots = knotsx, lower.tail = FALSE)
  v = flexsurv::psurvspline(q = Y, gamma = gammay, knots = knotsy, lower.tail = FALSE)
  #densities
  du = flexsurv::dsurvspline(x = X, gamma = gammax, knots = knotsx)
  dv = flexsurv::dsurvspline(x = Y, gamma = gammay, knots = knotsy)

  #copula
  C <- exp(-((-log(u))^(theta)+(-log(v))^(theta))^(1/theta))

  part1 <- ifelse(d1*d2==1,log(C)+(theta-1)*log(-log(u))+(theta-1)*log(-log(v))+
                    log(theta-1-log(C))-log(u)-log(v)-(2*theta-1)*log(-log(C))+
                    log(dv)+log(du),0)
  part2 <- ifelse(d1*(1-d2)==1,log(C)+(theta-1)*log(-log(u))-log(u)-(theta-1)*log(-log(C))+log(du),0)
  part3 <- ifelse(((1-d1)*(d2))==1,(log(C)+(theta-1)*log(-log(v))-log(v)-(theta-1)*log(-log(C))+log(dv)),0)
  part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)

  loglik = part1+part2+part3+part4
  if (sum_observations){
    loglik <- sum(loglik)
  }

  return(loglik)
}

#' @importFrom stats pnorm qnorm
normal_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy, sum_observations){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  # Last value in para is rho. We use the inverse Fisher transformation.
  rho <- (exp(2 * para[2*(k + 2) + 1]) - 1)/(exp(2 * para[2*(k + 2) + 1]) + 1)

  df.1 <- d1 & d2     #case part 1
  df.2 <- d1 & (!d2)  #case part 2
  df.3 <- (!d1)&d2;   #case part 3
  df.4 <- (!d1)&(!d2) #case part 4

  df = data.frame(X, Y)

  X.1 <- df[df.1,1]
  Y.1 <- df[df.1,2]
  X.2 <- df[df.2,1]
  Y.2 <- df[df.2,2]
  X.3 <- df[df.3,1]
  Y.3 <- df[df.3,2]
  X.4 <- df[df.4,1]
  Y.4 <- df[df.4,2]

  part1 <-
    ifelse((sum(df.1) > 0), sum(
      -0.5 * log(1 - rho ^ 2) + (((
        2 * rho * qnorm(flexsurv::psurvspline(
          q = X.1, gamma = gammax, knots = knotsx
        )) * qnorm(flexsurv::psurvspline(
          q = Y.1, gamma = gammay, knots = knotsy
        )) -
          rho ^ 2 * (qnorm(
            flexsurv::psurvspline(q = X.1, gamma = gammax, knots = knotsx)
          ) ^ 2 + qnorm(
            flexsurv::psurvspline(q = Y.1, gamma = gammay, knots = knotsy)
          ) ^ 2)
      )) / ((2 * (
        1 - rho ^ 2
      ))))
      + log(flexsurv::dsurvspline(
        x = X.1, gamma = gammax, knots = knotsx
      )) + log(flexsurv::dsurvspline(
        x = Y.1, gamma = gammay, knots = knotsy
      ))
    ), 0)

  part2 <-
    ifelse((sum(df.2) > 0), sum(log(
      pnorm(
        qnorm(
          flexsurv::psurvspline(
            q = Y.2,
            gamma = gammay,
            knots = knotsy,
            timescale = "log"
          )
        ),
        mean = rho * qnorm(
          flexsurv::psurvspline(
            q = X.2,
            gamma = gammax,
            knots = knotsx,
            timescale = "log"
          )
        ),
        sd = sqrt(1 - rho ^ 2),
        lower.tail = F
      )
    ) + log(
      flexsurv::dsurvspline(
        x = X.2,
        gamma = gammax,
        knots = knotsx,
        timescale = "log"
      )
    )), 0)

  part3 <-
    ifelse((sum(df.3) > 0), sum(log(
      pnorm(
        qnorm(
          flexsurv::psurvspline(
            q = X.3,
            gamma = gammax,
            knots = knotsx,
            timescale = "log"
          )
        ),
        mean = rho * qnorm(
          flexsurv::psurvspline(
            q = Y.3,
            gamma = gammay,
            knots = knotsy,
            timescale = "log"
          )
        ),
        sd = sqrt(1 - rho ^ 2),
        lower.tail = F
      )
    ) + log(
      flexsurv::dsurvspline(
        x = Y.3,
        gamma = gammay,
        knots = knotsy,
        timescale = "log"
      )
    )), 0)

  cov_matrix <- matrix(c(1,rho,rho,1),nrow=2)
  normal_cdf <- function(V,sigma){
    return(mvtnorm::pmvnorm(lower = V,upper=Inf, sigma=sigma,mean=c(0,0))[1])
  }

  part4 <- ifelse((sum(df.4) > 0), sum(log(apply(
    qnorm(cbind(
      flexsurv::psurvspline(
        q = X.4,
        gamma = gammax,
        knots = knotsx,
        timescale = "log"
      ),
      flexsurv::psurvspline(
        q = Y.4,
        gamma = gammay,
        knots = knotsy,
        timescale = "log"
      )
    )), 1, normal_cdf, cov_matrix
  ))), 0)

  loglik = part1+part2+part3+part4
  if (sum_observations){
    loglik <- sum(loglik)
  }

  return(loglik)
}

