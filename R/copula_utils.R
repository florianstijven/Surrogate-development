#' Loglikelihood on the Copula Scale for the Clayton Copula
#'
#' `clayton_loglik_copula_scale()` computes the loglikelihood on the copula
#' scale for the Clayton copula which is parameterized by `theta` as follows:
#' \deqn{C(u, v) = (u^{-\theta} + v^{-\theta} - 1)^{-\frac{1}{\theta}}}
#'
#' @param theta Copula parameter
#' @param u A numeric vector. Corresponds to first variable on the copula scale.
#' @param v A numeric vector. Corresponds to second variable on the copula
#'   scale.
#' @param d1 An integer vector. Indicates whether first variable is observed or
#'   right-censored,
#'  * `d1[i] = 1` if `u[i]` corresponds to non-censored value
#'  * `d1[i] = 0` if `u[i]` corresponds to right-censored value
#' @param d2 An integer vector. Indicates whether first variable is observed or
#'   right-censored,
#'  * `d1[i] = 1` if `u[i]` corresponds to non-censored value
#'  * `d1[i] = 0` if `u[i]` corresponds to right-censored value
#'
#' @return Value of the copula loglikelihood evaluated in `theta`.
clayton_loglik_copula_scale <- function(theta, u, v, d1, d2){
  # Natural logarithm of copula evaluated in u and v.
  log_C = - (1 / theta) * log(u**(-theta) + v**(-theta) - 1)
  # Log likelihood contribution for uncensored observations.
  part1 <-
    ifelse(
      d1 * d2 == 1,
      log(theta + 1) + (2 * theta + 1) * log_C - (theta + 1) * u * v,
      0
    )
  # Log likelihood contribution for second observation censored.
  part2 <-
    ifelse(d1 * (1 - d2) == 1,
           (theta + 1) * log_C - (theta + 1) * u,
           0)
  # Log likelihood contribution for first observation censored.
  part2 <-
    ifelse(d1 * (1 - d2) == 1,
           (theta + 1) * log_C - (theta + 1) * v,
           0)
  # Log likelihood contribution for both observations censored.
  part4 <- ifelse((1-d1)*(1-d2)==1,log_C,0)

  loglik_copula <- part1 + part2 + part3 + part4

  return(loglik_copula)
}

#' Loglikelihood on the Copula Scale for the Frank Copula
#'
#' `frank_loglik_copula_scale()` computes the loglikelihood on the copula
#' scale for the Frank copula which is parameterized by `theta` as follows:
#' \deqn{ C(u, v) = - \frac{1}{\theta} \log \left[ 1 - \frac{(1 - e^{-\theta u})(1 - e^{-\theta v})}{1 - e^{-\theta}} \right]}
#'
#' @inheritParams clayton_loglik_copula_scale
#'
#' @return Value of the copula loglikelihood evaluated in `theta`.
frank_loglik_copula_scale <- function(theta, u, v, d1, d2){

  # For efficiency purposes, some quantities that are needed multiple times are
  # first precomputed. In this way, we do not waste resources on computing the
  # same quantity multiple times.

  A_u = 1 - exp(-theta * u)
  A_v = 1 - exp(-theta * v)
  A_theta = 1 - exp(-theta)

  # Copula evaluated in (u, v)
  C = (-1 / theta) * log(1 - ((A_u * A_v) / A_theta))


  # Log likelihood contribution for uncensored observations.
  part1 <-
    ifelse(d1 * d2 == 1,
           log(theta) + theta * C + log(exp(theta * C) - 1) - log(A_u) - log(A_v),
           0)
  # Log likelihood contribution for second observation censored.
  part2 <-
    ifelse(d1 * (1 - d2) == 1,
           log(1 - exp(theta * C) - log(A_u)),
           0)
  # Log likelihood contribution for first observation censored.
  part3 <-
    ifelse(d1 * (1 - d2) == 1,
           log(1 - exp(theta * C) - log(A_v)),
           0)
  # Log likelihood contribution for both observations censored.
  part4 <- ifelse((1 - d1) * (1 - d2) == 1, log(C), 0)

  loglik_copula <- part1 + part2 + part3 + part4

  return(loglik_copula)
}

#' Loglikelihood on the Copula Scale for the Gumbel Copula
#'
#' `gumbel_loglik_copula_scale()` computes the loglikelihood on the copula
#' scale for the Gumbel copula which is parameterized by `theta` as follows:
#' \deqn{C(u, v) = \exp \left[ - \left\{ (-\log u)^{\theta} + (-\log v)^{\theta} \right\}^{\frac{1}{\theta}} \right]}
#'
#' @inheritParams clayton_loglik_copula_scale
#' @return Value of the copula loglikelihood evaluated in `theta`.
gumbel_loglik_copula_scale <- function(theta, u, v, d1, d2){
  # For efficiency purposes, some quantities that are needed multiple times are
  # first precomputed. In this way, we do not waste resources on computing the
  # same quantity multiple times.

  min_log_u = log(u)
  min_log_v = log(v)
  # Log Copula evaluated in (u, v)
  log_C = -(min_log_u**theta + min_log_v**theta)**(1 / theta)

  # Log likelihood contribution for uncensored observations.
  part1 <-
    ifelse(
      d1 * d2 == 1,
      log_C + (theta - 1) * (log(min_log_u) + log(min_log_v)) +
        log(theta - 1 - log_C) + min_log_u + min_log_v - (2 * theta - 1) * log(-log_C),
      0
    )
  # Log likelihood contribution for second observation censored.
  part2 <-
    ifelse(d1 * (1 - d2) == 1,
           log_C + (theta - 1) * log(min_log_u) + min_log_u - (theta - 1) * log(-log_C),
           0)
  # Log likelihood contribution for first observation censored.
  part2 <-
    ifelse(d1 * (1 - d2) == 1,
           log_C + (theta - 1) * log(min_log_v) + min_log_v - (theta - 1) * log(-log_C),
           0)
  # Log likelihood contribution for both observations censored.
  part4 <- ifelse((1 - d1) * (1 - d2) == 1, log_C, 0)

  loglik_copula <- part1 + part2 + part3 + part4

  return(loglik_copula)
}

#' @importFrom stats pnorm qnorm
normal_loglik <- function(para, X, Y, d1, d2, k = 2, knotsx, knotsy){
  #k is the number of knots in the model, this determines the length of para
  gammax <- para[1:(k + 2)]
  gammay <- para[(k + 3):(2*(k + 2))]
  #last value in para is rho
  rho <- (exp(para[2*(k + 2) + 1]) - 1)/(exp(para[2*(k + 2) + 1]) + 1)

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

  loglik <- (part1+part2+part3+part4)
  return(loglik)
}

