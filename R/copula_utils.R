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
#'  * `d1[i] = -1` if `u[i]` corresponds to left-censored value
#' @param d2 An integer vector. Indicates whether first variable is observed or
#'   right-censored,
#'  * `d2[i] = 1` if `v[i]` corresponds to non-censored value
#'  * `d2[i] = 0` if `v[i]` corresponds to right-censored value
#'  * `d2[i] = -1` if `v[i]` corresponds to left-censored value
#' @param return_sum Return the sum of the individual loglikelihoods? If `FALSE`,
#' a vector with the individual loglikelihood contributions is returned.
#'
#' @return Value of the copula loglikelihood evaluated in `theta`.
clayton_loglik_copula_scale <- function(theta, u, v, d1, d2, return_sum = TRUE) {
  # Replace entries that correspond to zero or one probabilities with NAs. This
  # avoids NaN warnings.
  zero_indicator = detect_zero_probs(u, v, d1, d2)
  one_indicator = detect_one_probs(u, v, d1, d2)
  zero_or_one = zero_indicator | one_indicator
  u[zero_or_one] = NA; v[zero_or_one] = NA
  d1[zero_or_one] = NA; d2[zero_or_one] = NA

  # Natural logarithm of copula evaluated in u and v.
  log_C = -(1 / theta) * log(u ** (-theta) + v ** (-theta) - 1)
  # Log likelihood contribution for uncensored observations.
  part1 <-
    ifelse((d1 == 1) & (d2 == 1),
           log(theta + 1) + (2 * theta + 1) * log_C - (theta + 1) * (log(u) + log(v)),
           0)
  # Log likelihood contribution for second observation left-censored.
  part2 <-
    ifelse((d1 == 1) & (d2 == -1),
           (theta + 1) * log_C - (theta + 1) * log(u),
           0)
  # Log likelihood contribution for first observation left-censored.
  part3 <-
    ifelse((d1 == -1) & (d2 == 1),
           (theta + 1) * log_C - (theta + 1) * log(v),
           0)
  # Log likelihood contribution for both observations left-censored.
  part4 <- ifelse((d1 == -1) & (d2 == -1), log_C, 0)

  # Log likelihood contribution for second observation right-censored.
  part5 <- ifelse((d1 == 1) & (d2 == 0),
                  log(1 - (exp(log_C) / u) ** (theta + 1)),
                  0)

  # Log likelihood contribution for first observation right-censored.
  part6 <- ifelse((d1 == 0) & (d2 == 1),
                  log(1 - (exp(log_C) / v) ** (theta + 1)),
                  0)

  # Log likelihood contribution for both observations right-censored.
  part7 <- ifelse((d1 == 0) & (d2 == 0),
                  log(1 + exp(log_C) - u - v),
                  0)

  # Log likelihood contribution for first observation left-censored and second
  # observation right-censored.
  part8 <- ifelse((d1 == -1) & (d2 == 0),
                  log(u - exp(log_C)),
                  0)

  # Log likelihood contribution for first observation right-censored and second
  # observation left-censored.
  part9 <- ifelse((d1 == 0) & (d2 == -1),
                  log(v - exp(log_C)),
                  0)

  loglik_copula <-
    part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8 + part9
  loglik_copula[zero_indicator] = -Inf
  loglik_copula[one_indicator] = 0
  if (return_sum)
    loglik_copula = sum(loglik_copula)

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
frank_loglik_copula_scale <- function(theta, u, v, d1, d2, return_sum = TRUE){
  # Replace entries that correspond to zero or one probabilities with NAs. This
  # avoids NaN warnings.
  zero_indicator = detect_zero_probs(u, v, d1, d2)
  one_indicator = detect_one_probs(u, v, d1, d2)
  zero_or_one = zero_indicator | one_indicator
  u[zero_or_one] = NA; v[zero_or_one] = NA
  d1[zero_or_one] = NA; d2[zero_or_one] = NA

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
    ifelse(
      (d1 == 1) & (d2 == 1),
      log(theta) + theta * C + log(exp(theta * C) - 1) - log((exp(theta * u) - 1) * (exp(theta * v) - 1)),
      0)
  # Log likelihood contribution for second observation left-censored.
  part2 <-
    ifelse((d1 == 1) & (d2 == -1),
           log((1 - exp(theta * C)) / (1 - exp(theta * u))),
           0)
  # Log likelihood contribution for first observation left-censored.
  part3 <-
    ifelse((d1 == -1) & (d2 == 1),
           log((1 - exp(theta * C)) / (1 - exp(theta * v))),
           0)
  # Log likelihood contribution for both observations left-censored.
  part4 <- ifelse((d1 == -1) & (d2 == -1), log(C), 0)
  # Log likelihood contribution for second observation right-censored.
  part5 <- ifelse((d1 == 1) & (d2 == 0),
                  log(1 - ((1 - exp(theta * C)) / (1 - exp(theta * u)))),
                  0)
  # Log likelihood contribution for first observation right-censored.
  part6 <- ifelse((d1 == 0) & (d2 == 1),
                  log(1 - ((1 - exp(theta * C)) / (1 - exp(theta * v)))),
                  0)
  # Log likelihood contribution for both observations right-censored.
  part7 <- ifelse((d1 == 0) & (d2 == 0),
                  log(1 + C - u - v),
                  0)
  # Log likelihood contribution for first observation left-censored and second
  # observation right-censored.
  part8 <- ifelse((d1 == -1) & (d2 == 0),
                  log(u - C),
                  0)
  # Log likelihood contribution for first observation right-censored and second
  # observation left-censored.
  part9 <- ifelse((d1 == 0) & (d2 == -1),
                  log(v - C),
                  0)

  loglik_copula <-
    part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8 + part9
  loglik_copula[zero_indicator] = -Inf
  loglik_copula[one_indicator] = 0
  if (return_sum)
    loglik_copula = sum(loglik_copula)

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
gumbel_loglik_copula_scale <- function(theta, u, v, d1, d2, return_sum = TRUE){
  # Replace entries that correspond to zero or one probabilities with NAs. This
  # avoids NaN warnings.
  zero_indicator = detect_zero_probs(u, v, d1, d2)
  one_indicator = detect_one_probs(u, v, d1, d2)
  zero_or_one = zero_indicator | one_indicator
  u[zero_or_one] = NA; v[zero_or_one] = NA
  d1[zero_or_one] = NA; d2[zero_or_one] = NA

  # For efficiency purposes, some quantities that are needed multiple times are
  # first precomputed. In this way, we do not waste resources on computing the
  # same quantity multiple times.

  min_log_u = -log(u)
  min_log_v = -log(v)
  # Log Copula evaluated in (u, v)
  log_C = -(min_log_u**theta + min_log_v**theta)**(1 / theta)

  # Log likelihood contribution for uncensored observations.
  part1 <-
    ifelse(
      (d1 == 1) & (d2 == 1),
      log_C + (theta - 1) * (log(min_log_u) + log(min_log_v)) +
        log(theta - 1 - log_C) + min_log_u + min_log_v - (2 * theta - 1) * log(-log_C),
      0
    )
  # Log likelihood contribution for second observation left-censored.
  part2 <-
    ifelse((d1 == 1) & (d2 == -1),
           log_C + (theta - 1) * log(min_log_u) + min_log_u - (theta - 1) * log(-log_C),
           0)
  # Log likelihood contribution for first observation left-censored.
  part3 <-
    ifelse((d1 == -1) & (d2 == 1),
           log_C + (theta - 1) * log(min_log_v) + min_log_v - (theta - 1) * log(-log_C),
           0)
  # Log likelihood contribution for both observations left-censored.
  part4 <- ifelse((d1 == -1) & (d2 == -1), log_C, 0)
  # Log likelihood contribution for second observation right-censored.
  part5 <- ifelse((d1 == 1) & (d2 == 0),
                  log(1 - (exp(log_C) / u) * (log(u - exp(log_C)) ** (theta - 1))),
                  0)
  # Log likelihood contribution for first observation right-censored.
  part6 <- ifelse((d1 == 0) & (d2 == 1),
                  log(1 - (exp(log_C) / v) * (log(v - exp(log_C)) ** (theta - 1))),
                  0)
  # Log likelihood contribution for both observations right-censored.
  part7 <- ifelse((d1 == 0) & (d2 == 0),
                  log(1 + exp(log_C) - u - v),
                  0)
  # Log likelihood contribution for first observation left-censored and second
  # observation right-censored.
  part8 <- ifelse((d1 == -1) & (d2 == 0),
                  log(u - exp(log_C)),
                  0)
  # Log likelihood contribution for first observation right-censored and second
  # observation left-censored.
  part9 <- ifelse((d1 == 0) & (d2 == -1),
                  log(v - exp(log_C)),
                  0)

  loglik_copula <-
    part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8 + part9
  loglik_copula[zero_indicator] = -Inf
  loglik_copula[one_indicator] = 0
  if (return_sum)
    loglik_copula = sum(loglik_copula)

  return(loglik_copula)
}


#' Loglikelihood on the Copula Scale for the Gaussian Copula
#'
#' `gaussian_loglik_copula_scale()` computes the loglikelihood on the copula
#' scale for the Gaussian copula which is parameterized by `theta` as follows:
#' \deqn{C(u, v) = \Psi \left[ \Phi^{-1} (u), \Phi^{-1} (v) | \rho \right]}
#'
#' @inheritParams clayton_loglik_copula_scale
#'
#' @return Value of the copula loglikelihood evaluated in `theta`.
gaussian_loglik_copula_scale <- function(theta, u, v, d1, d2, return_sum = TRUE){
  # Replace entries that correspond to zero or one probabilities with defaults.
  # This avoids NaN warnings. We don't use NAs here because the functions from
  # mvtnorm cannot handle NAs.
  zero_indicator = detect_zero_probs(u, v, d1, d2)
  one_indicator = detect_one_probs(u, v, d1, d2)
  zero_or_one = zero_indicator | one_indicator
  u[zero_or_one] = 0.5; v[zero_or_one] = 0.5
  d1[zero_or_one] = 0; d2[zero_or_one] = 0

  requireNamespace("mvtnorm")
  # For efficiency purposes, some quantities that are needed multiple times are
  # first precomputed. In this way, we do not waste resources on computing the
  # same quantity multiple times.

  # Covariance matrix
  Sigma = matrix(c(1, theta, theta, 1), nrow = 2)

  # Transform the variables on the copula scale to the standard normal scale.
  z_u = stats::qnorm(u)
  z_v = stats::qnorm(v)

  # Merge the variables on the standard normal scale into a matrix with one
  # column for each variable.
  z_UV = matrix(c(z_u, z_v), byrow = FALSE, ncol = 2)
  # Evaluate the Gaussian copula in (uv).
  C = apply(
    X = z_UV,
    FUN = function(x)
      mvtnorm::pmvnorm(
        lower = c(-Inf, -Inf),
        upper = x,
        sigma = Sigma,
        keepAttr = FALSE
      ),
    MARGIN = 1,
    simplify = TRUE
  )

  # Log likelihood contribution for uncensored observations.
  part1 <-
    ifelse((d1 == 1) & (d2 == 1),
           mvtnorm::dmvnorm(x = z_UV, sigma = Sigma, log = TRUE) -
             dnorm(z_u, log = TRUE) -
             dnorm(z_v, log = TRUE),
           0
    )
  # Log likelihood contribution for second observation left-censored.
  part2 <-
    ifelse((d1 == 1) & (d2 == -1),
           pnorm(
             q = z_v,
             mean = z_u * theta,
             sd =  sqrt(1 - theta ** 2),
             log.p = TRUE
           ),
           0)
  # Log likelihood contribution for first observation left-censored.
  part3 <-
    ifelse((d1 == -1) & (d2 == 1),
           pnorm(
             q = z_u,
             mean = z_v * theta,
             sd =  sqrt(1 - theta ** 2),
             log.p = TRUE
           ),
           0)
  # Log likelihood contribution for both observations left-censored.
  part4 <-
    ifelse((d1 == -1) & (d2 == -1),
           log(C),
           0)
  # Log likelihood contribution for second observation right-censored.
  part5 <- ifelse((d1 == 1) & (d2 == 0),
                  log(
                    1 - pnorm(
                      q = z_v,
                      mean = z_u * theta,
                      sd =  sqrt(1 - theta ** 2),
                      lower.tail = TRUE,
                      log.p = FALSE
                    )
                  ),
                  0)
  # Log likelihood contribution for first observation right-censored.
  part6 <- ifelse((d1 == 0) & (d2 == 1),
                  log(
                    1 - pnorm(
                      q = z_u,
                      mean = z_v * theta,
                      sd =  sqrt(1 - theta ** 2),
                      lower.tail = TRUE,
                      log.p = FALSE
                    )
                  ),
                  0)
  # Log likelihood contribution for both observations right-censored.
  part7 <- ifelse((d1 == 0) & (d2 == 0),
                  log(1 + C - u - v),
                  0)

  # Log likelihood contribution for first observation left-censored and second
  # observation right-censored.
  part8 <- ifelse((d1 == -1) & (d2 == 0),
                  log(u - C),
                  0)

  # Log likelihood contribution for first observation right-censored and second
  # observation left-censored.
  part9 <- ifelse((d1 == 0) & (d2 == -1),
                  log(v - C),
                  0)

  loglik_copula <-
    part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8 + part9
  loglik_copula[zero_indicator] = -Inf
  loglik_copula[one_indicator] = 0
  if (return_sum)
    loglik_copula = sum(loglik_copula)

  return(loglik_copula)
}


#' Loglikelihood on the Copula Scale
#'
#' `loglik_copula_scale()` computes the loglikelihood on the copula scale for
#' possibly right-censored data.
#'
#' @inheritParams clayton_loglik_copula_scale
#' @param copula_family Copula family, one of the following:
#' * `"clayton"`
#' * `"frank"`
#' * `"gumbel"`
#' * `"gaussian"`
#' @param r rotation parameter. Should be `0L`, `90L`, `180L`, or `270L`.
#'
#'   The parameterization of the respective copula families can be found in the
#'   help files of the dedicated functions named `copula_loglik_copula_scale()`.
#'
#' @return Value of the copula loglikelihood evaluated in `theta`.
loglik_copula_scale <- function(theta, u, v, d1, d2, copula_family, r = 0L, return_sum = TRUE){
  # The specific copula loglikelihood functions have been implemented for the
  # non-rotated case. By changing u, v, d1, and d2, we can compute the
  # corresponding loglikelihoods for the rotated copulas.
  switch(
    r,
    "180" = {
      u = 1 - u
      v = 1 - v
      d1 = -1 * d1
      d2 = -1 * d2
    },
    "90" = {
      u = 1 - u
      d1 = -1 * d1
    },
    "270" = {
      v = 1 - v
      d2 = -1 * d2
    }
  )
  loglik_copula = switch(
    copula_family,
    "clayton" = clayton_loglik_copula_scale(theta, u, v, d1, d2, return_sum = return_sum),
    "frank" = frank_loglik_copula_scale(theta, u, v, d1, d2, return_sum = return_sum),
    "gumbel" = gumbel_loglik_copula_scale(theta, u, v, d1, d2, return_sum = return_sum),
    "gaussian" = gaussian_loglik_copula_scale(theta, u, v, d1, d2, return_sum = return_sum)
  )

  return(loglik_copula)
}


detect_zero_probs = function(u, v, d1, d2) {
  # If U = 1 and there is right-censoring, the probability is zero.
  zero_probs_indicator = ifelse(((u == 1) &
                                   (d1 == 0)) | ((v == 1) & (d2 == 0)), TRUE, FALSE)
  # If U = 0 and there is left-censoring, the probability is zero.
  zero_probs_indicator = zero_probs_indicator |
    ifelse(((u == 0) &
              (d1 == -1)) | ((v == 0) & (d2 == -1)), TRUE, FALSE)
  return(zero_probs_indicator)
}

detect_one_probs = function(u, v, d1, d2) {
  # If U = 0 and V = 0 and there is right-censoring, the probability is one.
  zero_probs_indicator = ifelse(((u == 0) &
                                   (d1 == 0)) & ((v == 0) & (d2 == 0)), TRUE, FALSE)
  # If U = 1 and V = 1 and there is left-censoring, the probability is one.
  zero_probs_indicator = zero_probs_indicator |
    ifelse(((u == 1) &
              (d1 == -1)) & ((v == 1) & (d2 == -1)), TRUE, FALSE)
  return(zero_probs_indicator)
}

marginal_ord_constructor = function(param) {
  cum_probs = c(pnorm(param), 1)
  probs = cum_probs - c(0, cum_probs[-length(cum_probs)])
  K = length(param) + 1
  cdf = function(x) {
    cum_probs[x]
  }
  pmf = function(x) {
    probs[x]
  }
  inv_cdf = function(p) {
    sapply(p, function(p) {
      max((1:K)[c(0, cum_probs) < p])
    })
  }
  list(
    pmf = pmf,
    cdf = cdf,
    inv_cdf = inv_cdf
  )
}

marginal_cont_constructor = function(marginal_Y, param) {
  pdf = function(x) {
    marginal_Y[[1]](x, param)
  }
  cdf = function(x) {
    marginal_Y[[2]](x, param)
  }
  inv_cdf = function(p) {
    marginal_Y[[3]](p, param)
  }
  return(list(
    pdf = pdf,
    cdf = cdf,
    inv_cdf = inv_cdf
  ))
}

#' Estimate marginal distribution using ML
#'
#' @inheritParams fit_copula_submodel_OrdCont_twostep
#'
#' @return Estimated parameters
estimate_marginal = function(Y, marginal_Y, starting_values) {
  # Define log-likelihood function (as a function of the parameters).
  log_lik_function = function(para) {
    sum(log(marginal_Y[[1]](Y, para)))
  }

  coef(maxLik::maxLik(
    logLik = log_lik_function,
    start = starting_values,
    method = "NR"
  ))
}

#' Convert Ordinal Observations to Latent Cutpoints
#'
#' [ordinal_to_cutpoints()] converts the ordinal endpoints to the corresponding
#' cutpoints of the underlying latent continuous variable. Let
#' \eqn{P(x \le k) = G(c_k)} where \eqn{G} is the distribution function of the
#' latent variable. [ordinal_to_cutpoints()] converts \eqn{x} to \eqn{c_k} (or to
#' \eqn{c_{k - 1}} if `strict = TRUE`.
#'
#' @param x Integer vector with values in `1:(length(cutpoints) + 1)`.
#' @param cutpoints The cutpoints on the latent scale corresponding to
#' \eqn{\boldsymbol{c} = c(c_1, \cdots, c_{K - 1})}.
#'
#' @return Numeric vector with cutpoints corresponding to the values in `x`.
ordinal_to_cutpoints = function(x, cutpoints, strict) {
  if (strict) {
    # Add c_0 to the cutpoints
    cutpoints = c(-Inf, cutpoints)
  }
  else {
    # Add c_K to the cutpoints
    cutpoints = c(cutpoints, +Inf)
  }
  return(cutpoints[x])
}
