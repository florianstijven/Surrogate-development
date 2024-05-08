

#' Compute surrogacy measures for a time-to-event surrogate and a time-to-event true endpoint in the meta-analytic multiple-trial setting.
#'
#'The function 'MetaAnalyticSurvSurv()' fits the model for a time-to-event surrogate and time-to-event true endpoint developed by Burzykowski et al. (2001) in the meta-analytic multiple-trial setting.
#'
#' @details
#'
#' # Model
#'
#' In the model developed by Burzykowski et al. (2001), a copula-based model is used for the true time-to-event endpoint and the surrogate time-to-event endpoint.
#' More specifically, three copulas can be used: the Clayton copula, Hougaard copula and Plackett copula. The marginal model for the true and surrogate endpoint is the proportional hazard model.
#' The quality of the surrogate at the individual level can be evaluated by looking at the copula parameter.
#' The quality of the surrogate at the trial level can be evaluated by considering the correlation coefficient between the estimated treatment effects.
#'
#' # Data Format
#'
#' The data frame must contains the following columns:
#'
#' * a column with the observed time-to-event for the true endpoint
#' * a column with the time-to-event indicator for the true endpoint: 1 if true event is observed, 0 otherwise
#' * a column with the observed time-to-event for the surrogate endpoint
#' * a column with the time-to-event indicator for the surrogate endpoint: 1 if true event is observed, 0 otherwise
#' * a column with the treatment indicator: 0 or 1
#' * a column with the trial indicator
#' * a column with the center indicator. If there are no different centers within each trial, the center indicator is equal to the trial indicator
#' * a column with the patient indicator
#'
#' @references Burzykowski T, Molenberghs G, Buyse M, Geys H, Renard D (2001). “Validation of surrogate end points in multiple randomized clinical trials with failure time end points.” Journal of the Royal Statistical Society Series C: Applied Statistics, 50(4), 405–422
#'
#' @param data A data frame with the correct columns (See details).
#' @param true Observed time-to-event for true endpoint.
#' @param trueind Time-to-event indicator for the true endpoint.
#' @param surrog Observed time-to-event for surrogate endpoint.
#' @param surrogind Time-to-event indicator for the surrogate endpoint.
#' @param trt Treatment indicator.
#' @param center Center indicator (equal to trial if there are no different centers). This is the unit for which specific treatment effects are estimated.
#' @param trial Trial indicator. This is the unit for which common baselines are to be used.
#' @param patientid Patient indicator.
#' @param copula The copula that is used, either "Clayton", "Hougaard" or "Plackett"
#' @param adjustment The adjustment that should be made for the trial-level surrogacy, either "unadjusted", "weighted" or "adjusted"
#'
#' @return Returns an object of class "MetaAnalyticSurvSurv" that can be used to evaluate surrogacy and contains the following elements:
#'
#' * Indiv.Surrogacy: a data frame that contains the measure for the individual level surrogacy and 95% confidence interval.
#' * Trial.R2: a data frame that contains the correlation coefficient and 95% confidence interval to evaluate surrogacy at the trial level.
#' * EstTreatEffects: a data frame that contains the estimated treatment effects and sample size for each trial.
#' * fit_output: output of the maximization procedure (nlm) to maximize the likelihood.
#'
#' @export
#'
#' @author Dries De Witte
#'
#' @examples
#' \dontrun{
#' data("colorectal4")
#' fit <- MetaAnalyticSurvCat(data = colorectal4, true = truend, trueind = trueind, surrog = surrogend,
#'                            trt = treatn, center = center, trial = trialend, patientid = patid,
#'                            adjustment="unadjusted")
#' print(fit)
#' summary(fit)
#' plot(fit)
#' }
MetaAnalyticSurvSurv <- function(data, true, trueind, surrog, surrogind,
                     trt, center, trial, patientid, copula, adjustment) {
  data_m <- data
  data_m$surv <- data[[substitute(true)]]
  data_m$survind <- data[[substitute(trueind)]]
  data_m$pfs <- data[[substitute(surrog)]]
  data_m$pfsind <- data[[substitute(surrogind)]]
  data_m$treat <- data[[substitute(trt)]]
  data_m$center <- data[[substitute(center)]]
  data_m$trial <- data[[substitute(trial)]]
  data_m$patientid <- data[[substitute(patientid)]]

  data_m <- data_m %>%
    select(patientid, treat, center, survind, pfsind, trial, pfs, surv)

  #some functions
  MetaAnalyticSurvSurv_clayton <- function(data_m) {

    #clayton copula
    output <- 1
    lndelta <- 2
    beta0T <- 0
    beta0S <- 0
    pT <- 0.1
    pS <- 0.1
    llT <- 0
    llS <- 0
    betaT <- 0
    betaS <- 0

    # Read variables into vectors and matrices
    t1 <- as.vector(data_m$surv)
    delta1 <- as.vector(data_m$survind)
    t2 <- as.vector(data_m$pfs)
    delta2 <- as.vector(data_m$pfsind)
    x <- as.matrix(data_m$treat)
    center <- as.vector(data_m$center)
    trial <- as.vector(data_m$trial)
    xx <- 0
    z <- x

    cents <- unique(center)
    numcents <- length(cents)
    tri <- unique(trial)
    numtri <- length(tri)
    q <- ncol(x)
    r <- 0

    loglik <- function(param) {

      alpha <- exp(param[1]) + 1
      if (r > 0) {
        beta01 <- param[2:(1 + r)]
        beta02 <- param[(2 + r):(1 + 2 * r)]
      }

      lik <- 0
      for (i in 1:numcents) {
        places <- which(center == cents[i])
        hulp <- as.numeric(unique(trial[places]))
        stud <- which(tri == hulp)

        p1 <- param[1 + 2 * r + stud]
        p2 <- param[1 + 2 * r + numtri + stud]
        ll1 <- param[1 + 2 * r + 2 * numtri + stud]
        ll2 <- param[1 + 2 * r + 3 * numtri + stud]

        beta_1 <- param[(1 + 2 * r + 4 * numtri + q * (i - 1) + 1):(1 + 2 * r + 4 * numtri + q * (i - 1) + q)]
        beta_2 <- param[(1 + 2 * r + 4 * numtri + q * numcents + q * (i - 1) + 1):(1 + 2 * r + 4 * numtri + q * numcents + q * (i - 1) + q)]
        beta1= as.matrix(c(beta_1))
        beta2= as.matrix(c(beta_2))

        lambda1 <- exp(ll1)
        lambda2 <- exp(ll2)
        lt1 <- lambda1 * as.matrix(t1[places])
        lt2 <- lambda2 * as.matrix(t2[places])
        d1 <- as.matrix(delta1[places])
        d2 <- as.matrix(delta2[places])

        r1 <- lambda1 * (lt1^(p1 - 1)) * exp(z[places] * beta1) * p1
        s1 <- exp(-((lt1^p1) * exp(z[places,] * beta1)))
        r2 <- lambda2 * (lt2^(p2 - 1)) * exp(z[places] * beta2) * p2
        s2 <- exp(-((lt2^p2) * exp(z[places] * beta2)))

        s <- s1 ^ (1 - alpha) + s2 ^ (1 - alpha)

        bigf <- (s - 1) ^ (1 / (1 - alpha))
        f <- alpha * ((s - 1) ^ (alpha / (1 - alpha) - 1)) * r1 * (s1 ^ (1 - alpha)) * r2 * (s2 ^ (1 - alpha))
        f1 <- ((s - 1) ^ (alpha / (1 - alpha))) * r1 * (s1 ^ (1 - alpha))
        f2 <- ((s - 1) ^ (alpha / (1 - alpha))) * r2 * (s2 ^ (1 - alpha))

        lik <- lik + sum(d1 * d2 * log(f) + (1 - d1) * (1 - d2) * log(bigf) + d1 * (1 - d2) * log(f1) + d2 * (1 - d1) * log(f2))

      }

      return(-lik)
    }

    initp <- function(alpha0) {

      if (r > 0) {
        param0 <- c(alpha0, matrix(rep(c(beta0S, beta0T), times = r), nrow = 1),
                    matrix(rep(c(pS, pT, llS, llT), times = numtri), nrow = 1),
                    matrix(rep(c(betaS, betaT), times = q * numcents), nrow = 1))
      } else {

        param0 <- c(alpha0, matrix(pS, nrow = 1, ncol = numtri, byrow = TRUE),
                    matrix(pT, nrow = 1, ncol = numtri, byrow = TRUE),
                    matrix(llS, nrow = 1, ncol = numtri, byrow = TRUE),
                    matrix(llT, nrow = 1, ncol = numtri, byrow = TRUE),
                    matrix(betaS, nrow = 1, ncol = q * numcents, byrow = TRUE),
                    matrix(betaT, nrow = 1, ncol = q * numcents, byrow = TRUE))
      }

      return(param0)
    }

    lnalpha0 <- lndelta
    par0 <- initp(lnalpha0)

    nlm_output <- nlm(loglik, par0 , hessian=TRUE, iterlim = 1000)

    est <- nlm_output$estimate

    hess<-nlm_output$hessian
    fisher_info<-solve(hess)
    se <-sqrt(diag(fisher_info))

    lo <- est - qnorm(0.975) * se
    up <- est + qnorm(0.975) * se

    #### INDIVIDUAL LEVEL SURROGACY ####
    # For Log(alpha-1)
    alph <- est[1]
    se_alph <- se[1]
    lo_alph <- alph - qnorm(0.975) * se_alph
    up_alph <- alph + qnorm(0.975) * se_alph

    # For alpha
    alpha <- exp(est[1]) + 1

    # Calculate lower and upper bounds for alpha based on transforming 95% CI
    lo_al1 <- exp(lo[1]) + 1
    up_al1 <- exp(up[1]) + 1

    # Calculate standard error for alpha-
    se_alpha <- exp(est[1]) * se[1]

    # Calculate lower and upper bounds for alpha based on delta-method
    lo_al2 <- alpha - qnorm(0.975) * se_alpha
    up_al2 <- alpha + qnorm(0.975) * se_alpha

    alpha_df <- data.frame(cbind(alpha, lo_al1, up_al1), stringsAsFactors = TRUE)
    colnames(alpha_df) <- c("Copula parameter alpha", "CI lower limit",
                            "CI upper limit")
    rownames(alpha_df) <- c(" ")

    alpha_delta_df <- data.frame(cbind(alpha, lo_al2, up_al2), stringsAsFactors = TRUE)
    colnames(alpha_delta_df) <- c("Copula parameter alpha", "CI lower limit",
                                  "CI upper limit")
    rownames(alpha_delta_df) <- c(" ")

    # Calculate tau
    tau <- (alpha - 1) / (alpha + 1)

    # Calculate standard error for tau
    se_tau <- se[1] * 2 * (alpha - 1) / ((alpha + 1) ^ 2)

    # Calculate lower and upper bounds for tau
    lo_tau <- tau - qnorm(0.975) * se_tau
    up_tau <- tau + qnorm(0.975) * se_tau

    tau_delta_df <- data.frame(cbind(tau, lo_tau, up_tau), stringsAsFactors = TRUE)
    colnames(tau_delta_df) <- c("Kendall's tau", "CI lower limit",
                                "CI upper limit")
    rownames(tau_delta_df) <- c(" ")

    #### TRIAL LEVEL SURROGACY ####
    coef_se <- cbind(est, se)
    coef_se

    endp1 <- coef_se[(1 + 4 * numtri + 1):(1 + 4 * numtri + q * numcents),]
    endp2 <- coef_se[(1 + 4 * numtri + q * numcents + 1):(1 + 4 * numtri + 2 * q * numcents),]

    cova <- diag(fisher_info[(1 + 4 * numtri + 1):(1 + 4 * numtri + q * numcents),
                             (1 + 4 * numtri + q * numcents + 1):(1 + 4 * numtri + 2 * q * numcents)])

    weight <- matrix(0, nrow = numcents, ncol = 1)
    for (i in 1:numcents) {
      weight[i, ] <- nrow(as.matrix(which(center == cents[i])))
    }

    memo <- cbind(cents, endp1, endp2, cova, weight)
    colnames(memo) <- c("center", "surv_est", "surv_se", "pfs_est", "pfs_se", "cova", "weight")
    memo <- as.data.frame(memo)

    output_list <- list(copula_parameter = alpha_delta_df,
                        Indiv.Surrogacy = tau_delta_df, memo=memo, nlm_output=nlm_output)
    return(output_list)


  }

  MetaAnalyticSurvSurv_hougaard <- function(data_m) {

    #hougaard copula
    output <- 1
    lndelta <- -1 #smaller than 0
    beta0T <- 0
    beta0S <- 0
    pT <- 0.1
    pS <- 0.1
    llT <- 0
    llS <- 0
    betaT <- 0
    betaS <- 0

    # Read variables into vectors and matrices
    t1 <- as.vector(data_m$surv)
    delta1 <- as.vector(data_m$survind)
    t2 <- as.vector(data_m$pfs)
    delta2 <- as.vector(data_m$pfsind)
    x <- as.matrix(data_m$treat)
    center <- as.vector(data_m$center)
    trial <- as.vector(data_m$trial)
    xx <- 0
    z <- x

    cents <- unique(center)
    numcents <- length(cents)
    tri <- unique(trial)
    numtri <- length(tri)
    q <- ncol(x)
    r <- 0

    loglik <- function(param) {

      alpha <- exp(param[1])

      lik <- 0
      for (i in 1:numcents) {
        places <- which(center == cents[i])
        hulp <- as.numeric(unique(trial[places]))
        stud <- which(tri == hulp)
        p1 <- param[1 + 2 * r + stud]
        p2 <- param[1 + 2 * r + numtri + stud]
        ll1 <- param[1 + 2 * r + 2 * numtri + stud]
        ll2 <- param[1 + 2 * r + 3 * numtri + stud]
        beta_1 <- param[(1 + 2 * r + 4 * numtri + q * (i - 1) + 1):(1 + 2 * r + 4 * numtri + q * (i - 1) + q)]
        beta_2 <- param[(1 + 2 * r + 4 * numtri + q * numcents + q * (i - 1) + 1):(1 + 2 * r + 4 * numtri + q * numcents + q * (i - 1) + q)]
        lambda1 <- exp(ll1)
        lambda2 <- exp(ll2)
        beta1= as.matrix(c(beta_1))
        beta2= as.matrix(c(beta_2))

        zc <- as.vector(z[places])
        d1 <- as.matrix(delta1[places])
        d2 <- as.matrix(delta2[places])
        lt1 <- lambda1 * as.matrix(t1[places])
        lt2 <- lambda2 * as.matrix(t2[places])

        r1 <- lambda1 * (lt1^(p1/alpha - 1)) * exp(zc * beta1/alpha) * p1/alpha
        cr1 <- (lt1^(p1/alpha)) * exp(zc * beta1/alpha)
        r2 <- lambda2 * (lt2^(p2/alpha - 1)) * exp(zc * beta2/alpha) * p2/alpha
        cr2 <- (lt2^(p2/alpha)) * exp(zc * beta2/alpha)

        s <- cr1+cr2

        bigf <- exp(-s^alpha)

        f <- bigf*alpha*(1-alpha+alpha*(s^alpha))*r1*r2*(s^(alpha-2))
        f1 <- bigf*alpha*(s^(alpha-1))*r1
        f2 <- bigf*alpha*(s^(alpha-1))*r2

        lik <- lik + sum(d1 * d2 * log(f) + (1 - d1) * (1 - d2) * log(bigf) + d1 * (1 - d2) * log(f1) + d2 * (1 - d1) * log(f2))

      }

      return(-lik)
    }

    initp <- function(alpha0) {

      param0 <- c(alpha0, matrix(pS, nrow = 1, ncol = numtri, byrow = TRUE),
                  matrix(pT, nrow = 1, ncol = numtri, byrow = TRUE),
                  matrix(llS, nrow = 1, ncol = numtri, byrow = TRUE),
                  matrix(llT, nrow = 1, ncol = numtri, byrow = TRUE),
                  matrix(betaS, nrow = 1, ncol = q * numcents, byrow = TRUE),
                  matrix(betaT, nrow = 1, ncol = q * numcents, byrow = TRUE))

      return(param0)
    }

    lnalpha0 <- lndelta
    par0 <- initp(lnalpha0)

    nlm_output <- nlm(loglik, par0 , hessian=TRUE, iterlim = 1000)

    est <- nlm_output$estimate

    hess<-nlm_output$hessian
    fisher_info<-solve(hess)
    se <-sqrt(diag(fisher_info))

    lo <- est - qnorm(0.975) * se
    up <- est + qnorm(0.975) * se

    #### INDIVIDUAL LEVEL SURROGACY ####
    # For Log(alpha-1)
    alph <- est[1]
    se_alph <- se[1]
    lo_alph <- alph - qnorm(0.975) * se_alph
    up_alph <- alph + qnorm(0.975) * se_alph

    # For alpha
    alpha <- exp(est[1])

    # Calculate lower and upper bounds for alpha based on transforming 95% CI
    lo_al1 <- exp(lo[1])
    up_al1 <- exp(up[1])

    # Calculate standard error for alpha-
    se_alpha <- exp(est[1]) * se[1]

    # Calculate lower and upper bounds for alpha based on delta-method
    lo_al2 <- alpha - qnorm(0.975) * se_alpha
    up_al2 <- alpha + qnorm(0.975) * se_alpha

    alpha_df <- data.frame(cbind(alpha, lo_al1, up_al1), stringsAsFactors = TRUE)
    colnames(alpha_df) <- c("Copula parameter alpha", "CI lower limit",
                            "CI upper limit")
    rownames(alpha_df) <- c(" ")

    alpha_delta_df <- data.frame(cbind(alpha, lo_al2, up_al2), stringsAsFactors = TRUE)
    colnames(alpha_delta_df) <- c("Copula parameter alpha", "CI lower limit",
                                  "CI upper limit")
    rownames(alpha_delta_df) <- c(" ")

    # Calculate tau
    tau <- 1-alpha

    # Calculate standard error for tau
    se_tau <- se_alpha

    # Calculate lower and upper bounds for tau
    lo_tau <- tau - qnorm(0.975) * se_tau
    up_tau <- tau + qnorm(0.975) * se_tau

    tau_delta_df <- data.frame(cbind(tau, lo_tau, up_tau), stringsAsFactors = TRUE)
    colnames(tau_delta_df) <- c("Kendall's tau", "CI lower limit",
                                "CI upper limit")
    rownames(tau_delta_df) <- c(" ")

    #### TRIAL LEVEL SURROGACY ####
    coef_se <- cbind(est, se)
    coef_se

    endp1 <- coef_se[(1 + 4 * numtri + 1):(1 + 4 * numtri + q * numcents),]
    endp2 <- coef_se[(1 + 4 * numtri + q * numcents + 1):(1 + 4 * numtri + 2 * q * numcents),]

    cova <- diag(fisher_info[(1 + 4 * numtri + 1):(1 + 4 * numtri + q * numcents),
                             (1 + 4 * numtri + q * numcents + 1):(1 + 4 * numtri + 2 * q * numcents)])

    weight <- matrix(0, nrow = numcents, ncol = 1)
    for (i in 1:numcents) {
      weight[i, ] <- nrow(as.matrix(which(center == cents[i])))
    }

    memo <- cbind(cents, endp1, endp2, cova, weight)
    colnames(memo) <- c("center", "surv_est", "surv_se", "pfs_est", "pfs_se", "cova", "weight")
    memo <- as.data.frame(memo)

    output_list <- list(copula_parameter = alpha_delta_df,
                        Indiv.Surrogacy = tau_delta_df, memo=memo, nlm_output=nlm_output)
    return(output_list)


  }

  MetaAnalyticSurvSurv_plackett <- function(data_m) {

    #hougaard copula
    theta <- 1.5
    beta0T <- 0
    beta0S <- 0
    pT <- 0.1
    pS <- 0.1
    llT <- 0
    llS <- 0
    betaT <- 0
    betaS <- 0

    # Read variables into vectors and matrices
    t <- as.vector(data_m$surv)
    deltat <- as.vector(data_m$survind)
    s <- as.vector(data_m$pfs)
    deltas <- as.vector(data_m$pfsind)
    x <- as.matrix(data_m$treat)
    center <- as.vector(data_m$center)
    trial <- as.vector(data_m$trial)
    xx <- 0
    zt <- x
    zs <- x

    cents <- unique(center)
    numcents <- length(cents)
    tri <- unique(trial)
    numtri <- length(tri)
    qq <- ncol(x)
    r <- 0

    q1 <- function(u, v, theta) {
      quv <- sqrt((1 + (theta - 1) * (u + v))^2 - 4 * u * v * theta * (theta - 1))
      return(quv)
    }

    d1qth <- function(u, v, theta) {
      dq <- ((1 + (theta - 1) * (u + v)) * (u + v) - 2 * u * v * (2 * theta - 1)) / q1(u, v, theta)
      return(dq)
    }

    d1qu <- function(u, v, theta) {
      dq <- (theta - 1) * (1 + (theta - 1) * (u + v) - 2 * v * theta) / q1(u, v, theta)
      return(dq)
    }


    d1qv <- function(u, v, theta) {
      dq <- (theta - 1) * (1 + (theta - 1) * (u + v) - 2 * u * theta) / q1(u, v, theta)
      return(dq)
    }

    d2qth2 <- function(u, v, theta) {
      d2q <- ((u - v)^2 - d1qth(u, v, theta)^2) / q1(u, v, theta)
      return(d2q)
    }

    d2qu2 <- function(u, v, theta) {
      d2q <- ((theta - 1)^2 - d1qu(u, v, theta)^2) / q1(u, v, theta)
      return(d2q)
    }

    d2qv2 <- function(u, v, theta) {
      d2q <- ((theta - 1)^2 - d1qv(u, v, theta)^2) / q1(u, v, theta)
      return(d2q)
    }

    d2quv <- function(u, v, theta) {
      d2q <- -((theta - 1) * (theta + 1) + d1qu(u, v, theta) * d1qv(u, v, theta)) / q1(u, v, theta)
      return(d2q)
    }

    d2quth <- function(u, v, theta) {
      q <- q1(u, v, theta)
      d1qu <- d1qu(u, v, theta)
      d1qth <- d1qth(u, v, theta)
      d2q <- (theta - 1) * (u - v) / q + d1qu * (1 / (theta - 1) - d1qth / q)
      return(d2q)
    }

    d2qvth <- function(u, v, theta) {
      q <- q1(u, v, theta)
      d1qv <- d1qv(u, v, theta)
      d1qth <- d1qth(u, v, theta)
      d2q <- (theta - 1) * (v - u) / q + d1qv * (1 / (theta - 1) - d1qth / q)
      return(d2q)
    }

    C <- function(u, v, theta) {
      if (theta != 1)
        Cuv <- (1 + (theta - 1) * (u + v) - q1(u, v, theta)) / (2 * (theta - 1))
      else
        Cuv <- u * v
      return(Cuv)
    }

    d1Cth <- function(u, v, theta) {
      if (theta != 1)
        dC <- (-C(u, v, theta) + (u + v - d1qth(u, v, theta)) / 2) / (theta - 1)
      else
        dC <- matrix(0, nrow(u), 1)
      return(dC)
    }

    d1Cv <- function(u, v, theta) {
      if (theta != 1)
        dC <- (1 - d1qv(u, v, theta) / (theta - 1)) / 2
      else
        dC <- u
      return(dC)
    }

    d1Cu <- function(u, v, theta) {
      if (theta != 1)
        dC <- (1 - d1qu(u, v, theta) / (theta - 1)) / 2
      else
        dC <- v
      return(dC)
    }

    d2Cu2 <- function(u, v, theta) {
      if (theta != 1)
        d2C <- -(theta - 1) * (1 - (d1qu(u, v, theta) / (theta - 1))^2) / (2 * q1(u, v, theta))
      else
        d2C <- matrix(0, nrow(u), 1)
      return(d2C)
    }

    d2Cv2 <- function(u, v, theta) {
      if (theta != 1)
        d2C <- -(theta - 1) * (1 - (d1qv(u, v, theta) / (theta - 1))^2) / (2 * q1(u, v, theta))
      else
        d2C <- matrix(0, nrow(u), 1)
      return(d2C)
    }

    d2Cth2 <- function(u, v, theta) {
      if (theta != 1)
        d2C <- -(2 * d1Cth(u, v, theta) + d2qth2(u, v, theta) / 2) / (theta - 1)
      else
        d2C <- matrix(0, nrow(u), 1)
      return(d2C)
    }

    d2Cuv <- function(u, v, theta, ...) {
      if (theta != 1)
        d2C <- (theta + 1 + d1qu(u, v, theta) * d1qv(u, v, theta) / (theta - 1)) / (2 * q1(u, v, theta))
      else
        d2C <- matrix(1, nrow(u), 1)
      return(d2C)
    }

    d2Cuth <- function(u, v, theta) {
      if (theta != 1)
        d2C <- (v - u + d1qu(u, v, theta) * d1qth(u, v, theta) / (theta - 1)) / (2 * q1(u, v, theta))
      else
        d2C <- matrix(0, nrow(u), 1)
      return(d2C)
    }

    d2Cvth <- function(u, v, theta) {
      if (theta != 1)
        d2C <- (u - v + d1qv(u, v, theta) * d1qth(u, v, theta) / (theta - 1)) / (2 * q1(u, v, theta))
      else
        d2C <- matrix(0, nrow(u), 1)
      return(d2C)
    }

    d3Cu3 <- function(u, v, theta) {
      if (theta != 1) {
        q_val <- q1(u, v, theta)
        dqu <- d1qu(u, v, theta)
        d2qu2 <- d2qu2(u, v, theta)
        d2Cu2 <- d2Cu2(u, v, theta)
        d3C <- (-d2Cu2 + d2qu2 / (theta - 1)) * dqu / q_val
      } else {
        d3C <- matrix(0, nrow(u), 1)
      }
      return(d3C)
    }

    d3Cv3 <- function(u, v, theta) {
      if (theta != 1) {
        q_val <- q1(u, v, theta)
        dqv <- d1qv(u, v, theta)
        d2qv2 <- d2qv2(u, v, theta)
        d2Cv2 <- d2Cv2(u, v, theta)
        d3C <- (-d2Cv2 + d2qv2 / (theta - 1)) * dqv / q_val
      } else {
        d3C <- matrix(0, nrow(u), 1)
      }
      return(d3C)
    }

    d3Cu2v <- function(u, v, theta) {
      if (theta != 1) {
        q_val <- q1(u, v, theta)
        dqv <- d1qv(u, v, theta)
        dqu <- d1qu(u, v, theta)
        d2quv <- d2quv(u, v, theta)
        d2Cu2 <- d2Cu2(u, v, theta)
        d3C <- (-d2Cu2 * dqv + dqu * d2quv / (theta - 1)) / q_val
      } else {
        d3C <- matrix(0, nrow(u), 1)
      }
      return(d3C)
    }

    d3Cuv2 <- function(u, v, theta) {
      if (theta != 1) {
        q_val <- q1(u, v, theta)
        dqv <- d1qv(u, v, theta)
        dqu <- d1qu(u, v, theta)
        d2quv <- d2quv(u, v, theta)
        d2Cv2 <- d2Cv2(u, v, theta)
        d3C <- (-d2Cv2 * dqu + dqv * d2quv / (theta - 1)) / q_val
      } else {
        d3C <- matrix(0, nrow(u), 1)
      }
      return(d3C)
    }

    d3Cuth2 <- function(u, v, theta) {
      if (theta != 1) {
        q <- q1(u, v, theta)
        dqu <- d1qu(u, v, theta)
        dqth <- d1qth(u, v, theta)
        d2quth <- d2quth(u, v, theta)
        d2qth2 <- d2qth2(u, v, theta)
        d2Cuth <- d2Cuth(u, v, theta)
        d3C <- -d2Cuth * (dqth / q + 1 / (theta - 1)) + (dqu * d2qth2 + d2quth * dqth - (u - v)) / (2 * q * (theta - 1))
      } else {
        d3C <- matrix(0, nrow(u), 1)
      }
      return(d3C)
    }

    d3Cu2th <- function(u, v, theta) {
      if (theta != 1) {
        q <- q1(u, v, theta)
        dqu <- d1qu(u, v, theta)
        dqth <- d1qth(u, v, theta)
        d2quth <- d2quth(u, v, theta)
        d2Cu2 <- d2Cu2(u, v, theta)
        d3C <- d2Cu2 * (1 / (theta - 1) - dqth / q) + dqu * d2quth / (q * (theta - 1))
      } else {
        d3C <- matrix(0, nrow(u), 1)
      }
      return(d3C)
    }

    d3Cuvth <- function(u, v, theta) {
      if (theta != 1) {
        q <- q1(u, v, theta)
        dqu <- d1qu(u, v, theta)
        dqv <- d1qv(u, v, theta)
        dqth <- d1qth(u, v, theta)
        d2quth <- d2quth(u, v, theta)
        d2qvth <- d2qvth(u, v, theta)
        d2Cuv <- d2Cuv(u, v, theta)
        d3C <- -d2Cuv * (dqth / q + 1 / (theta - 1)) + (d2quth * dqv + d2qvth * dqu + 2 * theta) / (2 * q * (theta - 1))
      } else {
        d3C <- matrix(0, nrow(u), 1)
      }
      return(d3C)
    }

    Ft <- function(lambdat, pt, lambdas, ps, betat, betas, t, s, zt, zs) {
      lt <- lambdat * t
      Ft <- 1 - exp(-(lt^pt) * exp(zt %*% betat))
      return(Ft)
    }

    pft <- function(lambdat, pt, lambdas, ps, betat, betas, t, s, zt, zs) {
      lt <- lambdat * t
      f <- lambdat * pt * (lt^(pt - 1)) * exp(zt %*% betat) * exp(-(lt^pt) * exp(zt %*% betat))
      return(f)
    }

    Fs <- function(lambdat,pt,lambdas,ps,betat,betas,t,s,zt,zs) {
      ls <- lambdas*s
      Fs <- 1-exp(-(ls^ps)*exp(zs %*%betas))
      return(Fs)
    }

    pfs <- function(lambdat,pt,lambdas,ps,betat,betas,t,s,zt,zs) {
      ls <- lambdas*s
      f <- lambdas*ps*(ls^(ps-1))*exp(zs %*%betas)*exp(-(ls^ps)*exp(zs %*%betas))
    }


    loglik <- function(param) {
      theta <- param[1]

      lik <- 0
      for (i in 1:numcents) {
        Ti <- t[center == cents[i]]
        Si <- s[center == cents[i]]
        deltaTi <- deltat[center == cents[i]]
        deltaSi <- deltas[center == cents[i]]
        zTi <- as.matrix(zt[center == cents[i], ])
        zSi <- as.matrix(zs[center == cents[i], ])
        xi <- x[center == cents[i], ]

        hulp <- as.numeric(unique(trial[center == cents[i]]))
        stud <- which(tri == hulp)
        pTi <- param[1 + stud]
        pSi <- param[1 + numtri + stud]
        llTi <- param[1 + 2*numtri + stud]
        llSi <- param[1 + 3*numtri + stud]
        efftrtT <- param[(1 + 4 * numtri + qq*(i-1)+1):(1+4*numtri+qq*(i-1)+qq)]
        efftrtS <- param[(1 + 4 * numtri + qq*numcents + qq*(i-1)+1):(1+4*numtri+qq*numcents+qq*(i-1)+qq)]
        betaTi <- as.matrix(cbind(efftrtT))
        betaSi <- as.matrix(cbind(efftrtS))

        lambdaTi <- exp(llTi)
        lti <- lambdaTi * t[center == cents[i]]
        lambdaSi <- exp(llSi)
        lsi <- lambdaSi * s[center == cents[i]]

        Fti <- Ft(lambdaTi,pTi,lambdaSi,pSi,betaTi,betaSi,Ti,Si,zTi,zSi)
        pfti <- pft(lambdaTi,pTi,lambdaSi,pSi,betaTi,betaSi,Ti,Si,zTi,zSi)
        Sti <- 1-Fti
        Fsi <- Fs(lambdaTi,pTi,lambdaSi,pSi,betaTi,betaSi,Ti,Si,zTi,zSi)
        pfsi <- pfs(lambdaTi,pTi,lambdaSi,pSi,betaTi,betaSi,Ti,Si,zTi,zSi)
        Ssi <- 1-Fsi

        q <- q1(Sti,Ssi,theta)
        d1qth <- d1qth(Sti,Ssi,theta)
        d1qu <- d1qu(Sti,Ssi,theta)
        d1qv <- d1qv(Sti,Ssi,theta)
        d2qth2 <- d2qth2(Sti,Ssi,theta)

        d1Cth <- d1Cth(Sti,Ssi,theta)
        d1Cu <- d1Cu(Sti,Ssi,theta)
        d1Cv <- d1Cv(Sti,Ssi,theta)
        d2Cuv <- d2Cuv(Sti,Ssi,theta,q,d1qth,d1qu,d1qv,d2qth2,d1Cth)

        fi <- d2Cuv*pfti*pfsi
        Surv_i <- C(Sti,Ssi,theta)
        dtSi <- -d1Cu*pfti
        dsSi <- -d1Cv*pfsi

        lik <- lik + sum(deltaTi * deltaSi * log(fi) +
                           deltaTi * (1 - deltaSi) * log(-dtSi) +
                           (1 - deltaTi) * deltaSi * log(-dsSi) +
                           (1 - deltaTi) * (1 - deltaSi) * log(Surv_i))
      }

      return(-lik)
    }

    initp <- function(theta0) {

      param0 <- c(theta0, matrix(pS, nrow = 1, ncol = numtri, byrow = TRUE),
                  matrix(pT, nrow = 1, ncol = numtri, byrow = TRUE),
                  matrix(llS, nrow = 1, ncol = numtri, byrow = TRUE),
                  matrix(llT, nrow = 1, ncol = numtri, byrow = TRUE),
                  matrix(betaS, nrow = 1, ncol = qq * numcents, byrow = TRUE),
                  matrix(betaT, nrow = 1, ncol = qq * numcents, byrow = TRUE))

      return(param0)
    }

    theta0 <- theta
    par0 <- initp(theta0)

    nlm_output <- nlm(loglik, par0 , hessian=TRUE, iterlim = 1000)

    est <- nlm_output$estimate

    hess<-nlm_output$hessian
    fisher_info<-solve(hess)
    se <-sqrt(diag(fisher_info))

    lo <- est - qnorm(0.975) * se
    up <- est + qnorm(0.975) * se

    #### INDIVIDUAL LEVEL SURROGACY ####
    theta <- est[1]
    se_th <- se[1]
    lo_th <- theta - qnorm(0.975) * se_th
    up_th <- theta + qnorm(0.975) * se_th

    theta_df <- data.frame(cbind(theta, lo_th, up_th), stringsAsFactors = TRUE)
    colnames(theta_df) <- c("Copula parameter theta", "CI lower limit",
                            "CI upper limit")
    rownames(theta_df) <- c(" ")

    #### TRIAL LEVEL SURROGACY ####
    coef_se <- cbind(est, se)
    coef_se

    endp1 <- coef_se[(1 + 4 * numtri + 1):(1 + 4 * numtri + qq * numcents),]
    endp2 <- coef_se[(1 + 4 * numtri + qq * numcents + 1):(1 + 4 * numtri + 2 * qq * numcents),]

    cova <- diag(fisher_info[(1 + 4 * numtri + 1):(1 + 4 * numtri + qq * numcents),
                             (1 + 4 * numtri + qq * numcents + 1):(1 + 4 * numtri + 2 * qq * numcents)])

    weight <- matrix(0, nrow = numcents, ncol = 1)
    for (i in 1:numcents) {
      weight[i, ] <- nrow(as.matrix(which(center == cents[i])))
    }

    memo <- cbind(cents, endp1, endp2, cova, weight)
    colnames(memo) <- c("center", "surv_est", "surv_se", "pfs_est", "pfs_se", "cova", "weight")
    memo <- as.data.frame(memo)

    rho <- ((theta+1)/(theta-1))-2*theta*log(theta)/((theta-1)^2)
    drho <- 2*(2*(1-theta)+(theta+1)*log(theta))/((theta-1)^3)

    var_rho <- (drho^2)*(se_th^2)
    se_rho <- sqrt(var_rho)
    lo_rho <- rho-qnorm(0.975)*se_rho
    up_rho <- rho+qnorm(0.975)*se_rho

    rho_df <- data.frame(cbind(rho, lo_rho, up_rho), stringsAsFactors = TRUE)
    colnames(rho_df) <- c("Spearmans rho", "CI lower limit",
                          "CI upper limit")
    rownames(rho_df) <- c(" ")

    output_list <- list(copula_parameter = theta_df, Indiv.Surrogacy=rho_df, memo=memo, fit_output=nlm_output)
    return(output_list)


  }

  trial_level_secondstage_unadjusted <- function(memo) {
    surv_eff <- subset(memo, select = c(center, surv_est, weight))
    colnames(surv_eff)[2] <- "effect"
    surv_eff$endp <- "MAIN"

    pfs_eff <- subset(memo, select = c(center, pfs_est, weight))
    colnames(pfs_eff)[2] <- "effect"
    pfs_eff$endp <- "SURR"

    # Combine surv_eff and pfs_eff data frames
    shihco <- rbind(surv_eff, pfs_eff)
    shihco$center <- as.numeric(shihco$center)

    # Sort the combined data frame by center and endp
    shihco <- shihco[order(shihco$center, shihco$endp), ]
    shihco$endp <- as.factor(shihco$endp)
    shihco$effect <- as.numeric(shihco$effect)

    invisible(capture.output(model <- nlme::gls(effect ~ -1 + factor(endp), data = shihco,
                                                correlation = nlme::corCompSymm(form = ~ 1 | center),
                                                weights = nlme::varIdent(form = ~ 1 | endp),
                                                method = "ML",
                                                control = nlme::glsControl(maxIter = 25, msVerbose = TRUE))))

    rsquared <- nlme::intervals(model, which = "var-cov")$corStruct^2
    R2 <- as.vector(rsquared)[2]
    lo_R2 <- as.vector(rsquared)[1]
    up_R2 <- as.vector(rsquared)[3]

    Trial.R2 <- data.frame(cbind(R2, lo_R2, up_R2), stringsAsFactors = TRUE)
    colnames(Trial.R2) <- c("R2 Trial (unadjusted)", "CI lower limit",
                                  "CI upper limit")
    rownames(Trial.R2) <- c(" ")


    return(list(Trial.R2=Trial.R2))

  }

  trial_level_secondstage_weighted <- function(memo) {
    model <- lm(surv_est ~ pfs_est,
                data =memo, weights = weight)

    R2 <- summary(model)$r.squared

    ### based on cortinas
    Trial.ID <- memo$center
    N <- length(Trial.ID)
    R2.sd <- sqrt((4 * R2 * (1 - R2)^2)/(N - 3))
    lo_R2 <- max(0, R2 + qnorm(0.05/2) * (R2.sd))
    up_R2 <- min(1, R2 + qnorm(1 - 0.05/2) * (R2.sd))

    Trial.R2.weighted <- data.frame(cbind(R2, lo_R2, up_R2), stringsAsFactors = TRUE)
    colnames(Trial.R2) <- c("R2 Trial (weighted)", "CI lower limit",
                                     "CI upper limit")
    rownames(Trial.R2) <- c(" ")

    return(list(Trial.R2=Trial.R2))

  }

  trial_level_secondstage_adjusted <- function(memo) {
    #### some additional functions
    R2TrialFun <- function(D) {
      R2.trial <- D[1, 2] ^ 2 / prod(diag(D))
      return(R2.trial)
    }
    R2IndFun <- function(Sigma) {
      R2.ind <- Sigma[1, 2]^2/(Sigma[1, 1] * Sigma[2, 2])
      return(R2.ind)
    }
    pdDajustment = function(D) {
      EigenD <- eigen(D)
      E <- diag(EigenD$values)
      diag(E)[diag(E) <= 0] <- 1e-04
      L <- EigenD$vector
      DH <- L %*% tcrossprod(E, L)
      return(DH)
    }

    n <- as.numeric(memo$weight)
    Trial.ID <- memo$center
    N <- length(Trial.ID)

    BetaH_i = memo[,c(2,4)]
    a_i = n/sum(n)
    a_i.2 = (n - 2)/sum(n - 2)
    BetaH = apply(t(BetaH_i), 1, weighted.mean, w = a_i)
    N = length(a_i)
    b_i <- t(BetaH_i) - tcrossprod(BetaH, matrix(1, N, 1))
    Sb = tcrossprod(b_i)
    num = 1 - 2 * a_i + sum(a_i^2)
    denom = sum(num)


    R.list <- list()

    for (i in 1:nrow(memo)) {
      se1_sq <- memo$surv_se[i]^2
      se2_sq <- memo$pfs_se[i]^2
      cov <- memo$cova[i]

      var_cov_matrix <- matrix(c(se1_sq, cov, cov, se2_sq), nrow = 2, byrow = TRUE)

      R.list[[i]] <- var_cov_matrix
    }

    DH = (Sb - Reduce("+", mapply(function(X, x) {
      X * x
    }, R.list, num, SIMPLIFY = F)))/denom
    DH.pd = min(eigen(DH, only.values = T)$values) > 0
    if (!DH.pd) {
      DH = pdDajustment(DH)
      warning(paste("The estimate of D is non-positive definite. Adjustment for non-positive definiteness was required"))
    }
    R2 <- DH[1, 2] ^ 2 / prod(diag(DH))
    R2.sd <- sqrt((4 * R2 * (1 - R2)^2)/(N - 3))
    Alpha <- 0.05
    lo_R2 <- max(0, R2 + qnorm(Alpha/2) * (R2.sd))
    up_R2 <- min(1, R2 + qnorm(1 - Alpha/2) * (R2.sd))
    Trial.R2 <- data.frame(cbind(R2, lo_R2, up_R2), stringsAsFactors = TRUE)

    colnames(Trial.R2) <- c("R2 Trial (adjusted)", "CI lower limit",
                                "CI upper limit")
    rownames(Trial.R2) <- c(" ")

    return(list(Trial.R2=Trial.R2))

  }


  if(copula == "Clayton") {
    message("Fitting the model with Clayton copula...")
    stage1 <- MetaAnalyticSurvSurv_clayton(data_m = data_m)

    if(adjustment == "adjusted"){
      stage2 <- trial_level_secondstage_adjusted(memo=stage1$memo)
    }

    if(adjustment == "weighted"){
      stage2 <- trial_level_secondstage_weighted(memo=stage1$memo)
    }

    if(adjustment == "unadjusted"){
      stage2 <- trial_level_secondstage_unadjusted(memo=stage1$memo)
    }
  }

  if(copula == "Hougaard") {
    message("Fitting the model with Hougaard copula...")
    stage1 <- MetaAnalyticSurvSurv_hougaard(data_m = data_m)

    if(adjustment == "adjusted"){
      stage2 <- trial_level_secondstage_adjusted(memo=stage1$memo)
    }

    if(adjustment == "weighted"){
      stage2 <- trial_level_secondstage_weighted(memo=stage1$memo)
    }

    if(adjustment == "unadjusted"){
      stage2 <- trial_level_secondstage_unadjusted(memo=stage1$memo)
    }
  }

  if(copula == "Plackett") {
    message("Fitting the model with Plackett copula...")
    stage1 <- MetaAnalyticSurvSurv_plackett(data_m = data_m)

    if(adjustment == "adjusted"){
      stage2 <- trial_level_secondstage_adjusted(memo=stage1$memo)
    }

    if(adjustment == "weighted"){
      stage2 <- trial_level_secondstage_weighted(memo=stage1$memo)
    }

    if(adjustment == "unadjusted"){
      stage2 <- trial_level_secondstage_unadjusted(memo=stage1$memo)
    }
  }

  output_list <- list(Indiv.Surrogacy = stage1$Indiv.Surrogacy,
                      Trial.R2 = stage2$Trial.R2,
                      EstTreatEffects = stage1$memo,
                      fit_output = stage1$fit_output,
                      Call = match.call())

  class(output_list) <- "MetaAnalyticSurvSurv"

  return(output_list)

}

#' Provides a summary of the surrogacy measures for an object fitted with the 'MetaAnalyticSurvSurv()' function.
#'
#' @method summary MetaAnalyticSurvSurv
#'
#' @param object An object of class 'MetaAnalyticSurvSurv' fitted with the 'MetaAnalyticSurvSurv()' function.
#' @param ... ...
#'
#' @return The surrogacy measures with their 95% confidence intervals.
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal")
#' fit <- MetaAnalyticSurvSurv(data = colorectal4, true = truend, trueind = trueind, surrog = surrogend,
#'                trt = treatn, center = center, trial = trialend, patientid = patid, adjustment="unadjusted")
#' summary(fit)
#' }

summary.MetaAnalyticSurvSurv <- function(object,...){
  cat("Surrogacy measures with 95% confidence interval \n\n")
  cat("Individual level surrogacy: ", "\n\n")
  cat("Measure: ", sprintf("%.4f", object$Indiv.Surrogacy[1,1]), "[", sprintf("%.4f", object$Indiv.Surrogacy[1,2]),";", sprintf("%.4f", object$Indiv.Surrogacy[1,3]) , "]", "\n\n")
  cat("Trial level surrogacy: ", "\n\n")
  cat("R Square: ", sprintf("%.4f", object$Trial.R2[1,1]),"[", sprintf("%.4f", object$Trial.R2[1,2]),";", sprintf("%.4f", object$Trial.R2[1,3]) , "]", "\n\n")
}

#' Prints all the elements of an object fitted with the 'MetaAnalyticSurvSurv()' function.
#'
#' @method print MetaAnalyticSurvSurv
#'
#' @param x An object of class 'MetaAnalyticSurvSurv' fitted with the 'MetaAnalyticSurvSurv()' function.
#' @param ... ...
#'
#' @return The surrogacy measures with their 95% confidence intervals and the estimated treatment effect on the surrogate and true endpoint.
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal4")
#' fit <- MetaAnalyticSurvSurv(data = colorectal4, true = truend, trueind = trueind, surrog = surrogend,
#'                trt = treatn, center = center, trial = trialend, patientid = patid, adjustment="unadjusted")
#' print(fit)
#' }
print.MetaAnalyticSurvSurv <- function(x,...){
  cat("Surrogacy measures with 95% confidence interval \n\n")
  cat("Individual level surrogacy: ", "\n\n")
  cat("Global Odds: ", sprintf("%.4f", x$Indiv.GlobalOdds[1,1]), "[", sprintf("%.4f", x$Indiv.GlobalOdds[1,2]),";", sprintf("%.4f", x$Indiv.GlobalOdds[1,3]) , "]", "\n\n")
  cat("Trial level surrogacy: ", "\n\n")
  cat("R Square: ", sprintf("%.4f", x$Trial.R2[1,1]),"[", sprintf("%.4f", x$Trial.R2[1,2]),";", sprintf("%.4f", x$Trial.R2[1,3]) , "]", "\n\n")

  cat("Estimated treatment effects on surrogate (pfs_est) and survival (surv_est) endpoint: \n\n")
  print(x$EstTreatEffects[,c(1,2,3,4,5,7)], row.names = FALSE)
}

#' Generates a plot of the estimated treatment effects for the surrogate endpoint versus the estimated treatment effects for the true endpoint for an object fitted with the 'MetaAnalyticSurvSurv()' function.
#'
#' @method plot MetaAnalyticSurvSurv
#'
#' @param x An object of class 'MetaAnalyticSurvSurv' fitted with the 'MetaAnalyticSurvSurv()' function.
#' @param ... ...
#'
#' @return A plot of the type ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal4")
#' fit <- MetaAnalyticSurvSurv(data = colorectal4, true = truend, trueind = trueind, surrog = surrogend,
#'                trt = treatn, center = center, trial = trialend, patientid = patid, adjustment="unadjusted")
#' plot(fit)
#' }
#'
plot.MetaAnalyticSurvSurv <- function(x,...){
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    pfs_est <- surv_est <- sample_size <- NULL
    estimated_treatment_effects <- x$EstTreatEffects

    estimated_treatment_effects$sample_size <- as.numeric(estimated_treatment_effects$sample_size)
    estimated_treatment_effects$surv_est <- as.numeric(estimated_treatment_effects$surv_est)
    estimated_treatment_effects$pfs_est <- as.numeric(estimated_treatment_effects$pfs_est)

    # Create the scatter plot
    ggplot2::ggplot(data = estimated_treatment_effects, ggplot2::aes(x = pfs_est, y = surv_est, size = sample_size)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", se = FALSE, color = "royalblue3") +
      ggplot2::labs(x = "Treatment effect on surrogate", y = "Treatment effect on true") +
      ggplot2::ggtitle("Treatment effect on true endpoint vs. treatment effect on surrogate endpoint") +
      ggplot2::theme(legend.position="none")

  } else {
    stop("ggplot2 is not installed. Please install ggplot2 to use this function.")
  }
}



