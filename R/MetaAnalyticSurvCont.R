

#' Compute surrogacy measures for a continuous (normally-distributed) surrogate and a time-to-event true endpoint in the meta-analytic multiple-trial setting.
#'
#'The function 'MetaAnalyticSurvCont()' fits the model for a continuous surrogate and time-to-event true endpoint described by Alonso et al. (2016) in the meta-analytic multiple-trial setting.
#'
#' @details
#'
#' # Model
#'
#' In the model, a copula-based model is used for the true time-to-event endpoint and the surrogate continuous, normally distributed endpoint.
#' More specifically, three copulas can be used: the Clayton copula, Hougaard copula and Plackett copula. The marginal model for the true endpoint is the proportional hazard model.
#' The marginal model for the surrogate endpoint is the classical linear regression model.
#' The quality of the surrogate at the individual level can be evaluated by either Kendall's \eqn{\tau} or Spearman's \eqn{\rho}, depending on which copula function is used.
#' The quality of the surrogate at the trial level can be evaluated by considering the \eqn{R^2_{trial}} between the estimated treatment effects.
#'
#' # Data Format
#'
#' The data frame must contains the following columns:
#'
#' * a column with the observed time-to-event for the true endpoint
#' * a column with the time-to-event indicator for the true endpoint: 1 if the event is observed, 0 otherwise
#' * a column with the continuous surrogate endpoint
#' * a column with the treatment indicator: 0 or 1
#' * a column with the trial indicator
#' * a column with the center indicator. If there are no different centers within each trial, the center indicator is equal to the trial indicator
#' * a column with the patient indicator
#'
#' @references Alonso A, Bigirumurame T, Burzykowski T, Buyse M, Molenberghs G, Muchene L, Perualila NJ, Shkedy Z, Van der Elst W, et al. (2016). Applied surrogate endpoint evaluation methods with SAS and R. CRC Press New York
#'
#' @param data A data frame with the correct columns (See Data Format).
#' @param true Observed time-to-event for true endpoint.
#' @param trueind Time-to-event indicator for the true endpoint.
#' @param surrog Continuous surrogate endpoint.
#' @param trt Treatment indicator.
#' @param center Center indicator (equal to trial if there are no different centers). This is the unit for which specific treatment effects are estimated.
#' @param trial Trial indicator. This is the unit for which common baselines are to be used.
#' @param patientid Patient indicator.
#' @param copula The copula that is used, either "Clayton", "Hougaard" or "Plackett"
#' @param adjustment The adjustment that should be made for the trial-level surrogacy, either "unadjusted", "weighted" or "adjusted"
#'
#' @return Returns an object of class "MetaAnalyticSurvCont" that can be used to evaluate surrogacy and contains the following elements:
#'
#' * Indiv.Surrogacy: a data frame that contains the measure for the individual level surrogacy and 95% confidence interval.
#' * Trial.R2: a data frame that contains the \eqn{R^2_{trial}} and 95% confidence interval to evaluate surrogacy at the trial level.
#' * EstTreatEffects: a data frame that contains the estimated treatment effects and sample size for each trial.
#' * nlm.output: output of the maximization procedure (nlm) to maximize the likelihood.
#'
#' @export
#'
#' @author Dries De Witte
#'
#' @examples
#' \dontrun{
#' data("prostate")
#' fit <- MetaAnalyticSurvCont(data = prostate, true = SURVTIME, trueind = SURVIND, surrog = PSA,
#' trt = TREAT, center = TRIAL, trial = TRIAL, patientid = PATID,
#' copula = "Hougaard", adjustment = "weighted")
#' summary(fit)
#' print(fit)
#' plot(fit)
#' }
MetaAnalyticSurvCont <- function(data, true, trueind, surrog,
                                 trt, center, trial, patientid, copula, adjustment) {
  centnumb <- NULL
  data_m <- data
  data_m$surv <- data[[substitute(true)]]
  data_m$survind <- data[[substitute(trueind)]]
  data_m$cont <- data[[substitute(surrog)]]
  data_m$treat <- data[[substitute(trt)]]
  data_m$center <- data[[substitute(center)]]
  data_m$trial <- data[[substitute(trial)]]
  data_m$patientid <- data[[substitute(patientid)]]

  data_m <- data_m[, c("patientid", "treat", "center", "survind", "trial", "cont", "surv")]

  data_m <- data_m[order(data_m$center, data_m$trial),]

  data_m$centnumb <- cumsum(!duplicated(data_m$center))

  intercepts <- data.frame(centnumb = integer(), intercept = numeric())
  slopes <- data.frame(centnumb = integer(), treat = numeric())
  fit_statistics <- data.frame(centnumb = integer(), sigma2 = numeric())

  unique_centnumbs <- unique(data_m$centnumb)
  for (i in unique_centnumbs) {
    subset_data <- subset(data_m, centnumb == i)
    model <- lm(cont ~ treat, data = subset_data)

    coef_values <- coef(model)
    intercept <- coef_values["(Intercept)"]
    treat <- coef_values["treat"]

    sigma <- summary(model)$sigma
    sigma2 <- sigma^2

    intercepts <- rbind(intercepts, data.frame(centnumb = i, intercept = intercept))
    slopes <- rbind(slopes, data.frame(centnumb = i, treat = treat))
    fit_statistics <- rbind(fit_statistics, data.frame(centnumb = i, sigma2 = sigma2))
  }

  estS <- merge(intercepts, slopes, by = "centnumb")
  estS <- merge(estS, fit_statistics, by = "centnumb")

  estt2 <- data.frame(centnumb = numeric(),
                      intercept = numeric(),
                      treat_effect = numeric(),
                      scale = numeric(),
                      stringsAsFactors = FALSE)

  for (i in unique_centnumbs) {
    subset_data <- subset(data_m, centnumb == i)
    model <- survival::survreg(Surv(surv, survind) ~ treat, data = subset_data, dist = "weibull")
    intercept <- coef(model)["(Intercept)"]
    treat_effect <- coef(model)["treat"]
    scale_param <- model$scale
    estt2 <- tibble::add_row(estt2, centnumb = i,
                             intercept = intercept,
                             treat_effect = treat_effect,
                             scale = scale_param)
  }

  estT <- estt2

  estT$trt <- -estT$treat_effect / estT$scale
  estT <- estT[, c("centnumb", "intercept", "scale", "trt")]

  MetaAnalyticSurvCont_clayton <- function(data_m) {

    #clayton copula
    output <- 1
    theta <- 2

    # Read variables into vectors and matrices
    surv <- as.vector(data_m$surv)
    survind <- as.vector(data_m$survind)
    cont <- as.vector(data_m$cont)
    treat <- as.matrix(data_m$treat)
    centnumb <- as.vector(data_m$centnumb)
    trial <- as.vector(data_m$trial)

    cents <- unique(centnumb)
    numcents <- length(cents)

    loglik <- function(param) {

      lik <- 0
      theta <- exp(param[1]) + 1

      for (i in 1:numcents) {

        sigma <- exp(param[2+(i-1)*6])
        mu <- param[3+(i-1)*6]
        alpha <- param[4+(i-1)*6]

        lambda <- exp(param[5+(i-1)*6])
        p <- param[6+(i-1)*6]
        beta <- param[7+(i-1)*6]

        t <- as.matrix(surv[centnumb==i])
        delta <- as.matrix(survind[centnumb==i])
        s <- as.matrix(cont[centnumb==i])
        z <- as.matrix(treat[centnumb==i])

        fsurvT <- exp(-lambda*exp(beta*z)*(t^p))
        fdensT <- p*lambda*exp(beta*z)*(t^(p-1))*fsurvT


        fdensS <- dnorm(s, mean = mu + alpha * z, sd = sigma)
        fsurvS <- 1 - pnorm(s, mean = mu + alpha * z, sd = sigma)

        cop <- fsurvS^(1 - theta) + fsurvT^(1 - theta) - 1

        a <- log(theta) + ((2 * theta - 1) / (1 - theta)) * log(cop) - theta * log(fsurvS) - theta * log(fsurvT) + log(fdensS) + log(fdensT)
        b <- (theta / (1 - theta)) * log(cop) - theta * log(fsurvS) + log(fdensS)
        likc <- sum(delta * a + (1 - delta) * b)
        lik <- lik + likc
      }

      return(-lik)
    }

    #initial parameters
    intercS <- estS$intercept
    alpha <- estS$treat
    sigma2 <- estS$sigma2
    centnumb_a <- estS$centnumb

    # Read columns from estT
    intercT <- estT$intercept
    beta <- estT$trt
    scale <- estT$scale
    centnumb_b <- estT$centnumb

    parms <- c()

    for (c in 1:numcents) {
      intT_c <- intercT[centnumb_b == c]
      beta_c <- beta[centnumb_b == c]
      scale_c <- scale[centnumb_b == c]

      p_c <- 1 / scale_c
      llambda_c <- -intT_c / scale_c
      bet_c <- -beta_c / scale_c

      intS_c <- intercS[centnumb_a == c]
      alpha_c <- alpha[centnumb_a == c]
      sigma2_c <- sigma2[centnumb_a == c]

      parms <- c(parms, 0.5 * log(sigma2_c))
      parms <- c(parms, intS_c)
      parms <- c(parms, alpha_c)
      parms <- c(parms, llambda_c)
      parms <- c(parms, p_c)
      parms <- c(parms, bet_c)
    }

    thet0 <- theta;
    lthet0 <- log(thet0-1)

    p0 <- as.vector(c(lthet0, parms))

    suppressWarnings(nlm.output <- nlm(loglik, p0 , hessian=TRUE, iterlim = 1000))
    #return(nlm.output)

    est <- nlm.output$estimate

    hess<-nlm.output$hessian
    fisher_info<-solve(hess)
    var <- fisher_info
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
    cont_est <- c()
    surv_est <- c()
    cont_se <- c()
    surv_se <- c()
    cova <- c()
    weight <- c()
    for (c in 1:numcents) {
      cont_est <- c(cont_est, est[4 + (c - 1) * 6])
      surv_est <- c(surv_est, est[7 + (c - 1) * 6])
      cont_se <- c(cont_se, se[4 + (c - 1) * 6])
      surv_se <- c(surv_se, se[7 + (c - 1) * 6])
      cova <- rbind(cova, var[4 + (c - 1) * 6, 7 + (c - 1) * 6])
      hulp <- centnumb[centnumb==c]
      weight <- c(weight, length(hulp))
    }

    memo <- cbind(unique(data_m$center), surv_est, surv_se, cont_est, cont_se, cova, weight)
    colnames(memo) <- c("center", "surv_est", "surv_se", "cont_est", "cont_se", "cova", "weight")
    memo <- as.data.frame(memo)

    output_list <- list(copula_parameter = alpha_delta_df,
                        Indiv.Surrogacy = tau_delta_df, memo=memo, nlm.output=nlm.output)
    return(output_list)

  }

  MetaAnalyticSurvCont_hougaard <- function(data_m) {

    output <- 1
    theta <- 0.25

    # Read variables into vectors and matrices
    surv <- as.vector(data_m$surv)
    survind <- as.vector(data_m$survind)
    cont <- as.vector(data_m$cont)
    treat <- as.matrix(data_m$treat)
    centnumb <- as.vector(data_m$centnumb)
    trial <- as.vector(data_m$trial)

    cents <- unique(centnumb)
    numcents <- length(cents)

    loglik <- function(param) {

      lik <- 0
      theta <- exp(param[1])

      for (i in 1:numcents) {

        sigma <- exp(param[2+(i-1)*6])
        mu <- param[3+(i-1)*6]
        alpha <- param[4+(i-1)*6]

        lambda <- exp(param[5+(i-1)*6])
        p <- param[6+(i-1)*6]
        beta <- param[7+(i-1)*6]

        t <- as.matrix(surv[centnumb==i])
        delta <- as.matrix(survind[centnumb==i])
        s <- as.matrix(cont[centnumb==i])
        z <- as.matrix(treat[centnumb==i])

        fsurvT <- exp(-lambda*exp(beta*z)*(t^p))
        fdensT <- p*lambda*exp(beta*z)*(t^(p-1))*fsurvT


        fdensS <- dnorm(s, mean = mu + alpha * z, sd = sigma)
        fsurvS <- 1 - pnorm(s, mean = mu + alpha * z, sd = sigma)

        mcop <- (-log(fsurvS))^(1/theta)+(-log(fsurvT))^(1/theta)
        cop <- exp(-mcop^theta)

        a <- log(cop)+log(-log(fsurvS))*(1-theta)/theta-log(fsurvS)+log(-log(fsurvT))*(1-theta)/theta-log(fsurvT)+
          (theta-2)*log(mcop)+log(mcop^theta-(theta-1)/theta)+log(fdensS)+log(fdensT)
        b <- log(cop)+(theta-1)*log(mcop)+log(-log(fsurvS))*(1-theta)/theta-log(fsurvS)+log(fdensS)

        likc <- sum(delta * a + (1 - delta) * b)
        lik <- lik + likc
      }

      return(-lik)
    }

    #initial parameters
    intercS <- estS$intercept
    alpha <- estS$treat
    sigma2 <- estS$sigma2
    centnumb_a <- estS$centnumb

    # Read columns from estT
    intercT <- estT$intercept
    beta <- estT$trt
    scale <- estT$scale
    centnumb_b <- estT$centnumb

    parms <- c()

    for (c in 1:numcents) {
      intT_c <- intercT[centnumb_b == c]
      beta_c <- beta[centnumb_b == c]
      scale_c <- scale[centnumb_b == c]

      p_c <- 1 / scale_c
      llambda_c <- -intT_c / scale_c
      bet_c <- -beta_c / scale_c

      intS_c <- intercS[centnumb_a == c]
      alpha_c <- alpha[centnumb_a == c]
      sigma2_c <- sigma2[centnumb_a == c]

      parms <- c(parms, 0.5 * log(sigma2_c))
      parms <- c(parms, intS_c)
      parms <- c(parms, alpha_c)
      parms <- c(parms, llambda_c)
      parms <- c(parms, p_c)
      parms <- c(parms, bet_c)
    }

    thet0 <- theta;
    lthet0 <- log(thet0)

    p0 <- as.vector(c(lthet0, parms))

    suppressWarnings(nlm.output <- nlm(loglik, p0 , hessian=TRUE, iterlim = 1000))
    #return(nlm.output)

    est <- nlm.output$estimate

    hess<-nlm.output$hessian
    fisher_info<-solve(hess)
    var <- fisher_info
    se <-sqrt(diag(fisher_info))

    lo <- est - qnorm(0.975) * se
    up <- est + qnorm(0.975) * se

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
    cont_est <- c()
    surv_est <- c()
    cont_se <- c()
    surv_se <- c()
    cova <- c()
    weight <- c()
    for (c in 1:numcents) {
      cont_est <- c(cont_est, est[4 + (c - 1) * 6])
      surv_est <- c(surv_est, est[7 + (c - 1) * 6])
      cont_se <- c(cont_se, se[4 + (c - 1) * 6])
      surv_se <- c(surv_se, se[7 + (c - 1) * 6])
      cova <- rbind(cova, var[4 + (c - 1) * 6, 7 + (c - 1) * 6])
      hulp <- centnumb[centnumb==c]
      weight <- c(weight, length(hulp))
    }

    memo <- cbind(unique(data_m$center), surv_est, surv_se, cont_est, cont_se, cova, weight)
    colnames(memo) <- c("center", "surv_est", "surv_se", "cont_est", "cont_se", "cova", "weight")
    memo <- as.data.frame(memo)

    output_list <- list(copula_parameter = alpha_delta_df,
                        Indiv.Surrogacy = tau_delta_df, memo=memo, nlm.output=nlm.output)
    return(output_list)

  }

  MetaAnalyticSurvCont_plackett <- function(data_m) {
    theta <- 2
    output <- 1

    # Read variables into vectors and matrices
    surv <- as.vector(data_m$surv)
    survind <- as.vector(data_m$survind)
    cont <- as.vector(data_m$cont)
    treat <- as.matrix(data_m$treat)
    centnumb <- as.vector(data_m$centnumb)
    trial <- as.vector(data_m$trial)

    cents <- unique(centnumb)
    numcents <- length(cents)

    loglik <- function(param) {

      lik <- 0
      theta <- param[1]

      for (i in 1:numcents) {

        sigma <- exp(param[2+(i-1)*6])
        mu <- param[3+(i-1)*6]
        alpha <- param[4+(i-1)*6]

        lambda <- exp(param[5+(i-1)*6])
        p <- param[6+(i-1)*6]
        beta <- param[7+(i-1)*6]

        t <- as.matrix(surv[centnumb==i])
        delta <- as.matrix(survind[centnumb==i])
        s <- as.matrix(cont[centnumb==i])
        z <- as.matrix(treat[centnumb==i])

        fsurvT <- exp(-lambda*exp(beta*z)*(t^p))
        fdensT <- p*lambda*exp(beta*z)*(t^(p-1))*fsurvT


        fdensS <- dnorm(s, mean = mu + alpha * z, sd = sigma)
        fsurvS <- 1 - pnorm(s, mean = mu + alpha * z, sd = sigma)

        u <- fsurvS
        v <- fsurvT

        contr <- 1+(theta-1)*(u+v)
        hcop <- sqrt(contr^2+4*theta*(1-theta)*u*v)
        cop <- (contr-hcop)/(2*(theta-1))

        dh_u <- (theta - 1) * (contr - 2 * theta * v) / hcop
        dh_v <- (theta - 1) * (contr - 2 * theta * u) / hcop
        d2h_vu <- -(theta - 1) * ((theta + 1) * hcop + (contr - 2 * theta * v) * dh_v) / (hcop^2)

        dc_u <- 0.5 - dh_u / (2 * (theta - 1))
        d2c_vu <- -d2h_vu / (2 * (theta - 1))

        a <- log(d2c_vu) + log(fdensS) + log(fdensT)
        b <- log(dc_u) + log(fdensS)

        likc <- sum(delta * a + (1 - delta) * b)
        lik <- lik + likc
      }

      return(-lik)
    }

    #initial parameters
    intercS <- estS$intercept
    alpha <- estS$treat
    sigma2 <- estS$sigma2
    centnumb_a <- estS$centnumb

    # Read columns from estT
    intercT <- estT$intercept
    beta <- estT$trt
    scale <- estT$scale
    centnumb_b <- estT$centnumb

    parms <- c()

    for (c in 1:numcents) {
      intT_c <- intercT[centnumb_b == c]
      beta_c <- beta[centnumb_b == c]
      scale_c <- scale[centnumb_b == c]

      p_c <- 1 / scale_c
      llambda_c <- -intT_c / scale_c
      bet_c <- -beta_c / scale_c

      intS_c <- intercS[centnumb_a == c]
      alpha_c <- alpha[centnumb_a == c]
      sigma2_c <- sigma2[centnumb_a == c]

      parms <- c(parms, 0.5 * log(sigma2_c))
      parms <- c(parms, intS_c)
      parms <- c(parms, alpha_c)
      parms <- c(parms, llambda_c)
      parms <- c(parms, p_c)
      parms <- c(parms, bet_c)
    }

    thet0 <- theta
    p0 <- as.vector(c(thet0, parms))

    suppressWarnings(nlm.output <- nlm(loglik, p0 , hessian=TRUE, iterlim = 1000))

    est <- nlm.output$estimate

    hess<-nlm.output$hessian
    fisher_info<-solve(hess)
    var <- fisher_info
    se <-sqrt(diag(fisher_info))

    lo <- est - qnorm(0.975) * se
    up <- est + qnorm(0.975) * se

    theta <- est[1]
    se_th <- se[1]
    lo_th <- theta - qnorm(0.975) * se_th
    up_th <- theta + qnorm(0.975) * se_th

    theta_df <- data.frame(cbind(theta, lo_th, up_th), stringsAsFactors = TRUE)
    colnames(theta_df) <- c("Copula parameter theta", "CI lower limit",
                            "CI upper limit")
    rownames(theta_df) <- c(" ")


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

    #### TRIAL LEVEL SURROGACY ####
    cont_est <- c()
    surv_est <- c()
    cont_se <- c()
    surv_se <- c()
    cova <- c()
    weight <- c()
    for (c in 1:numcents) {
      cont_est <- c(cont_est, est[4 + (c - 1) * 6])
      surv_est <- c(surv_est, est[7 + (c - 1) * 6])
      cont_se <- c(cont_se, se[4 + (c - 1) * 6])
      surv_se <- c(surv_se, se[7 + (c - 1) * 6])
      cova <- rbind(cova, var[4 + (c - 1) * 6, 7 + (c - 1) * 6])
      hulp <- centnumb[centnumb==c]
      weight <- c(weight, length(hulp))
    }

    memo <- cbind(unique(data_m$center), surv_est, surv_se, cont_est, cont_se, cova, weight)
    colnames(memo) <- c("center", "surv_est", "surv_se", "cont_est", "cont_se", "cova", "weight")
    memo <- as.data.frame(memo)

    output_list <- list(copula_parameter = theta_df,
                        Indiv.Surrogacy = rho_df, memo=memo, nlm.output=nlm.output)
    return(output_list)

  }

  #### TRIAL LEVEL SURROGACY ####
  trial_level_secondstage_unadjusted <- function(memo) {
    surv_est <- weight <- cont_est <- NULL
    surv_eff <- subset(memo, select = c(center, surv_est, weight))
    colnames(surv_eff)[2] <- "effect"
    surv_eff$endp <- "MAIN"

    cont_eff <- subset(memo, select = c(center, cont_est, weight))
    colnames(cont_eff)[2] <- "effect"
    cont_eff$endp <- "SURR"

    # Combine surv_eff and cont_eff data frames
    shihco <- rbind(surv_eff, cont_eff)
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
    weight <- NULL
    model <- lm(surv_est ~ cont_est,
                data =memo, weights = weight)

    R2 <- summary(model)$r.squared

    ### based on cortinas
    Trial.ID <- memo$center
    N <- length(Trial.ID)
    R2.sd <- sqrt((4 * R2 * (1 - R2)^2)/(N - 3))
    lo_R2 <- max(0, R2 + qnorm(0.05/2) * (R2.sd))
    up_R2 <- min(1, R2 + qnorm(1 - 0.05/2) * (R2.sd))

    Trial.R2 <- data.frame(cbind(R2, lo_R2, up_R2), stringsAsFactors = TRUE)
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
      se2_sq <- memo$cont_se[i]^2
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
    stage1 <- MetaAnalyticSurvCont_clayton(data_m = data_m)

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
    stage1 <- MetaAnalyticSurvCont_hougaard(data_m = data_m)

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
    stage1 <- MetaAnalyticSurvCont_plackett(data_m = data_m)

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

  memo <- stage1$memo
  colnames(memo) <- c("center", "surv_est", "surv_se", "surro_est", "surro_se", "cova", "sample_size")

  output_list <- list(Indiv.Surrogacy = stage1$Indiv.Surrogacy,
                      Trial.R2 = stage2$Trial.R2,
                      EstTreatEffects = memo,
                      nlm.output = stage1$nlm.output,
                      Call = match.call())

  class(output_list) <- "MetaAnalyticSurvSurv"

  return(output_list)

}

#' Provides a summary of the surrogacy measures for an object fitted with the 'MetaAnalyticSurvCont()' function.
#'
#' @method summary MetaAnalyticSurvCont
#'
#' @param object An object of class 'MetaAnalyticSurvCont' fitted with the 'MetaAnalyticSurvCont()' function.
#' @param ... ...
#'
#' @return The surrogacy measures with their 95% confidence intervals.
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal")
#' data("prostate")
#' fit <- MetaAnalyticSurvCont(data = prostate, true = SURVTIME, trueind = SURVIND, surrog = PSA,
#' trt = TREAT, center = TRIAL, trial = TRIAL, patientid = PATID,
#' copula = "Hougaard", adjustment = "weighted")
#' summary(fit)
#' }

summary.MetaAnalyticSurvCont <- function(object,...){
  indiv.name <- colnames(object$Indiv.Surrogacy[1])
  cat("Surrogacy measures with 95% confidence interval \n\n")
  cat("Individual level surrogacy: ", "\n\n")
  cat(sprintf(indiv.name), ": ", sprintf("%.4f", object$Indiv.Surrogacy[1,1]), "[", sprintf("%.4f", object$Indiv.Surrogacy[1,2]),";", sprintf("%.4f", object$Indiv.Surrogacy[1,3]) , "]", "\n\n")
  cat("Trial level surrogacy: ", "\n\n")
  cat("R Square: ", sprintf("%.4f", object$Trial.R2[1,1]),"[", sprintf("%.4f", object$Trial.R2[1,2]),";", sprintf("%.4f", object$Trial.R2[1,3]) , "]", "\n\n")
}

#' Prints all the elements of an object fitted with the 'MetaAnalyticSurvCont()' function.
#'
#' @method print MetaAnalyticSurvCont
#'
#' @param x An object of class 'MetaAnalyticSurvCont' fitted with the 'MetaAnalyticSurvCont()' function.
#' @param ... ...
#'
#' @return The surrogacy measures with their 95% confidence intervals and the estimated treatment effect on the surrogate and true endpoint.
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal4")
#' data("prostate")
#' fit <- MetaAnalyticSurvCont(data = prostate, true = SURVTIME, trueind = SURVIND, surrog = PSA,
#' trt = TREAT, center = TRIAL, trial = TRIAL, patientid = PATID,
#' copula = "Hougaard", adjustment = "weighted")
#' print(fit)
#' }
print.MetaAnalyticSurvCont <- function(x,...){
  indiv.name <- colnames(x$Indiv.Surrogacy[1])
  cat("Surrogacy measures with 95% confidence interval \n\n")
  cat("Individual level surrogacy: ", "\n\n")
  cat(sprintf(indiv.name), ": ", sprintf("%.4f", x$Indiv.Surrogacy[1,1]), "[", sprintf("%.4f", x$Indiv.Surrogacy[1,2]),";", sprintf("%.4f", x$Indiv.Surrogacy[1,3]) , "]", "\n\n")
  cat("Trial level surrogacy: ", "\n\n")
  cat("R Square: ", sprintf("%.4f", x$Trial.R2[1,1]),"[", sprintf("%.4f", x$Trial.R2[1,2]),";", sprintf("%.4f", x$Trial.R2[1,3]) , "]", "\n\n")

  cat("Estimated treatment effects on surrogate (surro_est) and survival (surv_est) endpoint: \n\n")
  print(x$EstTreatEffects[,c(1,2,3,4,5,7)], row.names = FALSE)
}

#' Generates a plot of the estimated treatment effects for the surrogate endpoint versus the estimated treatment effects for the true endpoint for an object fitted with the 'MetaAnalyticSurvCont()' function.
#'
#' @method plot MetaAnalyticSurvCont
#'
#' @param x An object of class 'MetaAnalyticSurvCont' fitted with the 'MetaAnalyticSurvCont()' function.
#' @param ... ...
#'
#' @return A plot of the type ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal4")
#' data("prostate")
#' fit <- MetaAnalyticSurvCont(data = prostate, true = SURVTIME, trueind = SURVIND, surrog = PSA,
#' trt = TREAT, center = TRIAL, trial = TRIAL, patientid = PATID,
#' copula = "Hougaard", adjustment = "weighted")
#' plot(fit)
#' }
#'
plot.MetaAnalyticSurvCont <- function(x,...){
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    surro_est <- surv_est <- sample_size <- NULL
    estimated_treatment_effects <- x$EstTreatEffects

    estimated_treatment_effects$sample_size <- as.numeric(estimated_treatment_effects$sample_size)
    estimated_treatment_effects$surv_est <- as.numeric(estimated_treatment_effects$surv_est)
    estimated_treatment_effects$surro_est <- as.numeric(estimated_treatment_effects$surro_est)

    # Create the scatter plot
    suppressWarnings({p <- ggplot2::ggplot(data = estimated_treatment_effects, ggplot2::aes(x = surro_est, y = surv_est, size = sample_size)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", se = FALSE, color = "royalblue3") +
      ggplot2::labs(x = "Treatment effect on surrogate", y = "Treatment effect on true") +
      ggplot2::ggtitle("Treatment effect on true endpoint vs. treatment effect on surrogate endpoint") +
      ggplot2::theme(legend.position="none")
    })
    suppressWarnings(print(p))

  } else {
    stop("ggplot2 is not installed. Please install ggplot2 to use this function.")
  }
}



