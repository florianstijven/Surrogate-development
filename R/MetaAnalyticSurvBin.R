

#' Compute surrogacy measures for a binary surrogate and a time-to-event true endpoint in the meta-analytic multiple-trial setting.
#'
#'The function 'MetaAnalyticSurvBin()' fits the model for a binary surrogate and time-to-event true endpoint developed by Burzykowski et al. (2004) in the meta-analytic multiple-trial setting.
#'
#' @details
#'
#' # Model
#'
#' In the model developed by Burzykowski et al. (2004), a copula-based model is used for the true endpoint and a latent continuous variable, underlying the surrogate endpoint.
#' More specifically, the Plackett copula is used. The marginal model for the surrogate endpoint is a logistic regression model. For the true endpoint, the proportional hazard model is used.
#' The quality of the surrogate at the individual level can be evaluated by using the copula parameter Theta, which takes the form of a global odds ratio.
#' The quality of the surrogate at the trial level can be evaluated by considering the correlation coefficient between the estimated treatment effects.
#'
#' # Data Format
#'
#' The data frame must contains the following columns:
#'
#' * a column with the observed time-to-event (true endpoint)
#' * a column with the time-to-event indicator: 1 if true event is observed, 0 otherwise
#' * a column with the binary surrogate endpoint: 1 or 2
#' * a column with the treatment indicator: 0 or 1
#' * a column with the trial indicator
#' * a column with the center indicator. If there are no different centers within each trial, the center indicator can be equal to the trial indicator
#' * a column with the patient indicator
#'
#' @references Burzykowski, T., Molenberghs, G., & Buyse, M. (2004). The validation of surrogate end points by using data from randomized clinical trials: a case-study in advanced colorectal cancer. Journal of the Royal Statistical Society Series A: Statistics in Society, 167(1), 103-124.
#'
#' @param data A data frame with the correct columns (See details).
#' @param true Observed time-to-event (true endpoint).
#' @param trueind Time-to-event indicator.
#' @param surrog Binary surrogate endpoint, coded as 1 or 2.
#' @param trt Treatment indicator, coded as 0 or 1.
#' @param center Center indicator (equal to trial if there are no different centers). This is the unit for which specific treatment effects are estimated.
#' @param trial Trial indicator. This is the unit for which common baselines are to be used.
#' @param patientid Patient indicator.
#' @param adjustment The adjustment that should be made for the trial-level surrogacy, either "unadjusted", "weighted" or "adjusted"
#'
#' @return Returns an object of class "MetaAnalyticSurvBin" that can be used to evaluate surrogacy and contains the following elements:
#'
#' * Indiv.GlobalOdds: a data frame that contains the global odds ratio and 95% confidence interval to evaluate surrogacy at the individual level.
#' * Trial.R2: a data frame that contains the correlation coefficient and 95% confidence interval to evaluate surrogacy at the trial level.
#' * EstTreatEffects: a data frame that contains the estimated treatment effects and sample size for each trial.
#' * fit_output: output of the maximization procedure (nlm) to maximize the likelihood function.
#'
#' @export
#'
#' @author Dries De Witte
#'
#' @examples
#' \dontrun{
#' data("colorectal")
#' fit_bin <- MetaAnalyticSurvBin(data = colorectal, true = surv, trueind = SURVIND, surrog = responder,
#'                    trt = TREAT, center = CENTER, trial = TRIAL, patientid = patientid, adjustment="unadjusted")
#' print(fit_bin)
#' summary(fit_bin)
#' plot(fit_bin)
#' }
#'
MetaAnalyticSurvBin <- function(data, true, trueind, surrog,
                    trt, center, trial, patientid, adjustment) {
  resp_est <- surv_est <- sample_size <- Center_ID <- p <- shape <- variable <- NULL
  dataset1xxy <- data
  dataset1xxy$surv <- data[[substitute(true)]]
  dataset1xxy$survind <- data[[substitute(trueind)]]
  dataset1xxy$suro <- data[[substitute(surrog)]]
  dataset1xxy$Z <- data[[substitute(trt)]]
  dataset1xxy$Center_ID <- data[[substitute(center)]]
  dataset1xxy$Trial_ID <- data[[substitute(trial)]]
  dataset1xxy$Pat_ID <- data[[substitute(patientid)]]

  assign_data <- function(dataset1xxy) {
    dataset1xxy <- dataset1xxy[order(dataset1xxy$Center_ID), ]
    centnum <- unique(dataset1xxy$Center_ID)
    centid <- list()

    for (i in seq_along(centnum)) {
      centid[[paste0("centid", sprintf("%02d", i))]] <- centnum[i]
    }

    c <- sprintf("%02d", length(centnum))

    for (i in 1:length(centnum)) {
      const_var <- paste0("const", sprintf("%02d", i))
      treat_var <- paste0("treat", sprintf("%02d", i))

      dataset1xxy[[const_var]] <- as.integer(dataset1xxy$Center_ID == centid[[paste0("centid", sprintf("%02d", i))]])
      dataset1xxy[[treat_var]] <- as.numeric(dataset1xxy[[const_var]]) * as.numeric(dataset1xxy[["Z"]])
    }

    return(dataset1xxy)
  }

  dataset1xxy <- assign_data(dataset1xxy)

  c <- length(unique(dataset1xxy$Center_ID))
  covariate_names <- paste0("treat", sprintf("%02d", 1:c))

  dataset1xxy$suro2 <- dataset1xxy$suro-1
  dataset1xxy$suro2_bis <- factor(dataset1xxy$suro2, levels = c("1", "0"))
  formula <- as.formula(paste("suro2_bis ~", paste(covariate_names, collapse = " + ")))

  model <- stats::glm(formula, data = dataset1xxy,
               family = binomial(link = "logit"))
  ests <- coef(model)
  estimS <- as.matrix(ests)

  t <- as.matrix(dataset1xxy$surv)
  delta <- as.vector(dataset1xxy$survind)
  s <- as.vector(dataset1xxy$suro)
  x <- as.matrix(dataset1xxy$Z)
  center <- as.vector(dataset1xxy$Center_ID)
  trial <- as.vector(dataset1xxy$Trial_ID)

  qq <- ncol(x)
  qq
  r <- 0
  K <- length(unique(s))

  cents <- unique(center)
  numcents <- length(cents)
  tri <- unique(trial)
  numtri <- length(tri)

  par0s <- estimS

  intercS <- par0s[1:(K - 1), ]
  efftrtS <- par0s[(K - 1 + r + 1):(K - 1 + r + numcents * qq), ]

  if (numtri > 1) {
    hulp <- matrix(0, nrow = numtri - 1, ncol = 1)
    hulp <- as.vector(hulp)
    param0s <- c(intercS, hulp, efftrtS)
  } else {
    param0s <- c(intercS, efftrtS)
  }

  estimate <- as.data.frame(param0s)
  dataset1xxy <- dataset1xxy[order(dataset1xxy$Center_ID),]
  model_results <- list()

  estt2 <- data.frame(center = numeric(),
                      intercept = numeric(),
                      treat_effect = numeric(),
                      scale = numeric(),
                      stringsAsFactors = FALSE)

  centers <- unique(dataset1xxy$Center_ID)

  for (i in seq_along(centers)) {
    center_data <- subset(dataset1xxy, Center_ID == centers[i])
    model <- survival::survreg(Surv(surv, survind) ~ Z, data = center_data, dist = "weibull")

    intercept <- coef(model)["(Intercept)"]
    treat_effect <- coef(model)["Z"]
    scale_param <- model$scale

    estt2 <- estt2 %>%
      tibble::add_row(center = as.numeric(centers[i]),
              intercept = intercept,
              treat_effect = treat_effect,
              scale = scale_param)
  }
  estt <- estt2
  estt$p <- 1 / estt$scale
  estt$shape <- -estt$intercept
  estt$trt <- -estt$treat_effect / estt$scale
  estt <- estt[, c("center", "p", "shape", "trt")]

  estimT <- estt %>%
    tidyr::pivot_longer(cols = c(p, shape, trt),
                 names_to = "variable",
                 values_to = "value") %>%
    dplyr::mutate(variable = dplyr::case_when(
      variable == "p" ~ "p",
      variable == "shape" ~ "shape",
      variable == "trt" ~ "trt"
    ))

  estimT <- estimT %>%
    arrange(variable, center)
  estimT <- estimT[, "value", drop = FALSE]
  estim <- estimT

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

  d2quth(2, 3, 0.5)

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

  d2Cuv <- function(u, v, theta) {
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

  cFt <- function(lambda, p, beta, t, s, zt, zs) {
    lt <- lambda * t
    Ft <- 1 - exp(-(lt^p) * exp(zt %*% beta))
    return(Ft)
  }

  pft <- function(lambda, p, beta, t, s, zt, zs) {
    lt <- lambda * t
    f <- lambda * p * (lt^(p - 1)) * exp(zt %*% beta) * exp(-(lt^p) * exp(zt %*% beta))
    return(f)
  }

  cFs <- function(alpha, t, s, zt, zs) {
    Fs <- exp(zs %*% alpha) / (1 + exp(zs %*% alpha))
    return(Fs)
  }

  xsect <- function(...) {
    args <- list(...)
    result <- Reduce(intersect, args)
    result <- sort(result)
    return(result)
  }

  design <- function(x) {
    unique_vals <- sort(unique(x))
    design_matrix <- matrix(0, nrow = length(x), ncol = length(unique_vals))

    for (i in 1:length(unique_vals)) {
      design_matrix[, i] <- as.numeric(x == unique_vals[i])
    }

    colnames(design_matrix) <- unique_vals
    return(design_matrix)
  }

  pgts <- function(theta, lambda, p, beta, alpha, t, s, zt, zs, x) {
    pft <- pft(lambda, p, beta, t, s, zt, zs)
    cFt <- cFt(lambda, p, beta, t, s, zt, zs)

    K <- length(unique(s))

    dCu <- NULL
    cFs <- NULL

    for (i in 1:(K-1)) {
      zero <- matrix(0, nrow = nrow(as.matrix(x)), ncol = K-1)

      zero[, i] <- 1

      zsi <- cbind(zero, matrix(1, nrow = nrow(as.matrix(x)), ncol = 1), x)

      cFsi <- cFs(alpha, t, s, zt, zsi)

      cFs <- cbind(cFs, cFsi)

      dCu <- cbind(dCu, d1Cu(cFt, cFsi, theta))

    }

    dCui <- cbind(dCu, matrix(1, nrow = nrow(dCu), ncol = 1))
    dCui_1 <- cbind(matrix(0, nrow = nrow(as.matrix(x)), ncol = 1), dCu)

    yresp <- xsect(1:K,unique(s))
    des <- matrix(0, nrow = length(s), ncol = K)

    if (length(yresp) < K) {
      des[cbind(1:length(s), match(s, yresp))] <- 1
    } else {
      des <- design(s)
    }

    g <- as.matrix(rowSums(des*((dCui - dCui_1) * as.vector(pft))))

    return(g)
  }

  cGts <- function(theta, lambda, p, beta, alpha, t, s, zt, zs, x) {
    cFt <- cFt(lambda, p, beta, t, s, zt, zs)
    K <- length(unique(s))

    cFs <- NULL
    Fts <- NULL

    for (i in 1:(K-1)) {
      zero <- matrix(0, nrow = nrow(as.matrix(x)), ncol = K-1)
      zero[, i] <- 1
      zsi <- cbind(zero, matrix(1, nrow = nrow(as.matrix(x)), ncol = 1), x)
      cFsi <- cFs(alpha, t, s, zt, zsi)
      cFs <- cbind(cFs, cFsi)
      Fts <- cbind(Fts, C(cFt, cFsi, theta))
    }

    diffs1 <- cbind(cFs, matrix(1, nrow = nrow(as.matrix(cFs)), ncol = 1))-cbind(matrix(0, nrow = nrow(as.matrix(cFs)), ncol = 1),cFs)
    diffs2 <- cbind(Fts, cFt) - cbind(matrix(0, nrow(as.matrix(Fts)), 1), Fts)

    yresp <- xsect(1:K,unique(s))
    des <- matrix(0, nrow = length(s), ncol = K)
    if (length(yresp) < K) {
      des[cbind(1:length(s), match(s, yresp))] <- 1
    } else {
      des <- design(s)
    }

    g <- as.matrix(rowSums(des * (diffs1 - diffs2)))

    #g=(des#(diffs1-diffs2))[,+];

    return(g)
  }

  loglik <- function(param) {
    t <- as.matrix(dataset1xxy$surv)
    delta <- as.vector(dataset1xxy$survind)
    s <- as.vector(dataset1xxy$suro)
    x <- as.matrix(dataset1xxy$Z)
    center <- as.vector(dataset1xxy$Center_ID)
    trial <- as.vector(dataset1xxy$Trial_ID)

    cents <- unique(center)
    numcents <- length(cents)
    tri <- unique(trial)
    numtri <- length(tri)
    K <- length(unique(s))

    zt <- x
    hulp <- design(s)
    des <- hulp[,1:(K-1)]
    zs <- as.matrix(cbind(des, matrix(1, length(x), 1), x))

    theta <- param[1]
    eta <- param[(1 + 2 * numtri + numcents + 1):(1 + 2 * numtri + numcents + K - 1)]
    lik <- 0
    for (i in 1:numcents) {
      ti <- t[center == cents[i]]
      si <- s[center == cents[i]]
      deltai <- delta[center == cents[i]]
      zti <- as.matrix(zt[center == cents[i], ])
      zsi <- as.matrix(zs[center == cents[i], ])
      xi <- x[center == cents[i], ]
      hulp <- as.numeric(unique(trial[center == cents[i]]))
      stud <- which(tri == hulp)
      pi <- param[1 + stud]
      lli <- param[1 + numtri + stud]
      efftrtT <- param[1 + 2 * numtri + i]
      betai <- as.matrix(cbind(efftrtT))
      if (stud == 1) {
        interc <- 0
      } else {
        interc <- param[1 + 2 * numtri + numcents + (K - 1) + stud - 1]
      }
      efftrtS <- param[1 + 2 * numtri + numcents + (K - 1) + numtri - 1 + i]
      alphai <- as.matrix(c(eta, interc, efftrtS))

      lambdai <- exp(lli)
      lti <- lambdai * t[center == cents[i], ]
      pgi <- pgts(theta, lambdai, pi, betai, alphai, ti, si, zti, zsi, xi)
      cGi <- cGts(theta, lambdai, pi, betai, alphai, ti, si, zti, zsi, xi)

      lik <- lik + sum(deltai * log(pgi) + (1 - deltai) * log(cGi))
    }

    return(-lik)
  }

  initparm <- function(theta0) {
    s <- as.vector(dataset1xxy$suro)
    center <- as.vector(dataset1xxy$Center_ID)

    center <- unique(center)
    cents <- unique(center)
    numcents <- length(cents)
    K <- length(unique(s))

    param1 <- estim
    param2 <- estimate

    param0 <- cbind(theta0, t(param1), t(param2))

    return(param0)
  }

  initial_parameters <- initparm(5)

  par0 <- as.vector(initial_parameters)
  suppressWarnings(nlm_output <- nlm(loglik, par0 , hessian=TRUE, iterlim = 10000))

  hessian_m <-nlm_output$hessian
  fisher_info<-solve(hessian_m)
  se <-sqrt(diag(fisher_info))
  est <- nlm_output$estimate
  coef_se <- cbind(est, se)

  endp1 <- coef_se[(2 * numtri + 1 + 1):(2 * numtri + 1 + numcents), ]
  endp2 <- coef_se[(1 + 3 * numtri + numcents + (K - 1) - 1 + 1):nrow(coef_se), ]
  weight <- matrix(0, nrow = numcents, ncol = 1)
  for (i in 1:numcents) {
    weight[i, ] <- nrow(as.matrix(which(center == cents[i])))
  }
  cova <- diag(fisher_info[(2 * numtri + 1 + 1):(2 * numtri + 1 + numcents),
                           (1 + 3 * numtri + numcents + (K - 1) - 1 + 1):length(est)])
  memo <- cbind(cents, endp1, endp2, cova, weight)
  colnames(memo) <- c("center", "surv_est", "surv_se", "resp_est",
                      "resp_se", "cova", "weight")
  memo <- as.data.frame(memo)

  #Individual level surrogacy
  lo <- est - qnorm(0.975) * se
  up <- est + qnorm(0.975) * se
  theta <- est[1]
  se_th <- se[1]
  lo_th <- lo[1]
  up_th <- up[1]

  #Trial level surrogacy
  if(adjustment == "adjusted"){
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

    memo <- as.data.frame(apply(memo, 2, as.numeric))

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
      se2_sq <- memo$resp_se[i]^2
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
  }

  if(adjustment == "weighted"){
    memo <- as.data.frame(apply(memo, 2, as.numeric))

    model <- lm(surv_est ~ resp_est,
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
  }

  if(adjustment == "unadjusted"){
    memo <- as.data.frame(apply(memo, 2, as.numeric))

    surv_eff <- subset(memo, select = c(center, surv_est))
    colnames(surv_eff)[2] <- "effect"
    surv_eff$endp <- "MAIN"
    resp_eff <- subset(memo, select = c(center, resp_est))
    colnames(resp_eff)[2] <- "effect"
    resp_eff$endp <- "SURR"
    shihco <- rbind(surv_eff, resp_eff)
    shihco$center <- as.numeric(shihco$center)
    shihco <- shihco[order(shihco$center, shihco$endp), ]
    shihco$endp <- as.factor(shihco$endp)
    shihco$effect <- as.numeric(shihco$effect)
    invisible(capture.output(model <- nlme::gls(effect ~ -1 +
                                                  factor(endp), data = shihco, correlation = nlme::corCompSymm(form = ~1 |
                                                                                                                 center), weights = nlme::varIdent(form = ~1 | endp),
                                                method = "ML", control = nlme::glsControl(maxIter = 25,
                                                                                          msVerbose = TRUE))))
    rsquared <- nlme::intervals(model, which = "var-cov")$corStruct^2
    R2 <- as.vector(rsquared)[2]
    lo_R2 <- as.vector(rsquared)[1]
    up_R2 <- as.vector(rsquared)[3]
    Trial.R2 <- data.frame(cbind(R2, lo_R2, up_R2), stringsAsFactors = TRUE)
    colnames(Trial.R2) <- c("R2 Trial (unadjusted)", "CI lower limit", "CI upper limit")
    rownames(Trial.R2) <- c(" ")
  }

  Indiv.GlobalOdds <- data.frame(cbind(theta, lo_th, up_th),
                                 stringsAsFactors = TRUE)
  colnames(Indiv.GlobalOdds) <- c("Global Odds Ratio", "CI lower limit",
                                  "CI upper limit")
  rownames(Indiv.GlobalOdds) <- c(" ")
  EstTreatEffects <- memo
  rownames(EstTreatEffects) <- NULL
  colnames(EstTreatEffects) <- c("trial", "surv_est", "surv_se",
                                 "resp_est", "resp_se", "cova", "sample_size")
  output_list <- list(Indiv.GlobalOdds = Indiv.GlobalOdds,
                      Trial.R2 = Trial.R2, EstTreatEffects = EstTreatEffects, fit_output = nlm_output,
                      Call = match.call())
  class(output_list) <- "MetaAnalyticSurvBin"
  return(output_list)

}

#' Provides a summary of the surrogacy measures for an object fitted with the 'MetaAnalyticSurvBin()' function.
#'
#' @method summary MetaAnalyticSurvBin
#'
#' @param object An object of class 'MetaAnalyticSurvBin' fitted with the 'MetaAnalyticSurvBin()' function.
#' @param ... ...
#'
#' @return The surrogacy measures with their 95% confidence intervals.
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal")
#' fit <- MetaAnalyticSurvBin(data = colorectal, true = surv, trueind = SURVIND, surrog = responder,
#'                    trt = TREAT, center = CENTER, trial = TRIAL, patientid = patientid, adjustment="unadjusted")
#' summary(fit)
#' }

summary.MetaAnalyticSurvBin <- function(object,...){
  cat("Surrogacy measures with 95% confidence interval \n\n")
  cat("Individual level surrogacy: ", "\n\n")
  cat("Global Odds Ratio: ", sprintf("%.4f", object$Indiv.GlobalOdds[1,1]), "[", sprintf("%.4f", object$Indiv.GlobalOdds[1,2]),";", sprintf("%.4f", object$Indiv.GlobalOdds[1,3]) , "]", "\n\n")
  cat("Trial level surrogacy: ", "\n\n")
  cat("R Square: ", sprintf("%.4f", object$Trial.R2[1,1]),"[", sprintf("%.4f", object$Trial.R2[1,2]),";", sprintf("%.4f", object$Trial.R2[1,3]) , "]", "\n\n")
}

#' Prints all the elements of an object fitted with the 'MetaAnalyticSurvBin()' function.
#'
#' @method print MetaAnalyticSurvBin
#'
#' @param x An object of class 'MetaAnalyticSurvBin' fitted with the 'MetaAnalyticSurvBin()' function.
#' @param ... ...
#'
#' @return The surrogacy measures with their 95% confidence intervals and the estimated treament effect on the surrogate and true endpoint.
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal")
#' fit_bin <- MetaAnalyticSurvBin(data = colorectal, true = surv, trueind = SURVIND, surrog = responder,
#'                    trt = TREAT, center = CENTER, trial = TRIAL, patientid = patientid, adjustment="unadjusted")
#' print(fit_bin)
#' }
print.MetaAnalyticSurvBin <- function(x,...){
  cat("Surrogacy measures with 95% confidence interval \n\n")
  cat("Individual level surrogacy: ", "\n\n")
  cat("Global Odds Ratio: ", sprintf("%.4f", x$Indiv.GlobalOdds[1,1]), "[", sprintf("%.4f", x$Indiv.GlobalOdds[1,2]),";", sprintf("%.4f", x$Indiv.GlobalOdds[1,3]) , "]", "\n\n")
  cat("Trial level surrogacy: ", "\n\n")
  cat("R Square: ", sprintf("%.4f", x$Trial.R2[1,1]),"[", sprintf("%.4f", x$Trial.R2[1,2]),";", sprintf("%.4f", x$Trial.R2[1,3]) , "]", "\n\n")

  cat("Estimated treatment effects on surrogate (resp_est) and survival (surv_est) endpoint: \n\n")
  print(x$EstTreatEffects[,c(1,2,3,4,5,7)], row.names = FALSE)
}

#' Generates a plot of the estimated treatment effects for the surrogate endpoint versus the estimated treatment effects for the true endpoint for an object fitted with the 'MetaAnalyticSurvBin()' function.
#'
#' @method plot MetaAnalyticSurvBin
#'
#' @param x An object of class 'MetaAnalyticSurvBin' fitted with the 'MetaAnalyticSurvBin()' function.
#' @param ... ...
#'
#' @return A plot of the type ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' data("colorectal")
#' fit_bin <- MetaAnalyticSurvBin(data = colorectal, true = surv, trueind = SURVIND, surrog = responder,
#'                                trt = TREAT, center = CENTER, trial = TRIAL, patientid = patientid,
#'                                adjustment="unadjusted")
#' plot(fit_bin)
#' }
#'
plot.MetaAnalyticSurvBin <- function(x,...){
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    resp_est <- surv_est <- sample_size <- NULL
    estimated_treatment_effects <- x$EstTreatEffects

    estimated_treatment_effects$sample_size <- as.numeric(estimated_treatment_effects$sample_size)
    estimated_treatment_effects$surv_est <- as.numeric(estimated_treatment_effects$surv_est)
    estimated_treatment_effects$resp_est <- as.numeric(estimated_treatment_effects$resp_est)

    # Create the scatter plot
    suppressWarnings(ggplot2::ggplot(data = estimated_treatment_effects, ggplot2::aes(x = resp_est, y = surv_est, size = sample_size)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", se = FALSE, color = "royalblue3") +
      ggplot2::labs(x = "Treatment effect on surrogate", y = "Treatment effect on true") +
      ggplot2::ggtitle("Treatment effect on true endpoint vs. treatment effect on surrogate endpoint") +
      ggplot2::theme(legend.position="none"))

  } else {
    stop("ggplot2 is not installed. Please install ggplot2 to use this function.")
  }
}


