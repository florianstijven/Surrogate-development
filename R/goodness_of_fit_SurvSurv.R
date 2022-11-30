#' Marginal survival function goodness of fit
#'
#' The [marginal_gof_scr()] function plots the estimated marginal survival
#' functions for the fitted model. This results in four plots of survival
#' functions, one for each of \eqn{S_0}, \eqn{S_1}, \eqn{T_0}, \eqn{T_1}.
#'
#' @param fitted_model Returned value from [fit_model_SurvSurv()]. This object
#'   essentially contains the estimated identifiable part of the joint
#'   distribution for the potential outcomes.
#' @param data data that was supplied to [fit_model_SurvSurv()].
#' @param grid grid of time-points for which to compute the estimated survival
#'   functions.
#' @param time_unit character vector that reflects the time unit of the
#'   endpoints, defaults to `"years"`.
#'
#' @export
#'
#' @examples
#' library(Surrogate)
#' data("Ovarian")
#' #For simplicity, data is not recoded to semi-competing risks format, but is
#' #left in the composite event format.
#' data = data.frame(
#'   Ovarian$Pfs,
#'   Ovarian$Surv,
#'   Ovarian$Treat,
#'   Ovarian$PfsInd,
#'   Ovarian$SurvInd
#' )
#' ovarian_fitted =
#'   fit_model_SurvSurv(data = data,
#'                      copula_family = "clayton",
#'                      nknots = 1)
#' grid = seq(from = 0, to = 2, length.out = 200)
#' marginal_gof_scr(ovarian_fitted, data, grid)
#'
#' @importFrom graphics lines par
marginal_gof_scr = function(fitted_model,
                            data, grid, time_unit = "years"){
  nknots = length(fitted_model$knots0)
  par(mfrow = c(2, 2))
  colnames(data) = c("Pfs", "Surv", "Treat", "PfsInd", "SurvInd")
  #time unit label for x-axis
  label_time_unit = paste0("Time (", time_unit, ")")

  #goodness of fit of marginal survival function of OS
  Surv_t0 = 1 - flexsurv::psurvspline(q = grid,
                                      gamma = fitted_model$parameters0[(nknots + 1):(2 *
                                                                                       nknots)],
                                      knots = fitted_model$knott0)
  plot(
    survival::survfit(
      survival::Surv(Surv, SurvInd) ~ 1,
      data = data,
      subset = data$Treat == 0
    ),
    xlim = c(min(grid), max(grid)),
    main = "T(0)",
    xlab = label_time_unit,
    ylab = "S(t)"
  )
  lines(grid, Surv_t0, col = "red")

  Surv_t1 = 1 - flexsurv::psurvspline(q = grid,
                                      gamma = fitted_model$parameters1[(nknots + 1):(2 *
                                                                                       nknots)],
                                      knots = fitted_model$knott1)
  plot(
    survival::survfit(
      survival::Surv(Surv, SurvInd) ~ 1,
      data = data,
      subset = data$Treat == 1
    ),
    xlim = c(min(grid), max(grid)),
    main = "T(1)",
    xlab = label_time_unit,
    ylab = "S(t)"
  )
  lines(grid, Surv_t1, col = "red")

  #goodness of fit for marginal distribution of PFS
  pfs_surv = function(s,
                      gammas,
                      gammat,
                      knots,
                      knott,
                      theta,
                      copula_family) {
    u = 1 - flexsurv::psurvspline(q = s, gamma = gammas, knots = knots)
    v = 1 - flexsurv::psurvspline(q = s, gamma = gammat, knots = knott)
    uv = matrix(data = c(u, v), ncol = 2, byrow = FALSE)
    if (copula_family == "frank") {
      C = copula::pCopula(copula = copula::frankCopula(theta),
                          u = uv)
      # C = (-1/theta)*log(((1-exp(-theta)-(1-exp(-theta*u))*(1-exp(-theta*v))))/(1-exp(-theta)))
    }
    else if (copula_family == "gaussian") {
      rho = (exp(theta) - 1) / (exp(theta) + 1)
      C = copula::pCopula(copula = copula::normalCopula(rho),
                          u = uv)
      # Sigma = matrix(data = c(1, rho, rho, 1), ncol = 2)
      # V = qnorm(c(u, v))
      # C = pmvnorm(lower = -Inf,upper=V, sigma=Sigma, mean=c(0,0))[1]
    }
    else if (copula_family == "clayton") {
      C = copula::pCopula(copula = copula::claytonCopula(theta),
                          u = uv)
      # C = (u^(-theta) + v^(-theta) - 1)^(-1/theta)
    }
    else if (copula_family == "gumbel") {
      C = copula::pCopula(copula = copula::gumbelCopula(theta),
                          u = uv)
      # C = exp(-((-log(u))^(theta)+(-log(v))^(theta))^(1/theta))
    }
    return(C)
  }

  probs0 = pfs_surv(s = grid,
                    gammas = fitted_model$parameters0[1:nknots],
                    gammat = fitted_model$parameters0[(nknots + 1):(2 * nknots)],
                    knots = fitted_model$knots0,
                    knott = fitted_model$knott0,
                    theta = fitted_model$parameters0[2 * nknots + 1],
                    copula_family = fitted_model$copula_family)

  # probs0 = sapply(
  #   grid,
  #   pfs_surv,
  #   gammas = fitted_model$parameters0[1:nknots],
  #   gammat = fitted_model$parameters0[(nknots + 1):(2 * nknots)],
  #   knots = fitted_model$knots0,
  #   knott = fitted_model$knott0,
  #   theta = fitted_model$parameters0[2 * nknots + 1],
  #   copula_family = fitted_model$copula_family
  # )
  plot(
    survival::survfit(
      survival::Surv(data$Pfs, pmax(data$PfsInd, data$SurvInd)) ~ 1,
      data = data,
      subset = data$Treat == 0
    ),
    xlim = c(min(grid), max(grid)),
    main = "PFS (0)",
    xlab = label_time_unit,
    ylab = "S(t)"
  )
  lines(grid, probs0, col = "red")

  probs1 = pfs_surv(s = grid,
                    gammas = fitted_model$parameters1[1:nknots],
                    gammat = fitted_model$parameters1[(nknots + 1):(2 * nknots)],
                    knots = fitted_model$knots1,
                    knott = fitted_model$knott1,
                    theta = fitted_model$parameters1[2 * nknots + 1],
                    copula_family = fitted_model$copula_family)

  # probs1 = sapply(
  #   grid,
  #   pfs_surv,
  #   gammas = fitted_model$parameters1[1:nknots],
  #   gammat = fitted_model$parameters1[(nknots + 1):(2 * nknots)],
  #   knots = fitted_model$knots1,
  #   knott = fitted_model$knott1,
  #   theta = fitted_model$parameters1[2 * nknots + 1],
  #   copula_family = fitted_model$copula_family
  # )
  plot(
    survival::survfit(
      survival::Surv(data$Pfs, pmax(data$PfsInd, data$SurvInd)) ~ 1,
      data = data,
      subset = data$Treat == 1
    ),
    xlim = c(min(grid), max(grid)),
    main = "PFS (1)",
    xlab = label_time_unit,
    ylab = "S(t)"
  )
  lines(grid, probs1, col = "red")
}


