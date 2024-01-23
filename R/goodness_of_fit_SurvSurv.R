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
#'                      n_knots = 1)
#' grid = seq(from = 0, to = 2, length.out = 200)
#' marginal_gof_scr(ovarian_fitted, data, grid)
#'
#' @importFrom graphics lines
#' @importFrom flexsurv psurvspline
#' @importFrom survival survfit Surv
#' @importFrom copula pcopula
marginal_gof_scr = function(fitted_model,
                            data, grid, time_unit = "years"){
  n_knots = length(fitted_model$knots0)
  colnames(data) = c("Pfs", "Surv", "Treat", "PfsInd", "SurvInd")
  #time unit label for x-axis
  label_time_unit = paste0("Time (", time_unit, ")")

  #goodness of fit of marginal survival function of OS
  Surv_t0 = 1 - flexsurv::psurvspline(
    q = grid,
    gamma = coef(fitted_model$fit_0)[(n_knots + 1):(2 * n_knots)],
    knots = fitted_model$knott0
  )
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

  Surv_t1 = 1 - flexsurv::psurvspline(
    q = grid,
    gamma = coef(fitted_model$fit_1)[(n_knots + 1):(2 * n_knots)],
    knots = fitted_model$knott1
  )
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
                    gammas = coef(fitted_model$fit_0)[1:n_knots],
                    gammat = coef(fitted_model$fit_0)[(n_knots + 1):(2 * n_knots)],
                    knots = fitted_model$knots0,
                    knott = fitted_model$knott0,
                    theta = coef(fitted_model$fit_0)[2 * n_knots + 1],
                    copula_family = fitted_model$copula_family)

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
                    gammas = coef(fitted_model$fit_1)[1:n_knots],
                    gammat = coef(fitted_model$fit_1)[(n_knots + 1):(2 * n_knots)],
                    knots = fitted_model$knots1,
                    knott = fitted_model$knott1,
                    theta = coef(fitted_model$fit_1)[2 * n_knots + 1],
                    copula_family = fitted_model$copula_family)

  # probs1 = sapply(
  #   grid,
  #   pfs_surv,
  #   gammas = fitted_model$parameters1[1:n_knots],
  #   gammat = fitted_model$parameters1[(n_knots + 1):(2 * n_knots)],
  #   knots = fitted_model$knots1,
  #   knott = fitted_model$knott1,
  #   theta = fitted_model$parameters1[2 * n_knots + 1],
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

mean_S_before_T_plots_scr = function(fitted_model,
                                     plot_method = scatter.smooth,
                                     grid,
                                     ...)
{
  # Select appropriate subpopulation.
  selected_data = fitted_model$data[fitted_model$data$PfsInd == 1 & fitted_model$data$SurvInd == 1, ]

  # Compute the expected value at the grid points.
  model_cond_means0 = sapply(
    X = grid,
    FUN = mean_S_before_T,
    fitted_model = fitted_model,
    treated = 0
  )
  model_cond_means1 = sapply(
    X = grid,
    FUN = mean_S_before_T,
    fitted_model = fitted_model,
    treated = 1
  )

  # Plot the smooth and model-based estimates
  plot_method(
    x = selected_data$Surv[selected_data$Treat == 0],
    y = selected_data$Pfs[selected_data$Treat == 0],
    xlab = "t",
    ylab = "E(S | T = t, S < T)",
    col = "gray",
    main = "Control Treatment"
  )
  lines(x = grid, y = model_cond_means0, col = "red")
  legend(
    x = "topleft",
    lty = 1,
    col = c("red", "black"),
    legend = c("Model-Based", "Smooth Estimated")
  )

  plot_method(
    x = selected_data$Surv[selected_data$Treat == 1],
    y = selected_data$Pfs[selected_data$Treat == 1],
    xlab = "t",
    ylab = "E(S | T = t, S < T)",
    col = "gray",
    main = "Experimental Treatment", degree = 2
  )
  lines(x = grid, y = model_cond_means1, col = "red")
  legend(
    x = "topleft",
    lty = 1,
    col = c("red", "black"),
    legend = c("Model-Based", "Smooth Estimated")
  )

}

mean_S_before_T = function(t, fitted_model, treated) {
  if (treated == 0) {
    fitted_submodel = fitted_model$fit_0
    knots = fitted_model$knots0
    knott = fitted_model$knott0
    }  else {
    fitted_submodel = fitted_model$fit_1
    knots = fitted_model$knots1
    knott = fitted_model$knott1
  }
  # Extract estimated model parameters
  para = fitted_submodel$estimate
  k = length(knott)

  # Helper function that computes the density for T = t.
  para_t = para[(1 + k):(2 * k)]
  dens_t = function(x) {
    flexsurv::dsurvspline(x = x, gamma = para_t, knots = knott)
  }
  # Compute probability of S_Tilde < T. This is the probability of progressing
  # before dying, given that the patien died at t. We re-use the log likelihood
  # functions to compute this probability.

  # Probability of (S_tilde > t, T = t).
  prob1 = exp(
    survival_survival_loglik(
      para = para,
      X = t,
      delta_X = 0,
      Y = t,
      delta_Y = 1,
      copula_family = fitted_model$copula_family,
      knotsx = knots,
      knotsy = knott
    )
  )
  prob_progression_t =  1 - (prob1 / dens_t(t))

  # Compute the expected value through numerical integration. First, define a
  # helper function that compute the probability of (S_tilde = s, T = t).
  # Second, define the integrand.
  dens_joint = function(s, t) {
    exp(
      survival_survival_loglik(
        para = para,
        X = s,
        delta_X = 1,
        Y = t,
        delta_Y = 1,
        copula_family = fitted_model$copula_family,
        knotsx = knots,
        knotsy = knott
      )
    )
  }
  integrand = function(x) {
    (1 / prob_progression_t) *
      x * sapply(x, dens_joint, t = t) / dens_t(t)
  }

  mean_S_given_T = stats::integrate(
    f = integrand,
    lower = 0,
    upper = t
  )$value

  return(mean_S_given_T)
}




