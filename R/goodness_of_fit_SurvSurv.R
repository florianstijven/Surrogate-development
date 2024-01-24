#' Marginal survival function goodness of fit
#'
#' The [marginal_gof_plots_scr()] function plots the estimated marginal survival
#' functions for the fitted model. This results in four plots of survival
#' functions, one for each of \eqn{S_0}, \eqn{S_1}, \eqn{T_0}, \eqn{T_1}.
#'
#' @param fitted_model Returned value from [fit_model_SurvSurv()]. This object
#'   essentially contains the estimated identifiable part of the joint
#'   distribution for the potential outcomes.
#' @param grid grid of time-points for which to compute the estimated survival
#'   functions.
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
#' grid = seq(from = 0, to = 2, length.out = 50)
#' marginal_gof_plots_scr(ovarian_fitted, grid)
#'
#' @importFrom survival survfit Surv
marginal_gof_plots_scr = function(fitted_model,
                                  grid) {
  marginal_gof_scr_T_plot(fitted_model = fitted_model,
                          grid = grid,
                          treated = 0,
                          xlim = c(0, max(grid)),
                          main = "Survival Function T - Control Treatment",
                          xlab = "t",
                          ylab = "P(T > t)")
  marginal_gof_scr_T_plot(fitted_model = fitted_model,
                          grid = grid,
                          treated = 1,
                          xlim = c(0, max(grid)),
                          main = "Survival Function T - Active Treatment",
                          xlab = "t",
                          ylab = "P(T > t)")
  marginal_gof_scr_S_plot(fitted_model = fitted_model,
                          grid = grid,
                          treated = 0,
                          main = "Survival Function S - Control Treatment",
                          xlab = "s",
                          ylab = "P(S > s)")
  marginal_gof_scr_S_plot(fitted_model = fitted_model,
                          grid = grid,
                          treated = 1,
                          main = "Survival Function S - Active Treatment",
                          xlab = "s",
                          ylab = "P(S > s)")
}


marginal_gof_scr_S_plot = function(fitted_model,
                                   grid,
                                   treated,
                                   ...) {
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

  # Helper function that computes P(S_tilde > x, T > x).
  para_t = para[(1 + k):(2 * k)]
  surv_joint = function(x) {
    exp(
      survival_survival_loglik(
        para = para,
        X = x,
        delta_X = 0,
        Y = x,
        delta_Y = 0,
        copula_family = fitted_model$copula_family[treated + 1],
        knotsx = knots,
        knotsy = knott
      )
    )
  }
  #goodness of fit of marginal survival function of OS
  Surv_probs = sapply(
    X = grid,
    FUN = surv_joint
  )
  plot(survival::survfit(
    survival::Surv(fitted_model$data$Pfs, pmax(fitted_model$data$PfsInd, fitted_model$data$SurvInd)) ~ 1,
    data = fitted_model$data,
    subset = fitted_model$data$Treat == treated
  ),
  ...)
  lines(grid, Surv_probs, col = "red")
  legend(
    x = "topright",
    lty = 1,
    col = c("red", "black"),
    legend = c("Model-Based", "Kaplan-Meier")
  )
}

marginal_gof_scr_T_plot = function(fitted_model,
                                   grid,
                                   treated,
                                   ...) {
  if (treated == 0) {
    fitted_submodel = fitted_model$fit_0
    knott = fitted_model$knott0
  }  else {
    fitted_submodel = fitted_model$fit_1
    knott = fitted_model$knott1
  }
  # Extract estimated model parameters
  para = fitted_submodel$estimate
  k = length(knott)

  # Helper function that computes the survival function for T.
  para_t = para[(1 + k):(2 * k)]
  surv_t = function(x) {
    flexsurv::psurvspline(q = x, gamma = para_t, knots = knott)
  }
  #goodness of fit of marginal survival function of OS
  Surv_probs = 1 - surv_t(grid)
  plot(
    survival::survfit(
      survival::Surv(Surv, SurvInd) ~ 1,
      data = fitted_model$data,
      subset = fitted_model$data$Treat == treated
    ),
    ...
  )
  lines(grid, Surv_probs, col = "red")
  legend(
    x = "topright",
    lty = 1,
    col = c("red", "black"),
    legend = c("Model-Based", "Kaplan-Meier")
  )
}


mean_S_before_T_plots_scr = function(fitted_model,
                                     plot_method = scatter.smooth,
                                     grid,
                                     ...)
{
  mean_S_before_T_plot_scr(
    fitted_model = fitted_model,
    grid = grid,
    treated = 0,
    xlab = "t",
    ylab = "E(S | T = t, S < T)",
    col = "gray",
    main = "Control Treatment"
  )
  mean_S_before_T_plot_scr(
    fitted_model = fitted_model,
    grid = grid,
    treated = 1,
    xlab = "t",
    ylab = "E(S | T = t, S < T)",
    col = "gray",
    main = "Active Treatment"
  )

}

mean_S_before_T_plot_scr = function(fitted_model,
                                    plot_method = scatter.smooth,
                                    grid,
                                    treated,
                                    ...){
  # Select appropriate subpopulation.
  selected_data = fitted_model$data[fitted_model$data$PfsInd == 1 & fitted_model$data$SurvInd == 1, ]

  # Compute the expected value at the grid points.
  model_cond_means = sapply(
    X = grid,
    FUN = mean_S_before_T,
    fitted_model = fitted_model,
    treated = treated
  )

  # Plot the smooth and model-based estimates
  plot_method(
    x = selected_data$Surv[selected_data$Treat == treated],
    y = selected_data$Pfs[selected_data$Treat == treated],
    ...
  )
  lines(x = grid, y = model_cond_means, col = "red")
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
      copula_family = fitted_model$copula_family[treated + 1],
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
        copula_family = fitted_model$copula_family[treated + 1],
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




