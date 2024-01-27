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
#' Surrogate:::marginal_gof_plots_scr(ovarian_fitted, grid)
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


#' Goodness-of-fit plot for the marginal survival functions
#'
#' The [marginal_gof_scr_S_plot()] and [marginal_gof_scr_T_plot()] functions
#' plot the estimated marginal survival functions for the surrogate and true
#' endpoints. In these plots, it is assumed that the copula model has been
#' fitted for \eqn{(T_0, \tilde{S}_0, \tilde{S}_1, T_1)'} where \deqn{S_k =
#' \min(\tilde{S_k}, T_k)} is the (composite) surrogate of interest. In these
#' plots, the model-based survival functions for \eqn{(T_0, S_0, S_1, T_1)'} are
#' plotted together with the corresponding Kaplan-Meier etimates.
#'
#' @param fitted_model Returned value from [fit_model_SurvSurv()]. This object
#'   essentially contains the estimated identifiable part of the joint
#'   distribution for the potential outcomes.
#' @param grid Grid of time-points at which the model-based estimated regression
#'   functions, survival functions, or probabilities are evaluated.
#' @param treated (numeric) Treatment group. Should be `0` or `1`.
#' @param ... Additional arguments to pass to [plot()].
#'
#' @details # True Endpoint
#'
#'   The marginal goodness-of-fit plots for the true endpoint, build by
#'   [marginal_gof_scr_T_plot()], is simply a comparison of the model-based
#'   estimate of \eqn{P(T_k > t)} with the Kaplan-Meier (KM) estimate obtained
#'   with [survival::survfit()]. A pointwise 95% confidence interval for the KM
#'   estimate is also plotted.
#'
#'   # Surrogate Endpoint
#'
#'   The model-based estimate of \eqn{P(S_k > s)} follows indirectly from the
#'   fitted copula model because the copula model has been fitted for
#'   \eqn{\tilde{S}_k} instead of \eqn{S_k}. However, the model-based estimate
#'   still follows easily from the copula model as follows,
#'   \deqn{P(S_k > s) = P(\min(\tilde{S}_k, T_k)) = P(\tilde{S}_k > s, T_k > s).}
#'
#'   The [marginal_gof_scr_T_plot()] function plots the model-based estimate for
#'   \eqn{P(\tilde{S}_k > s, T_k > s)} together with the KM estimate (see above).
#'
#'
#' @return `NULL`
#' @export
#' @inherit mean_S_before_T_plot_scr examples
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

#' @rdname marginal_gof_scr_S_plot
#' @export
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
                                     plot_method = NULL,
                                     grid,
                                     ...)
{
  mean_S_before_T_plot_scr(
    fitted_model = fitted_model,
    plot_method = plot_method,
    grid = grid,
    treated = 0,
    xlab = "t",
    ylab = "E(S | T = t, S < T)",
    col = "gray",
    main = "Control Treatment"
  )
  mean_S_before_T_plot_scr(
    fitted_model = fitted_model,
    plot_method = plot_method,
    grid = grid,
    treated = 1,
    xlab = "t",
    ylab = "E(S | T = t, S < T)",
    col = "gray",
    main = "Active Treatment"
  )

}

#' Goodness of fit plot for the fitted copula
#'
#' The [mean_S_before_T_plot_scr()] and [prob_dying_without_progression_plot()]
#' functions build plots to assess the goodness-of-fit of the copula model
#' fitted by [fit_model_SurvSurv()]. Specifically, these two functions focus on
#' the appropriateness of the copula. Note that to assess the appropriateness of
#' the marginal functions, two other functions are available:
#' [marginal_gof_scr_S_plot()] and [marginal_gof_scr_T_plot()].
#'
#' @param plot_method Defaults to `NULL`. Should not be modified.
#' @inheritParams marginal_gof_scr_S_plot
#'
#' @details # Progression Before Death
#'
#' If a patient progresses before death, this means that \eqn{S_k < T_k}. For
#' these patients, we can look at the expected progression time given that the
#' patient has died at \eqn{T_k = t}: \deqn{E(S_k | T_k = t, S_k < T_k).} The
#' [mean_S_before_T_plot_scr()] function plots the model-based estimate of this
#' regression function together with a non-parametric estimate.
#'
#' This regression function can also be estimated non-parametrically by
#' regressing \eqn{S_k} onto \eqn{T_k} in the subset of uncensored patients.
#' This non-parametric estimate is obtained via `mgcv::gam(y~s(x))` with
#' additionally `family = stats::quasi(link = "log", variance = "mu")` because
#' this tends to describe survival data better. The 95% confidence intervals are
#' added for this non-parametric estimate; although, they should be interpreted
#' with caution because the Poisson mean-variance relation may be wrong.
#'
#' # Death Before Progression
#'
#' If a patient dies before progressing, this means that \eqn{S_k = T_k}. This
#' probability can be modeled as a function of time, i.e., \deqn{	\pi_k(t) & =
#' P(S_k = t  \, | \, T_k = t).} The [prob_dying_without_progression_plot()]
#' function plots the model-based estimate of this regression function together
#' with a non-parametric estimate.
#'
#' This regression function can also be estimated non-parametrically by
#' regressing the censoring indicator for \eqn{S_k}, \eqn{\delta_{S_k}},
#' onto \eqn{T_k} in the subset of patients with uncensored \eqn{T_k}.
#'
#' @return `NULL`
#' @export
#'
#' @examples
#' # Load Ovarian data
#' data("Ovarian")
#' # Recode the Ovarian data in the semi-competing risks format.
#' data_scr = data.frame(
#'   ttp = Ovarian$Pfs,
#'   os = Ovarian$Surv,
#'   treat = Ovarian$Treat,
#'   ttp_ind = ifelse(
#'     Ovarian$Pfs == Ovarian$Surv &
#'       Ovarian$SurvInd == 1,
#'     0,
#'     Ovarian$PfsInd
#'   ),
#'   os_ind = Ovarian$SurvInd
#' )
#' # Fit copula model.
#' fitted_model = fit_model_SurvSurv(data = data_scr,
#'                                   copula_family = "clayton",
#'                                   n_knots = 1)
#' # Define grid for GoF plots.
#' grid = seq(from = 1e-3,
#'            to = 2.5,
#'            length.out = 30)
#' # Assess marginal goodness-of-fit in the control group.
#' marginal_gof_scr_S_plot(fitted_model, grid = grid, treated = 0)
#' marginal_gof_scr_T_plot(fitted_model, grid = grid, treated = 0)
#' # Assess goodness-of-fit of the association structure, i.e., the copula.
#' prob_dying_without_progression_plot(fitted_model, grid = grid, treated = 0)
#' mean_S_before_T_plot_scr(fitted_model, grid = grid, treated = 0)
mean_S_before_T_plot_scr = function(fitted_model,
                                    plot_method = NULL,
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
  if (is.null(plot_method)) {
    if(requireNamespace("mgcv", quietly = TRUE)){
      plot_method = function(x, y, ...) {
        # Fit GAM
        fit_gam = mgcv::gam(y~s(x), family = stats::quasi(link = "log", variance = "mu"))
        # Plot smooth estimated.
        predictions = mgcv::predict.gam(fit_gam, newdata = data.frame(x = grid), type = "link", se.fit = TRUE)
        plot(
          x = grid,
          y = exp(predictions$fit),
          type = "l",
          ylim = c(0, max(selected_data$Surv)),
          ...
        )
        lines(
          x = grid,
          y = exp(predictions$fit + 1.96 * predictions$se.fit),
          lty = 2
        )
        lines(
          x = grid,
          y = exp(predictions$fit - 1.96 * predictions$se.fit),
          lty = 2
        )
        points(x = x, y = y, col = "gray")
      }
    } else {
      warning("The R package mgcv is not installed. Some plots may not be displayed.")
    }
  }
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
  # Compute P(S_Tilde < T | t = t). This is the probability of progressing
  # before dying, given that the patient died at t.
  prob_progression_t = prob_progression_before_dying(fitted_model = fitted_model,
                                                     t = t,
                                                     treated = treated)

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

prob_progression_before_dying = function(t, fitted_model, treated) {
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
  return(prob_progression_t)
}

#' @rdname mean_S_before_T_plot_scr
#' @export
prob_dying_without_progression_plot = function(fitted_model,
                                              plot_method = NULL,
                                              grid,
                                              treated,
                                              ...) {
  # Select appropriate subpopulation.
  selected_data = fitted_model$data[fitted_model$data$SurvInd == 1, ]

  # Compute the expected value at the grid points.
  model_prob = 1 - sapply(
    X = grid,
    FUN = prob_progression_before_dying,
    fitted_model = fitted_model,
    treated = treated
  )

  # If no plot_method has been specified, use a generalized additive linear
  # model with penalized splines.
  if (is.null(plot_method)) {
    if(requireNamespace("mgcv", quietly = TRUE)){
      plot_method = function(x, y, ...) {
        # Fit GAM
        fit_gam = mgcv::gam(y~s(x), family = stats::binomial())
        expit = function(x) 1 / (1 + exp(-x))
        # Plot smooth estimated.
        predictions = mgcv::predict.gam(fit_gam, newdata = data.frame(x = grid), type = "link", se.fit = TRUE)
        plot(
          x = grid,
          y = expit(predictions$fit),
          type = "l",
          ylim = c(0, 1),
          ...
        )
        lines(
          x = grid,
          y = expit(predictions$fit + 1.96 * predictions$se.fit),
          lty = 2
        )
        lines(
          x = grid,
          y = expit(predictions$fit - 1.96 * predictions$se.fit),
          lty = 2
        )
        points(x = x, y = y, col = "gray")
      }
    } else {
      warning("The R package mgcv is not installed. Some plots may not be displayed.")
    }
  }
  # Plot the smooth and model-based estimates.
  plot_method(
    x = selected_data$Surv[selected_data$Treat == treated],
    y = 1 - selected_data$PfsInd[selected_data$Treat == treated],
    ...
  )
  lines(x = grid, y = model_prob, col = "red")
  legend(
    x = "topright",
    lty = 1,
    bty = "n",
    col = c("red", "black"),
    legend = c("Model-Based", "Smooth Estimated")
  )
}

prob_dying_without_progression_plots = function(fitted_model,
                                                plot_method = NULL,
                                                grid,
                                                treated,
                                                ...) {
  prob_dying_without_progression_plot(
    fitted_model = fitted_model,
    plot_method = plot_method,
    grid = grid,
    treated = 0,
    xlab = "t",
    ylab = "P(S = T | T = t)",
    col = "gray",
    main = "Control Treatment"
  )
  prob_dying_without_progression_plot(
    fitted_model = fitted_model,
    plot_method = plot_method,
    grid = grid,
    treated = 1,
    xlab = "t",
    ylab = "P(S = T | T = t)",
    col = "gray",
    main = "Active Treatment"
  )
}




