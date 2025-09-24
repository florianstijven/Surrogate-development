#' Produce marginal GoF plot
#'
#' @details
#' # Return Plotting Data
#'
#' If `return_data` is `TRUE`, this function will return a data frame that can be
#' used to create customized plots. The following variables are present in the
#' returned data frame:
#' * `observed`: The empirical proportions (`type = "ordinal"`). `NA` for
#' `type = "continuous"`.
#' * `upper_ci`, `lower_ci`: Upper limit of the 95% confidence interval for the empirical
#' proportions. Defaults to `NA` if `type = "continuous"`.
#' * `value`: Value for the continuous or ordinal variable.
#' * `model_based`: Estimated model-based density (`type = "continuous"`) or
#' proportions (`type = "ordinal"`)
#'
#'
#' @param marginal Estimated marginal distribution represented by a list with
#'   three elements in the following order: the estimated cdf, pdf, and inverse
#'   cdf.
#' @param observed Observed values. These are used for the histogram.
#' @param name Name of the endpoint (used in the plot title).
#' @param type Type of endpoint: `"ordinal"` or `"continuous"`
#' @param treat Value for the treatment indicator.
#' @param return_data (boolean) Return the data used in the goodness-of-fit plot
#' (without the plot itself). This is useful when the user wants to customize the
#' plots, e.g., using `ggplot2`. See Details.
#' @param grid (numeric) vector of values for the endpoint at which the
#'   model-based density is computed.
#' @param ... Extra arguments passed onto [plot()] or [hist()] for an ordinal
#'  and continuous endpoint, respectively.
#'
#' @seealso [plot.vine_copula_fit()]
#' @return NULL
marginal_gof_copula = function(marginal,
                               observed,
                               name,
                               type,
                               treat,
                               return_data = FALSE,
                               grid = NULL,
                               ...) {
  # Number of observations.
  n = length(observed)

  # The type of plot depends on the type of endpoint.
  if (type == "ordinal") {
    K = length(unique(observed))

    gof_data = data.frame(value = 1:K)
    gof_data$observed = sapply(1:K, function(x)
      mean(observed == x))
    gof_data$model_based = marginal$pmf(1:K)
    gof_data$upper_ci = sapply(1:K, function(x) {
      pi = mean(observed == x)
      logit = log(pi / (1 - pi))
      se_logit = sqrt(1 / (pi * (1 - pi))) / sqrt(n)
      logit_upper = log(pi / (1 - pi)) + 1.96 * se_logit
      return(exp(logit_upper) / (1 + exp(logit_upper)))
    })
    gof_data$lower_ci = sapply(1:K, function(x) {
      pi = mean(observed == x)
      logit = log(pi / (1 - pi))
      se_logit = sqrt(1 / (pi * (1 - pi))) / sqrt(n)
      logit_upper = log(pi / (1 - pi)) - 1.96 * se_logit
      return(exp(logit_upper) / (1 + exp(logit_upper)))
    })

    if (return_data) {
      return(gof_data)
    }
    else {
      plot(
        gof_data$value,
        gof_data$model_based,
        xlab = "Category",
        ylab = "Probability Mass Function",
        main = paste0(name, ", Treat = ", treat),
        ylim = c(0, 1),
        xaxp = c(1, K, K - 1),
        ...
      )
      points(1:K, gof_data$observed, col = "red")
      arrows(
        1:K,
        gof_data$upper_ci,
        1:K,
        gof_data$lower_ci,
        length = 0.05,
        angle = 90,
        code = 3,
        col = "red"
      )
      # legend(
      #   x = "topright",
      #   lty = 1,
      #   col = c("red", "black"),
      #   legend = c("Model-Based", "Empirical")
      # )
    }

  }
  if (type == "continuous") {
    if (is.null(grid)) {
      grid = seq(
        from = min(observed),
        to = max(observed),
        length.out = 2e2
      )
    }
    gof_data = data.frame(grid = grid)
    gof_data$observed = rep(NA, length(grid))
    gof_data$model_based = marginal$pdf(grid)
    gof_data$upper_ci = rep(NA, length(grid))
    gof_data$lower_ci = rep(NA, length(grid))

    if (return_data) {
      return(gof_data)
    }
    else {
      hist(observed,
           main = paste0(name, ", Treat = ", treat),
           freq = FALSE,
           ...)
      lines(grid, gof_data$model_based, col = "red")
      # legend(
      #   x = "topright",
      #   lty = 1,
      #   col = c("red"),
      #   legend = c("Model-Based Density"),
      # )
    }

  }
}

#' Produce Associational GoF plot
#'
#' @details
#' # Semi-Parametric Regression estimates
#'
#' See the documentation of [plot.vine_copula_fit()] for the default
#' semi-parametric estimators.
#'
#' # Return Plotting Data
#'
#' If `return_data` is `TRUE`, this function will return a data frame that can
#' be used to create customized plots. The following variables are present in
#' the returned data frame:
#' * `observed`: The semi-parametric estimate of the regression function
#' \deqn{E(T | S)}.
#' * `upper_ci`, `lower_ci`: Upper and lower limit of the pointwise 95%
#' confidence interval for the semi-parametric estimate of the regression
#' function.
#' * `value`: Value for the surrogate endpoint at which the estimates for the
#' regression function are evaluated.
#' * `model_based`: Model-based estimate of the regression function.
#'
#' @param fitted_submodel List returned by [fit_copula_submodel_OrdCont()],
#'   [fit_copula_submodel_ContCont()], or [fit_copula_submodel_OrdOrd()].
#' @param grid (numeric) vector of values for the (surrogate) endpoint at which
#'   the regression function is evaluated.
#' @inheritParams new_vine_copula_fit
#' @inheritParams marginal_gof_copula
#' @param ... Extra argument passed onto [plot()].
#'
#' @return NULL
#' @seealso [plot.vine_copula_fit()]
association_gof_copula = function(fitted_submodel, treat, endpoint_types, return_data = FALSE, grid = NULL, ...) {
  if (all(endpoint_types == c("ordinal", "continuous"))) {
    # Compute model-based regression function
    if (is.null(grid)) {
      grid = seq(from = min(fitted_submodel$data$Y), to = max(fitted_submodel$data$Y), length.out = 2e2)
    }
    gof_data = data.frame(grid = grid)
    gof_data$model_based = conditional_mean_copula_OrdCont(fitted_submodel, grid)
    # Compute semi-parametric regression function.
    K = length(unique(fitted_submodel$data$X))
    y = fitted_submodel$data$X - 1
    x = fitted_submodel$data$Y
    if (K == 2) family = stats::binomial()
    else family = stats::quasi()
    fit_gam = mgcv::gam(y~s(x), family = family)
    predictions = mgcv::predict.gam(fit_gam, type = "response", newdata = data.frame(x = grid), se.fit = TRUE)

    gof_data$observed = predictions$fit + 1
    gof_data$upper_ci = gof_data$observed + 1.96 * predictions$se.fit
    gof_data$lower_ci = gof_data$observed - 1.96 * predictions$se.fit

    if (return_data) {
      return(gof_data)
    }
    else {
      # Make plot.
      plot(
        x = x,
        y = y + 1,
        col = "gray",
        ylim = c(0.5, K + 0.5),
        xlab = "S",
        ylab = "E(T | S)",
        main = paste0("Treat = ", treat),
        ...
      )
      lines(gof_data$grid, gof_data$model_based, col = "red")
      lines(
        x = gof_data$grid,
        y = gof_data$observed
      )
      lines(
        x = grid,
        y = gof_data$upper_ci,
        lty = 2
      )
      lines(
        x = grid,
        y = gof_data$lower_ci,
        lty = 2
      )
    }
  }
  else if (all(endpoint_types == c("continuous", "continuous"))) {
    # Compute model-based regression function
    if (is.null(grid)) {
      grid = seq(from = min(fitted_submodel$data$Y), to = max(fitted_submodel$data$Y), length.out = 2e2)
    }
    gof_data = data.frame(grid = grid)
    gof_data$model_based = conditional_mean_copula_ContCont(fitted_submodel, grid)

    # Compute semi-parametric regression function.
    x = fitted_submodel$data$X
    y = fitted_submodel$data$Y
    fit_gam = mgcv::gam(y~s(x), family = stats::quasi())
    predictions = mgcv::predict.gam(fit_gam, type = "response", newdata = data.frame(x = grid), se.fit = TRUE)

    gof_data$observed = predictions$fit
    gof_data$upper_ci = gof_data$observed + 1.96 * predictions$se.fit
    gof_data$lower_ci = gof_data$observed - 1.96 * predictions$se.fit

    if (return_data) {
      return(gof_data)
    }
    else {
      # Make plot.
      plot(
        x = x,
        y = y,
        col = "gray",
        xlab = "S",
        ylab = "E(T | S)",
        main = paste0("Treat = ", treat)
      )
      lines(gof_data$grid, gof_data$model_based, col = "red")
      lines(
        x = gof_data$grid,
        y = gof_data$observed
      )
      lines(
        x = grid,
        y = gof_data$upper_ci,
        lty = 2
      )
      lines(
        x = grid,
        y = gof_data$lower_ci,
        lty = 2
      )
    }
  }
  else if (all(endpoint_types == c("ordinal", "ordinal"))) {
    # Compute semi-parametric regression function.
    K_X = length(unique(fitted_submodel$data$X))
    K_Y = length(unique(fitted_submodel$data$Y))
    x = fitted_submodel$data$X
    y = fitted_submodel$data$Y
    n = length(x)
    # Construct the estimated joint pmf.
    joint_pmf = function(x, y) {
      exp(
        ordinal_ordinal_loglik(
          para = coef(fitted_submodel$ml_fit),
          X = x,
          Y = y,
          copula_family = fitted_submodel$copula_family,
          K_X = K_X,
          K_Y = K_Y,
          return_sum = FALSE
        )
      )
    }
    marginal_prob_x = sapply(1:K_X, function(x) {
      sum(joint_pmf(rep(x, K_Y), 1:K_Y))
    })
    marginal_prob_x_emp = sapply(1:K_X, function(x) mean(fitted_submodel$data$X == x))
    cond_pmf = function(x, y) {
      joint_pmf(x, y) / marginal_prob_x[x]
    }

    # Create data frame to which new goodness-of-fit data will be added in a for
    # loop.
    gof_data = data.frame()

    # Plot the conditional probabilities P(T = t | S). We have one plot for
    # each possible value of T. In each plot, the conditional probabilities are
    # plotted for all possible values of S.
    for (y_value in seq_along(1:K_Y)) {
      gof_temp = data.frame(grid = 1:K_X)
      gof_temp$model_based = cond_pmf(1:K_X, rep(y_value, K_X))
      gof_temp$observed = sapply(1:K_X, function(x) {
        mean(fitted_submodel$data$X == x &
               fitted_submodel$data$Y == y_value) / marginal_prob_x_emp[x]
      })
      # number of observations for each conditional probability.
      n_X = marginal_prob_x_emp * n

      gof_temp$upper_ci = sapply(1:K_X, function(x) {
        pi = gof_temp$observed[x]
        logit = log(pi / (1 - pi))
        se_logit = sqrt(1 / (pi * (1 - pi))) / sqrt(n_X[x])
        logit_upper = log(pi / (1 - pi)) + 1.96 * se_logit
        return(exp(logit_upper) / (1 + exp(logit_upper)))
      })
      gof_temp$lower_ci = sapply(1:K_X, function(x) {
        pi = gof_temp$observed[x]
        logit = log(pi / (1 - pi))
        se_logit = sqrt(1 / (pi * (1 - pi))) / sqrt(n_X[x])
        logit_upper = log(pi / (1 - pi)) - 1.96 * se_logit
        return(exp(logit_upper) / (1 + exp(logit_upper)))
      })

      if (return_data) {
        # Add the surrogate-value specific goodness-of-fit data to the global
        # data frame that will eventually be returned.
        gof_data = rbind(gof_data, gof_temp)
      }
      else {
        plot(
          gof_temp$grid,
          gof_temp$model_based,
          xlab = "S",
          ylab = paste0("P(T = ", y_value, "|S)"),
          main = paste0("Treat = ", treat),
          ylim = c(0, 1),
          xaxp = c(1, K_X, K_X - 1)
        )
        points(gof_temp$grid, gof_temp$observed, col = "red")

        # Add whiskers representing 95% CIs based on the empirical proportions.
        arrows(
          gof_temp$grid,
          gof_temp$upper_ci,
          gof_temp$grid,
          gof_temp$lower_ci,,
          length = 0.05,
          angle = 90,
          code = 3,
          col = "red"
        )
      }
    }
    # Return the GoF data if asked.
    if (return_data) {
      return(gof_data)
    }
  }
}
