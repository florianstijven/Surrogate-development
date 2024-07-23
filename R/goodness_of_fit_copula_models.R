#' Produce marginal GoF plot
#'
#' @param marginal Estimated marginal distribution represented by a list with
#'   three elements in the following order: the estimated cdf, pdf, and inverse
#'   cdf.
#' @param observed Observed values. These are used for the histogram.
#' @param name Name of the endpoint (used in the plot title).
#' @param type Type of endpoint: `"ordinal"` or `"continuous"`
#' @param treat Value for the treatment indicator.
#'
#' @seealso [plot.vine_copula_fit()]
#' @return NULL
marginal_gof_copula = function(marginal, observed, name, type, treat) {
  # The type of plot depends on the type of endpoint.
  if (type == "ordinal") {
    K = length(unique(observed))
    plot(
      1:K,
      marginal$pmf(1:K),
      xlab = "Category",
      ylab = "Probability Mass Function",
      main = paste0(name, ", Treat = ", treat),
      ylim = c(0, 1),
      xaxp = c(1, K, K - 1)
    )
    points(1:K, sapply(1:K, function(x) mean(observed == x)), col = "red")
    n = length(observed)
    arrows(
      1:K,
      sapply(1:K, function(x) {
        pi = mean(observed == x)
        logit = log(pi / (1 - pi))
        se_logit = sqrt(1 / (pi * (1 - pi))) / sqrt(n)
        logit_upper = log(pi / (1 - pi)) + 1.96 * se_logit
        return(exp(logit_upper) / (1 + exp(logit_upper)))
      }),
      1:K,
      sapply(1:K, function(x) {
        pi = mean(observed == x)
        logit = log(pi / (1 - pi))
        se_logit = sqrt(1 / (pi * (1 - pi))) / sqrt(n)
        logit_lower = log(pi / (1 - pi)) - 1.96 * se_logit
        return(exp(logit_lower) / (1 + exp(logit_lower)))
      }),
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
  if (type == "continuous") {
    grid = seq(from = min(observed), to = max(observed), length.out = 2e2)
    hist(observed, main = paste0(name, ", Treat = ", treat), freq = FALSE)
    lines(grid, marginal$pdf(grid), col = "red")
    # legend(
    #   x = "topright",
    #   lty = 1,
    #   col = c("red"),
    #   legend = c("Model-Based Density"),
    # )
  }
}

#' Produce Associational GoF plot
#'
#'
#' @param fitted_submodel List returned by [fit_copula_submodel_OrdCont()],
#' [fit_copula_submodel_ContCont()], or [fit_copula_submodel_OrdOrd()].
#' @param endpoint_types Character vector with 2 elements indicating the type of endpoints. Each element is either
#' `"ordinal"` or `"continuous"`.
#' @inheritParams marginal_gof_copula
#'
#' @return NULL
#' @seealso [plot.vine_copula_fit()]
association_gof_copula = function(fitted_submodel, treat, endpoint_types) {
  if (all(endpoint_types == c("ordinal", "continuous"))) {
    # Compute model-based regression function
    grid = seq(from = min(fitted_submodel$data$Y), to = max(fitted_submodel$data$Y), length.out = 2e2)
    cond_mean = conditional_mean_copula_OrdCont(fitted_submodel, grid)
    # Compute semi-parametric regression function.
    K = length(unique(fitted_submodel$data$X))
    y = fitted_submodel$data$X - 1
    x = fitted_submodel$data$Y
    if (K == 2) family = stats::binomial()
    else family = stats::quasi()
    fit_gam = mgcv::gam(y~s(x), family = family)
    predictions = mgcv::predict.gam(fit_gam, type = "response", newdata = data.frame(x = grid), se.fit = TRUE)

    # Make plot.

    plot(
      x = x,
      y = y + 1,
      col = "gray",
      ylim = c(0.5, K + 0.5),
      xlab = "S",
      ylab = "E(T | S)",
      main = paste0("Treat = ", treat)
    )
    lines(grid, cond_mean, col = "red")
    lines(
      x = grid,
      y = predictions$fit + 1
    )
    lines(
      x = grid,
      y = predictions$fit + 1.96 * predictions$se.fit + 1,
      lty = 2
    )
    lines(
      x = grid,
      y = predictions$fit - 1.96 * predictions$se.fit + 1,
      lty = 2
    )
  }
  else if (all(endpoint_types == c("continuous", "continuous"))) {
    # Compute model-based regression function
    grid = seq(from = min(fitted_submodel$data$X), to = max(fitted_submodel$data$X), length.out = 2e2)
    cond_mean = conditional_mean_copula_ContCont(fitted_submodel, grid)
    # Compute semi-parametric regression function.
    x = fitted_submodel$data$X
    y = fitted_submodel$data$Y
    fit_gam = mgcv::gam(y~s(x), family = stats::quasi())
    predictions = mgcv::predict.gam(fit_gam, type = "response", newdata = data.frame(x = grid), se.fit = TRUE)
    # Make plot.
    plot(
      x = x,
      y = y,
      col = "gray",
      xlab = "S",
      ylab = "E(T | S)",
      main = paste0("Treat = ", treat)
    )
    lines(grid, cond_mean, col = "red")
    lines(
      x = grid,
      y = predictions$fit
    )
    lines(
      x = grid,
      y = predictions$fit + 1.96 * predictions$se.fit,
      lty = 2
    )
    lines(
      x = grid,
      y = predictions$fit - 1.96 * predictions$se.fit,
      lty = 2
    )
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
    # Plote the conditional probabilities P(T = t | S). We have one plot for
    # each possible value of T. In each plot, the conditional probabilities are
    # plotted for all possible values of S.
    for (y_value in seq_along(1:K_Y)) {
      plot(
        1:K_X,
        cond_pmf(1:K_X, rep(y_value, K_X)),
        xlab = "S",
        ylab = paste0("P(T = ", y_value, "|S)"),
        main = paste0("Treat = ", treat),
        ylim = c(0, 1),
        xaxp = c(1, K_X, K_X - 1)
      )
      pi_cond = sapply(1:K_X, function(x) {
        mean(fitted_submodel$data$X == x &
               fitted_submodel$data$Y == y_value) / marginal_prob_x_emp[x]
      })
      points(1:K_X, pi_cond, col = "red")
      # number of observations for each conditional probality.
      n_X = marginal_prob_x_emp * n
      # Add whiskers representing 95% CIs based on the empirical proportions.
      arrows(
        1:K_X,
        sapply(1:K_X, function(x) {
          pi = pi_cond[x]
          logit = log(pi / (1 - pi))
          se_logit = sqrt(1 / (pi * (1 - pi))) / sqrt(n_X[x])
          logit_upper = log(pi / (1 - pi)) + 1.96 * se_logit
          return(exp(logit_upper) / (1 + exp(logit_upper)))
        }),
        1:K_X,
        sapply(1:K_X, function(x) {
          pi = pi_cond[x]
          logit = log(pi / (1 - pi))
          se_logit = sqrt(1 / (pi * (1 - pi))) / sqrt(n_X[x])
          logit_upper = log(pi / (1 - pi)) - 1.96 * se_logit
          return(exp(logit_upper) / (1 + exp(logit_upper)))
        }),
        length = 0.05,
        angle = 90,
        code = 3,
        col = "red"
      )
    }
  }
}
