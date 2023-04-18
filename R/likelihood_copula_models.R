#' Computes loglikelihood for a given copula model
#'
#' `log_likelihood_copula_model()` computes the loglikelihood for a given
#' bivariate copula model and data set while allowin for right-censoring of both
#' outcome variables.
#'
#' @param X Numeric vector corresponding to first outcome variable.
#' @param Y Numeric vector corresponding to second outcome variable.
#' @param cdf_X Distribution function for the first outcome variable.
#' @param cdf_Y Distribution function for the second outcome variable.
#' @param pdf_X Density function for the first outcome variable.
#' @param pdf_Y Density function for the second outcome variable.
#' @inheritParams loglik_copula_scale
#'
#' @return Loglikelihood of the bivariate copula model evaluated in the observed
#'   data.
log_likelihood_copula_model = function(theta,
                                       X,
                                       Y,
                                       d1,
                                       d2,
                                       copula,
                                       cdf_X,
                                       cdf_Y,
                                       pdf_X,
                                       pdf_Y) {
  # The loglikelihood contribution in a copula model can be separated
  # into a part that depends only on the copula parameters, and a part
  # that only depends on the marginal distribution parameters. We can compute
  # these two parts separately and then sum them.

  # Transform observations on original scale to the copula scale.
  u = cdf_X(X)
  v = cdf_Y(Y)

  # loglikelihood contribution for the copula part.
  loglik_copula = loglik_copula_scale(theta, u, v, d1, d2, copula)


  # Log likelihood contribution for the marginal distribution part.
  # uncensored observations.
  part1 <-
    ifelse(d1 * d2 == 1,
           log(pdf_X(X)) + log(pdf_Y(Y)),
           0)
  # second observation censored.
  part2 <-
    ifelse(d1 * (1 - d2) == 1,
           log(pdf_X(X)),
           0)
  # first observation censored.
  part3 <-
    ifelse((1 - d1) * d2 == 1,
           log(pdf_Y(Y)),
           0)
  # When both observations are censored the likelihood does not depend
  # (directly) on the marginal part. The dependence on the marginal part is only
  # through the copula itself.
  loglik = sum(part1 + part2 + part3) + loglik_copula
}

#' Function factory for distribution functions
#'
#' @param mean Mean parameter.
#' @param extra_par Additional parameters. The specific meaning of these extra
#'   parameters depends on `family`.
#' @param family Distributional family, one of the following:
#' * `"normal"`: normal distribution where `mean` is the mean and `extra_par` is
#'   the standard deviation.
#' * `"logistic"`: logistic distribution as parameterized in `stats::plogis()` where
#'   `mean` and `extra_par` correspond to `location` and `scale`, respectively.
#' * `"t"`: t distribution as parameterized in `stats::pt()` where `mean` and
#'   `extra_par` correspond to `ncp` and `df`, respectively.
#'
#' @details
#'
#' @return A distribution function that has a single argument. This is the
#'   vector of values in which the distribution function is evaluated.
cdf_fun = function(mean, extra_par, family){
  cdf_function = function(x){
    # Distribution function evaluated in x.
    cdf_x = switch(
      family,
      normal = stats::pnorm(q = x, mean = mean, sd = extra_par),
      logistic = stats::plogis(
        q = x,
        location = mean,
        scale = extra_par
      ),
      t = stats::pt(q = x, ncp = mean, df = extra_par)
    )
    return(cdf_x)
  }
  # Return the appropriate pdf function.
  return(cdf_function)
}

#' Function factory for density functions
#'
#' @inheritParams cdf_fun
#'
#' @return A density function that has a single argument. This is the vector of
#'   values in which the density function is evaluated.
pdf_fun = function(mean, extra_par, family){
  pdf_function = function(x){
    # Density function evaluated in x.
    pdf_x = switch(
      family,
      normal = stats::dnorm(x = x, mean = mean, sd = extra_par),
      logistic = stats::dlogis(
        x = x,
        location = mean,
        scale = extra_par
      ),
      t = stats::dt(x = x, ncp = mean, df = extra_par)
    )
    return(pdf_x)
  }
  # Return the appropriate pdf function.
  return(pdf_function)
}
