#' Computes loglikelihood for a given copula model
#'
#' `log_likelihood_copula_model()` computes the loglikelihood for a given
#' bivariate copula model and data set while allowin for right-censoring of both
#' outcome variables.
#'
#' @param X Numeric vector corresponding to first outcome variable.
#' @param Y Numeric vector corresponding to second outcome variable.
#' @param cdf_X
#' @param cdf_Y
#' @param pdf_X
#' @param pdf_Y
#' @inheritParams loglik_copula_scale
#'
#' @return
#' @export
#'
#' @examples
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
