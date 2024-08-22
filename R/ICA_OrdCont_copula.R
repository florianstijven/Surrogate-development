estimate_mutual_information_OrdCont = function(delta_S, delta_T) {

  requireNamespace("cubature")
  # Estimate conditional densities for all possible values of delta T.
  support_delta_T = unique(delta_T)

  lower_S = min(delta_S)
  upper_S = max(delta_S)
  range_S = upper_S - lower_S
  lower_S = lower_S - 0.2 * range_S
  upper_S = upper_S + 0.2 * range_S

  # Estimate the conditional densities for every possible value of Delta T and
  # immediately convert the estimated densities to R functions.
  dens_delta_S_given_delta_T =
    lapply(X = support_delta_T,
           FUN = function(x) {
             densfun = density(delta_S[delta_T == x],
                     from = lower_S,
                     to = upper_S,
                     n = 1024)
             approxfun(densfun$x, densfun$y)
           })

  # Compute marginal probabilities for distribution of Delta T.
  props_delta_T = sapply(
    X = support_delta_T,
    FUN = function(x) mean(delta_T == x)
  )

  # Construct marginal density function of S.
  densfun_marg = function(x) {
    mixture_components = mapply(
      cond_dens = dens_delta_S_given_delta_T,
      prop = props_delta_T,
      FUN = function(cond_dens, prop) cond_dens(x) * prop
    )
    sum(mixture_components)
  }

  # Compute integral for each conditional density.
  part1_integrals = mapply(
    cond_dens = dens_delta_S_given_delta_T,
    prop = props_delta_T,
    FUN = function(cond_dens, prop) {
      prop * cubature::cubintegrate(
        f = function(x) {
          y = cond_dens(x)
          ifelse(y == 0,
                 0,
                 y * log(y))
        },
        lower = lower_S,
        upper = upper_S
      )$integral
    }
  )

  # Compute differential entropy of Delta S
  diff_entropy = cubature::cubintegrate(
    f = function(x) {
      y = densfun_marg(x)
      ifelse(y == 0,
             0,
             y * log(y))
    },
    lower = lower_S,
    upper = upper_S
  )$integral

  mutual_information = sum(part1_integrals) - diff_entropy
  return(mutual_information)
}

#' Estimate ICA in Ordinal-Continuous Setting
#'
#' `estimate_ICA_OrdCont()` estimates the individual causal association (ICA)
#' for a sample of individual causal treatment effects with a continuous
#' surrogate and an ordinal true endpoint. The ICA in this setting is defined as
#' follows, \deqn{R^2_H = \frac{I(\Delta S; \Delta T)}{H(\Delta T)}} where
#' \eqn{I(\Delta S; \Delta T)} is the mutual information and \eqn{H(\Delta T)}
#' the entropy.
#'
#' @details
#' # Individual Causal Association
#'
#' Many association measures can operationalize the ICA. For each setting, we
#' consider one 'main' definition for the ICA which follows from the mutual
#' information.
#'
#' ## Continuous-Continuous
#'
#' The ICA is defined as the squared informational coefficient of correlation
#' (SICC or \eqn{R^2_H}), which is a transformation of the mutual information
#' to the unit interval: \deqn{R^2_h = 1 - e^{-2 \cdot I(\Delta S; \Delta T)}}
#' where 0 indicates independence, and 1 a functional relationship between
#' \eqn{\Delta S} and \eqn{\Delta T}. If \eqn{(\Delta S, \Delta T)'} is bivariate
#' normal, the ICA equals the Pearson correlation between \eqn{\Delta S} and
#' \eqn{\Delta T}.
#'
#' ## Ordinal-Continuous
#'
#' The ICA is defined as the following transformation of the mutual information:
#' \deqn{R^2_H = \frac{I(\Delta S; \Delta T)}{H(\Delta T)},}
#' where \eqn{I(\Delta S; \Delta T)} is the mutual information and \eqn{H(\Delta T)}
#' the entropy.
#'
#' ## Ordinal-Ordinal
#'
#' The ICA is defined as the following transformation of the mutual information:
#' \deqn{R^2_H = \frac{I(\Delta S; \Delta T)}{\min \{H(\Delta S), H(\Delta T) \}},}
#' where \eqn{I(\Delta S; \Delta T)} is the mutual information, and \eqn{H(\Delta S)}
#' and \eqn{H(\Delta T)} the entropy of \eqn{\Delta S} and \eqn{\Delta T},
#' respectively.
#'
#'
#' @param delta_S (numeric) Vector of individual causal treatment effects on the
#'   surrogate.
#' @param delta_T (integer) Vector of individual causal treatment effects on the true
#'   endpoint.
#'
#' @return (numeric) Estimated ICA
estimate_ICA_OrdCont = function(delta_S, delta_T) {
  # Compute marginal probabilities for distribution of Delta T.
  support_delta_T = unique(delta_T)
  props_delta_T = sapply(
    X = support_delta_T,
    FUN = function(x) mean(delta_T == x)
  )
  # Compute ICA
  ICA = estimate_mutual_information_OrdCont(delta_S, delta_T) /
    compute_entropy(props_delta_T)
  return(ICA)
}




#' Compute Individual Causal Association for a given D-vine copula model in the
#' Ordinal-Continuous Setting
#'
#' The [compute_ICA_OrdCont()] function computes the individual causal
#' association for a fully identified D-vine copula model in the setting with a
#' continuous surrogate endpoint and an ordinal true endpoint.
#'
#' @inheritParams compute_ICA_ContCont
#'
#' @inherit compute_ICA_ContCont return
compute_ICA_OrdCont = function(copula_par,
                               rotation_par,
                               copula_family1,
                               copula_family2 = copula_family1,
                               n_prec,
                               q_S0,
                               q_T0,
                               q_S1,
                               q_T1,
                               marginal_sp_rho = TRUE,
                               seed = 1,
                               mutinfo_estimator = NULL)
{
  if (is.null(mutinfo_estimator)) mutinfo_estimator = function(delta_S, delta_T) {
    ICA = estimate_ICA_OrdCont(delta_S, delta_T)
    return(-0.5 * log(1 - ICA))
  }
  # We can use the ICA function for the continuous-continuous setting with an
  # alternative mutual information estimator.
  compute_ICA_ContCont(copula_par = copula_par,
                       rotation_par = rotation_par,
                       copula_family1 = copula_family1,
                       copula_family2 = copula_family2,
                       n_prec = n_prec,
                       q_S0 = q_S0,
                       q_T0 = q_T0,
                       q_S1 = q_S1,
                       q_T1 = q_T1,
                       marginal_sp_rho = marginal_sp_rho,
                       seed = seed,
                       mutinfo_estimator = mutinfo_estimator,
                       plot_deltas = FALSE)
}
