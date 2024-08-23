estimate_mutual_information_OrdOrd = function(delta_S, delta_T) {
  # Estimate marginal pmf for Delta S
  support_delta_S = unique(delta_S)
  props_delta_S = sapply(
    X = support_delta_S,
    FUN = function(x) mean(delta_S == x)
  )
  pmf_delta_S = function(x) {
    sapply(
      X = x,
      FUN = function(x)
        props_delta_S[support_delta_S == x]
    )
  }
  # Estimate marginal pmf for Delta T
  support_delta_T = unique(delta_T)
  props_delta_T = sapply(
    X = support_delta_T,
    FUN = function(x) mean(delta_T == x)
  )
  pmf_delta_T = function(x) {
    sapply(
      X = x,
      FUN = function(x)
        props_delta_T[support_delta_T == x]
    )
  }
  # Estimated joint pmf
  joint_pmf = function(x, y) {
    pmf = mapply(x = x, y = y, function(x, y) {
      mean((delta_S == x) & (delta_T == y))
    })
    return(as.numeric(pmf))
  }


  # Compute entropy of joint distribution
  x = expand.grid(support_delta_T, support_delta_S)[, 1]
  y = expand.grid(support_delta_T, support_delta_S)[, 2]
  integrand = ifelse(
    joint_pmf(x, y) > 0,
    -1 * log(joint_pmf(x, y)) * joint_pmf(x, y),
    0
  )
  joint_entropy = sum(as.numeric(integrand))
  # Compute entropy of Delta S
  integrand = ifelse(
    pmf_delta_S(support_delta_S) > 0,
    -1 * log(pmf_delta_S(support_delta_S)) * pmf_delta_S(support_delta_S),
    0
  )
  entropy_delta_S = sum(integrand)
  # Compute entropy of Delta T
  integrand = ifelse(
    pmf_delta_T(support_delta_T) > 0,
    -1 * log(pmf_delta_T(support_delta_T)) * pmf_delta_T(support_delta_T),
    0
  )
  entropy_delta_T = sum(integrand)

  mutual_information = entropy_delta_T + entropy_delta_S - joint_entropy
  return(mutual_information)
}

#' Estimate ICA in Ordinal-Ordinal Setting
#'
#' `estimate_ICA_OrdOrd()` estimates the individual causal association (ICA) for
#' a sample of individual causal treatment effects with an ordinal surrogate and
#' true endpoint. The ICA in this setting is defined as follows: \deqn{R^2_H =
#' \frac{I(\Delta S; \Delta T)}{\min \{H(\Delta S), H(\Delta T) \}}} where
#' \eqn{I(\Delta S; \Delta T)} is the mutual information, and \eqn{H(\Delta S)}
#' and \eqn{H(\Delta T)} the entropy of \eqn{\Delta S} and \eqn{\Delta T},
#' respectively.
#'
#' @param delta_S (integer) Vector of individual causal treatment effects on the
#'   surrogate.
#' @param delta_T (integer) Vector of individual causal treatment effects on the
#'   true endpoint.
#'
#' @return (numeric) Estimated ICA
estimate_ICA_OrdOrd = function(delta_S, delta_T) {
  # Compute marginal probabilities for distribution of Delta T.
  support_delta_T = unique(delta_T)
  props_delta_T = sapply(
    X = support_delta_T,
    FUN = function(x) mean(delta_T == x)
  )
  # Compute marginal probabilities for distribution of Delta S.
  support_delta_S = unique(delta_S)
  props_delta_S = sapply(
    X = support_delta_S,
    FUN = function(x) mean(delta_S == x)
  )
  # Compute ICA
  ICA = estimate_mutual_information_OrdOrd(delta_S, delta_T) /
    min(c(
      compute_entropy(props_delta_T),
      compute_entropy(props_delta_S)
    ))
  return(ICA)
}

estimate_entropy = function(x) {
  support = unique(x)
  props_x = sapply(
    X = support,
    FUN = function(y) mean(x == y)
  )
  compute_entropy(props_x)
}




#' Compute Individual Causal Association for a given D-vine copula model in the
#' Ordinal-Ordinal Setting
#'
#' The [compute_ICA_OrdOrd()] function computes the individual causal
#' association for a fully identified D-vine copula model in the setting with an
#' ordinal surrogate and true endpoint.
#'
#' @param ICA_estimator Function that estimates the ICA between the first two
#'   arguments which are numeric vectors. Defaults to `NULL` which corresponds
#'   to using [estimate_ICA_OrdOrd()].
#' @inheritParams compute_ICA_ContCont
#'
#' @inherit compute_ICA_ContCont return
compute_ICA_OrdOrd = function(copula_par,
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
                              ICA_estimator = NULL)
{
  if (is.null(ICA_estimator)) {
    ICA_estimator = estimate_ICA_OrdOrd
  }

  # We can use the ICA function for the continuous-continuous setting with an
  # alternative mutual information estimator.
  compute_ICA_ContCont(
    copula_par = copula_par,
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
    ICA_estimator = ICA_estimator,
    plot_deltas = FALSE
  )
}
