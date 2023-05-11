compute_entropy = function(probs) {
  return(-1 * sum(probs * log(probs)))
}

compute_mutual_information_BinCont = function(delta_S, delta_T) {
  # Estimate three conditional densities for the three possible values of delta
  # T.
  lower_S = min(delta_S)
  upper_S = max(delta_S)
  range_S = upper_S - lower_S
  lower_S = lower_S - 0.2 * range_S
  upper_S = upper_S + 0.2 * range_S

  dens_min1 = density(delta_S[delta_T == -1],
                      from = lower_S,
                      to = upper_S,
                      n = 1024)
  dens_0 = density(delta_S[delta_T == 0],
                   from = lower_S,
                   to = upper_S,
                   n = 1024)
  dens_plus1 = density(delta_S[delta_T == 1],
                       from = lower_S,
                       to = upper_S,
                       n = 1024)

  # Convert the estimated densities to R functions.
  densfun_min1 = approxfun(dens_min1$x, dens_min1$y)
  densfun_0 = approxfun(dens_0$x, dens_0$y)
  densfun_plus1 = approxfun(dens_plus1$x, dens_plus1$y)

  # Compute marginal probabilities for distribution of Delta T.
  pi_min1 = mean(delta_T == -1)
  pi_0 = mean(delta_T == 0)
  pi_plus1 = 1 - pi_min1 - pi_0

  # Construct marginal density function of S.
  densfun_marg = function(x) {
    pi_min1 * densfun_min1(x) +
      pi_0 * densfun_0(x) +
      pi_plus1 * densfun_plus1(x)
  }

  # Compute integral for each conditional density using the trapezoidal rule.
  integral_min1 = pi_min1 * cubature::cubintegrate(
    f = function(x) {
      y = densfun_min1(x)
      ifelse(y == 0,
             0,
             y * log(y))
    },
    lower = lower_S,
    upper = upper_S
  )$integral
  integral_0 = pi_0 * cubature::cubintegrate(
    f = function(x) {
      y = densfun_0(x)
      ifelse(y == 0,
             0,
             y * log(y))
    },
    lower = lower_S,
    upper = upper_S
  )$integral
  integral_plus1 = pi_plus1 * cubature::cubintegrate(
    f = function(x) {
      y = densfun_plus1(x)
      ifelse(y == 0,
             0,
             y * log(y))
    },
    lower = lower_S,
    upper = upper_S
  )$integral

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

  mutual_information = integral_min1 + integral_0 + integral_plus1 - diff_entropy
  return(mutual_information)
}

compute_ICA_BinCont = function(delta_S, delta_T) {
  # Compute marginal probabilities for distribution of Delta T.
  pi_min1 = mean(delta_T == -1)
  pi_0 = mean(delta_T == 0)
  pi_plus1 = 1 - pi_min1 - pi_0
  # Compute ICA
  ICA = compute_mutual_information_BinCont(delta_S, delta_T) /
    compute_entropy(c(pi_min1, pi_0, pi_plus1))
  return(ICA)
}


#' Sample copula data from a given four-dimensional D-vine copula
#'
#' @param copula_par
#' @param rotation_par
#' @param copula_family1
#' @param copula_family2
#' @param n
#'
#' @return
#' @export
#'
#' @examples
sample_dvine = function(copula_par,
                        rotation_par,
                        copula_family1,
                        copula_family2 = copula_family1,
                        n) {
  # D-vine copula parameters
  c12 = copula_par[1]
  c34 = copula_par[2]
  c23 = copula_par[3]
  c13_2 = copula_par[4]
  c24_3 = copula_par[5]
  c14_23 = copula_par[6]

  # D-vine copula rotations
  r12 = rotation_par[1]
  r34 = rotation_par[2]
  r23 = rotation_par[3]
  r13_2 = rotation_par[4]
  r24_3 = rotation_par[5]
  r14_23 = rotation_par[6]

  pair_copulas = list(
    list(
      rvinecopulib::bicop_dist(
        family = fitted_model$copula_family,
        rotation = r12,
        parameters = c12
      ),
      rvinecopulib::bicop_dist(
        family = copula_family2,
        rotation = r23,
        parameters = c23
      ),
      rvinecopulib::bicop_dist(
        family = fitted_model$copula_family,
        rotation = r34,
        parameters = c34
      )
    ),
    list(
      rvinecopulib::bicop_dist(
        family = copula_family2,
        rotation = r13_2,
        parameters = c13_2
      ),
      rvinecopulib::bicop_dist(
        family = copula_family2,
        rotation = r24_3,
        parameters = c24_3
      )
    ),
    list(
      rvinecopulib::bicop_dist(
        family = copula_family2,
        rotation = r14_23,
        parameters = c14_23
      )
    )
  )
  copula_structure = rvinecopulib::dvine_structure(order = 1:4)
  vine_cop = rvinecopulib::vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)

  U = rvinecopulib::rvinecop(n = n, vinecop = vine_cop, cores = 1)


  return(U)
}
