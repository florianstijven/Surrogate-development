#' Sample Unidentifiable Copula Parameters
#'
#' The [sample_copula_parameters()] function samples the unidentifiable copula
#' parameters for the partly identifiable D-vine copula model, see for example
#' [fit_copula_model_BinCont()] and [fit_model_SurvSurv()] for more information
#' regarding the D-vine copula model.
#'
#' @details # Sampling
#'
#'   In the D-vine copula model in the Information-Theoretic Causal Inference
#'   (ITCI) framework, the following copulas are not identifiable: \eqn{c_{23}},
#'   \eqn{c_{13;2}}, \eqn{c_{24;3}}, \eqn{c_{14;23}}. Let the corresponding
#'   copula
#' parameters be \deqn{\boldsymbol{\theta}_{unid} = (\theta_{23}, \theta_{13;2},
#' \theta_{24;3}, \theta_{14;23})'.}
#'   The allowable range for this parameter vector depends on the corresponding
#'   copula families. For parsimony and comparability across different copula
#'   families, the sampling procedure consists of two steps:
#'  1. Sample Spearman's rho parameters from a uniform distribution,
#'   \deqn{\boldsymbol{\rho}_{unid} = (\rho_{23}, \rho_{13;2}, \rho_{24;3},
#'   \rho_{14;23})' \sim U(\boldsymbol{a}, \boldsymbol{b}).}
#'  1. Transform the sampled Spearman's rho parameters to the copula parameter
#'  scale, \eqn{\boldsymbol{\theta}_{unid}}.
#'
#'  These two steps are repeated `n_sim` times.
#'
#'   # Conditional Independence
#'
#' In addition to range restrictions through the `lower` and `upper` arguments,
#' we allow for so-called conditional independence assumptions.
#' These assumptions entail that \eqn{\rho_{13;2} = 0} and \eqn{\rho_{24;3} =
#' 0}. Or in other words, \eqn{U_1 \perp U_3 \, | \, U_2} and \eqn{U_2 \perp U_4 \, | \, U_3}.
#' In the context of a surrogate evaluation trial (where \eqn{(U_1, U_2, U_3,
#' U_4)'} corresponds to the probability integral transformation of \eqn{(T_0,
#' S_0, S_1, T_1)'}) this assumption could be justified by subject-matter knowledge.
#'
#'
#' @param copula_family2 Copula family of the unidentifiable bivariate copulas.
#'   For the possible options, see [loglik_copula_scale()].
#' @param n_sim Number of copula parameter vectors to be sampled.
#' @param cond_ind (boolean) Indicates whether conditional independence is
#'   assumed, see Conditional Independence. Defaults to `FALSE`.
#' @param lower (numeric) Vector of length 4 that provides the lower limit,
#'   \eqn{\boldsymbol{a} = (a_{23}, a_{13;2}, a_{24;3},
#'   a_{14;23})'}. Defaults to `c(-1, -1, -1, -1)`. If the provided lower limit
#'   is smaller than what is allowed for a particular copula family, then the
#'   copula family's lowest possible value is used instead.
#' @param upper (numeric) Vector of length 4 that provides the upper limit,
#'   \eqn{\boldsymbol{b} = (b_{23}, b_{13;2}, b_{24;3},
#'   b_{14;23})'}. Defaults to `c(1, 1, 1, 1)`.
#'
#' @return A `n_sim` by `4` numeric matrix where each row corresponds to a
#'   sample for \eqn{\boldsymbol{\theta}_{unid}}.
sample_copula_parameters = function(copula_family2,
                                    n_sim,
                                    cond_ind = FALSE,
                                    lower = c(-1,-1,-1,-1),
                                    upper = c(1, 1, 1, 1)) {
  # Enforce restrictions of copula family on lower and upper limits that were
  # provided by the user.
  switch(
    copula_family2,
    frank = {
      lower = ifelse(lower > -1, lower, -1)
      upper = ifelse(upper < 1, lower, 1)
    },
    gaussian = {
      lower = ifelse(lower > -1, lower, -1)
      upper = ifelse(upper < 1, lower, 1)
    },
    clayton = {
      lower = ifelse(lower > 0, lower, 0)
      upper = ifelse(upper < 1, lower, 1)
    },
    gumbel = {
      lower = ifelse(lower > 0, lower, 0)
      upper = ifelse(upper < 1, lower, 1)
    }
  )
  # Sample Spearman's rho parameters from the specified uniform distributions.
  # This code chunk returns a list.
  u = purrr::map2(
    .x = lower,
    .y = upper,
    .f = function(x, y)
      runif(n = n_sim, min = x, max = y)
  )
  # Impose conditional independence assumption if necessary.
  if (cond_ind) {
    u[[2]] = rep(0, n_sim)
    u[[3]] = rep(0, n_sim)
  }
  # Convert sampled Spearman's rho parameter to the copula parameter scale.
  switch(
    copula_family2,
    frank = {
      c = purrr::map(
        .x = u,
        .f = function(x) {
          purrr::map_dbl(
            .x = x,
            .f = function(x) {
              c_vec = copula::iRho(copula::frankCopula(), rho = x)
              # Functions for the frank copula cannot handle parameters larger
              # than abs(35).
              ifelse(abs(c_vec) > 35 | is.na(c_vec), sign(c_vec) * 35, c_vec)
            }
          )
        }
      )
    },
    gaussian = {
      c = purrr::map(
        .x = u,
        .f = function(x) {
          purrr::map_dbl(
            .x = x,
            .f = function(x) {
              copula::iRho(copula::ellipCopula(family = "normal"), rho = x)
            }
          )
        }
      )
    },
    clayton = {
      c = purrr::map(
        .x = u,
        .f = function(x) {
          purrr::map_dbl(
            .x = x,
            .f = function(x) {
              c_vec = ifelse(x != 0, copula::iRho(copula::claytonCopula(), rho = x), 1e-5)
              # Functions for the Clayton copula cannot handle parameter values
              # larger than 28.
              c_vec = ifelse(c_vec > 28 | is.na(c_vec), 28, c_vec)
            }
          )
        }
      )
    },
    gumbel = {
      c = purrr::map(
        .x = u,
        .f = function(x) {
          purrr::map_dbl(
            .x = x,
            .f = function(x) {
              c_vec = copula::iRho(copula::gumbelCopula(), rho = x)
              # Functions for the Gumbel copula cannot handle parameter values
              # larger than 50.
              c_vec = ifelse(c_vec > 50 | is.na(c_vec), 50, c_vec)
            }
          )
        }
      )
    }
  )
  # Convert list to a data frame.
  c = as.data.frame(c, col.names = c("theta_23", "theta_13;2", "theta_24;3", "theta_14;23"))

  return(c)
}

sample_rotation_parameters = function(n_sim, degrees = c(0, 90, 180, 270)) {
  r = sample(x = degrees, size = 4 * n_sim, replace = TRUE)
  r = matrix(r, ncol = 4)
  r = as.data.frame(r, col.names = c("r_23", "r_13;2", "r_24;3", "r_14;23"))
  return(r)
}

sensitivity_analysis_BinCont_copula = function(fitted_model,
                                               n_sim,
                                               cond_ind,
                                               lower = c(-1, -1, -1, -1),
                                               upper = c(1, 1, 1, 1),
                                               minfo_prec = 1e4) {
  # Extract relevant estimated parameters/objects for the fitted copula model.
  copula_family = fitted_model$copula_family
  copula_family2 = fitted_model$copula_family
  q_S0 = function(p) {
    q = quantile(fitted_model$submodel0$marginal_S_dist, probs = p)
    q = as.numeric(t(q$quantiles))
  }
  q_S1 = function(p) {
    q = quantile(fitted_model$submodel1$marginal_S_dist, probs = p)
    q = as.numeric(t(q$quantiles))
  }
  c_12 = coef(fitted_model$submodel0$ml_fit)[length(coef(fitted_model$submodel0$ml_fit))]
  c_34 = coef(fitted_model$submodel1$ml_fit)[length(coef(fitted_model$submodel0$ml_fit))]
  # Sample unidentifiable copula parameters.
  c = sample_copula_parameters(copula_family2 = copula_family2,
                               n_sim = n_sim,
                               cond_ind = cond_ind,
                               lower = lower,
                               upper = upper
                               )
  c = cbind(rep(c_12, n_sim), rep(c_34, n_sim), c)
  c_list = purrr::map(.x = split(c, seq(nrow(c))), .f = as.double)
  # Sample rotation parameters of the unidentifiable copula's family does not
  # allow for negative associations.
  if (copula_family2 %in% c("gumbel", "clayton")) {
    r = sample_rotation_parameters(n_sim)
  } else {
    r = sample_rotation_parameters(n_sim, degrees = 0)
  }
  # Add rotation parameters for identifiable copulas.
  r = cbind(rep(0, n_sim), rep(0, n_sim), r)
  r_list = purrr::map(.x = split(r, seq(nrow(r))), .f = as.double)
  # For every set of sampled unidentifiable parameters, compute the
  # required quantities.

  # Put all other arguments in a list for the apply function.
  MoreArgs = list(
    copula_family1 = copula_family,
    copula_family2 = copula_family2,
    minfo_prec = minfo_prec,
    q_S0 = q_S0,
    q_S1 = q_S1
  )

  # Use multicore computing if asked.
  if(ncores > 1){
    cl  <- parallel::makeCluster(ncores)
    print("Starting parallel simulations")
    temp = parallel::clusterMap(
      cl = cl,
      fun = compute_ICA_BinCont,
      copula_par = c_list,
      rotation_par = r_list,
      seed = 1:n_sim,
      MoreArgs = MoreArgs
    )
    on.exit(expr = {parallel::stopCluster(cl)
      rm("cl")})
  }
  else if (ncores == 1){
    temp = mapply(
      FUN = compute_ICA_BinCont,
      c = c_list,
      r = r_list,
      MoreArgs = MoreArgs,
      SIMPLIFY = FALSE
    )
  }

  measures_df = t(as.data.frame(temp))
  if(get_marg_tau){
    colnames(measures_df) = c("kendall", "sp_rho", "minfo",
                              "sp_s0s1", "sp_s0t0", "sp_s0t1",
                              "sp_s1t0", "sp_s1t1",
                              "sp_t0t1",
                              "prop_harmed", "prop_protected",
                              "prop_always", "prop_never")
  }
  else{
    colnames(measures_df) = c("kendall", "sp_rho", "minfo")
  }
  rownames(measures_df) = NULL

  c = as.data.frame(x = c_matrix)
  colnames(c) = c("c23", "c13_2", "c24_3", "c14_23")
  r = as.data.frame(x = r_matrix)
  colnames(r) = c("r23", "r13_2", "r24_3", "r14_23")
  return(dplyr::bind_cols(as.data.frame(measures_df), c, r))

  # Return a data frame with the results of the sensiviity analysis.

}
