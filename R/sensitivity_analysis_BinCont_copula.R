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
#' @param n_sim Number of copula parameter vectors to be sampled.
#' @param eq_cond_association (boolean) Indicates whether \eqn{\rho_{13;2}} and
#' \eqn{\rho_{24;3}} are set equal.
#' @param lower (numeric) Vector of length 4 that provides the lower limit,
#'   \eqn{\boldsymbol{a} = (a_{23}, a_{13;2}, a_{24;3},
#'   a_{14;23})'}. Defaults to `c(-1, -1, -1, -1)`. If the provided lower limit
#'   is smaller than what is allowed for a particular copula family, then the
#'   copula family's lowest possible value is used instead.
#' @param upper (numeric) Vector of length 4 that provides the upper limit,
#'   \eqn{\boldsymbol{b} = (b_{23}, b_{13;2}, b_{24;3},
#'   b_{14;23})'}. Defaults to `c(1, 1, 1, 1)`.
#' @inheritParams sample_dvine
#'
#' @return A `n_sim` by `4` numeric matrix where each row corresponds to a
#'   sample for \eqn{\boldsymbol{\theta}_{unid}}.
sample_copula_parameters = function(copula_family2,
                                    n_sim,
                                    eq_cond_association = FALSE,
                                    lower = c(-1,-1,-1,-1),
                                    upper = c(1, 1, 1, 1)) {
  requireNamespace("copula")
  # If copula_family2 contains only 1 element, this vector is appended to
  # the correct length.
  if(length(copula_family2) == 1) copula_family2 = rep(copula_family2, 4)

  # Families that are restricted to positive associations.
  positive_only = c("clayton", "gumbel")
  all_associations = c("gaussian", "frank")
  # Enforce restrictions of copula family on lower and upper limits that were
  # provided by the user.
  lower = purrr::map2_dbl(
    .x = copula_family2,
    .y = lower,
    .f = function(x, y) {
      if (x %in% positive_only) {
        ifelse(y > 0, y, 0)
      }
      else if (x %in% all_associations) {
        ifelse(y > -1, y, -1)
      }
    }
  )
  # Regardless of the copula family, the upper bound should be (at max) 1. In
  # principle, there also exist copula families with an upper bound for
  # Spearman's rho strictly smaller than 1; however, they should probably not be
  # used in a sensitivity analysis (and have not been implemented).
  upper = ifelse(upper < 1, upper, 1)

  # Sample Spearman's rho parameters from the specified uniform distributions.
  # This code chunk returns a list.
  u = purrr::map2(
    .x = lower,
    .y = upper,
    .f = function(x, y)
      runif(n = n_sim, min = x, max = y)
  )
  # Impose conditional independence assumption if necessary.
  if (eq_cond_association) {
    u[[2]] = u[[3]]
  }

  # Convert sampled Spearman's rho parameter to the copula parameter scale.
  c = purrr::map2(
    .x = copula_family2,
    .y = u,
    .f = function(x, y) {
      switch(
        x,
        frank = {
          copula_fun = copula::frankCopula()
          upper_limit = 35
        },
        gaussian = {
          copula_fun = copula::ellipCopula(family = "normal")
          upper_limit = 1
        },
        clayton = {
          copula_fun = copula::claytonCopula()
          upper_limit = 28
        },
        gumbel = {
          copula_fun = copula::gumbelCopula()
          upper_limit = 50
        }
      )
      c_vec = purrr::map_dbl(
        .x = y,
        .f = function(x) copula::iRho(copula_fun, rho = x)
      )
      # Clayton copulas cannot handle independence (c_vec will contain NAs). The
      # ad-hoc solution employed here is to let the copula parameter be close
      # enough to the independence copula.
      if(x == "clayton") c_vec = ifelse(y == 0, 1e-5, c_vec)

      # Functions for the chosen copulas cannot handle parameter values larger
      # than upper_limit as defined above.
      c_vec = ifelse(c_vec > upper_limit | is.na(c_vec), upper_limit, c_vec)
      return(c_vec)
    }
  )

  # Convert list to a data frame.
  c = as.data.frame(c, col.names = c("c23", "c13_2", "c24_3", "c14_23"))

  return(c)
}

sample_rotation_parameters = function(n_sim, degrees = c(0, 90, 180, 270)) {
  r = sample(x = degrees, size = 4 * n_sim, replace = TRUE)
  r = matrix(r, ncol = 4)
  r = as.data.frame(r, col.names = c("r23", "r13_2", "r24_3", "r14_23"))
  return(r)
}

#' Perform Sensitivity Analysis for the Individual Causal Association with a
#' Continuous Surrogate and Binary True Endpoint

#' @details
#' # Information-Theoretic Causal Inference Framework
#'
#' The information-theoretic causal inference (ITCI) is a general framework to
#' evaluate surrogate endpoints in the single-trial setting (Alonso et al.,
#' 2015). In this framework, we focus on the individual causal effects,
#' \eqn{\Delta S = S_1 - S_0} and \eqn{\Delta T = T_1 - T_0} where \eqn{S_z}
#' and \eqn{T_z} are the potential surrogate end true endpoint under treatment
#' \eqn{Z = z}.
#'
#' In the ITCI framework, we say that \eqn{S} is a good surrogate for \eqn{T}
#' if
#' *\eqn{\Delta S} conveys a substantial amount of information on \eqn{\Delta T}*
#' (Alonso, 2018). This amount of shared information can generally be quantified
#' by the mutual information between \eqn{\Delta S} and \eqn{\Delta T},
#' denoted by \eqn{I(\Delta S; \Delta T)}. However, the mutual information lies
#' in \eqn{[0, + \infty]} which complicates the interpretation. In addition,
#' the mutual information may not be defined in specific scenarios where
#' absolute continuity of certain probability measures fails. Therefore, the
#' mutual information is transformed, and possibly modified, to enable a simple
#' interpretation in light of the definition of surrogacy. The resulting measure
#' is termed the individual causal association (ICA). This is explained in
#' the next sections.
#'
#' While the definition of surrogacy in the ITCI framework rests on information
#' theory, shared information is closely related to statistical association. Hence,
#' we can also define the ICA in terms of statistical association measures, like
#' Spearman's rho and Kendall's tau. The advantage of the latter are that they
#' are well-known, simple and rank-based measures of association.
#'
#' # Quantifying Surrogacy
#'
#' Alonso et al. (na) proposed to the following measure for the ICA: \deqn{R^2_H
#' = \frac{I(\Delta S; \Delta T)}{H(\Delta T)}} where \eqn{H(\Delta T)} is the
#' entropy of \eqn{\Delta T}. By token of that transformation of the mutual
#' information, \eqn{R^2_H} is restricted to the unit interval where 0 indicates
#' independence, and 1 a functional relationship between \eqn{\Delta S} and
#' \eqn{\Delta T}.
#'
#' The association between \eqn{\Delta S} and \eqn{\Delta T} can also be
#' quantified by Spearman's \eqn{\rho} (or Kendall's \eqn{\tau}). This quantity
#' requires appreciably less computing time than the mutual information. This
#' quantity is therefore always returned for every replication of the
#' sensitivity analysis.
#'
#' @inheritSection sensitivity_analysis_SurvSurv_copula Sensitivity Analysis
#' @inheritSection sensitivity_analysis_SurvSurv_copula Additional Assumptions
#'
#'
#' @param fitted_model Returned value from [fit_copula_model_BinCont()]. This object
#'   contains the estimated identifiable part of the joint distribution for the
#'   potential outcomes.
#' @inheritParams sensitivity_analysis_SurvSurv_copula
#' @inheritParams sample_copula_parameters
#'
#' @return A data frame is returned. Each row represents one replication in the
#'   sensitivity analysis. The returned data frame always contains the following
#'   columns:
#' * `R2H`, `sp_rho`, `minfo`: ICA as quantified by \eqn{R^2_H}, Spearman's rho, and
#'   Kendall's tau, respectively.
#' * `c12`, `c34`: estimated copula parameters.
#' * `c23`, `c13_2`, `c24_3`, `c14_23`: sampled copula parameters of the
#'   unidentifiable copulas in the D-vine copula. The parameters correspond to
#'   the parameterization of the `copula_family2` copula as in the `copula`
#'   R-package.
#' * `r12`, `r34`: Fixed rotation parameters for the two identifiable copulas.
#' * `r23`, `r13_2`, `r24_3`, `r14_23`: Sampled rotation parameters of the
#'   unidentifiable copulas in the D-vine copula. These values are constant for
#'   the Gaussian copula family since that copula is invariant to rotations.
#'
#' The returned data frame also contains the following columns when
#' `marg_association` is `TRUE`:
#' * `sp_s0s1`, `sp_s0t0`, `sp_s0t1`, `sp_s1t0`, `sp_s1t1`, `sp_t0t1`:
#' Spearman's rho between the corresponding potential outcomes. Note that these
#' associations refer to the observable potential outcomes. In contrary, the
#' estimated association parameters from [fit_copula_model_BinCont()] refer to
#' associations on a latent scale.
#'
#' @inherit fit_copula_model_BinCont examples
sensitivity_analysis_BinCont_copula = function(fitted_model,
                                               n_sim,
                                               eq_cond_association = TRUE,
                                               lower = c(-1, -1, -1, -1),
                                               upper = c(1, 1, 1, 1),
                                               marg_association = TRUE,
                                               n_prec = 1e4,
                                               ncores = 1) {
  if(!requireNamespace("fitdistrplus")) {
    errorCondition("Please install the 'fitdistrplus' package.")
  }
  # Extract relevant estimated parameters/objects for the fitted copula model.
  copula_family = fitted_model$copula_family
  copula_family2 = fitted_model$copula_family
  q_S0 = function(p) {
    q = quantile(fitted_model$submodel0$marginal_S_dist, probs = p)
    q = as.numeric(t(q$quantiles))
    return(q)
  }
  q_S1 = function(p) {
    q = quantile(fitted_model$submodel1$marginal_S_dist, probs = p)
    q = as.numeric(t(q$quantiles))
    return(q)
  }
  c12 = coef(fitted_model$submodel0$ml_fit)[length(coef(fitted_model$submodel0$ml_fit))]
  c34 = coef(fitted_model$submodel1$ml_fit)[length(coef(fitted_model$submodel0$ml_fit))]
  # For the Gaussian copula, fisher's Z transformation was applied. We have to
  # backtransform to the correlation scale in that case.
  if (copula_family == "gaussian") {
    c12 = (exp(2 * c12) - 1) / (exp(2 * c12) + 1)
    c34 = (exp(2 * c34) - 1) / (exp(2 * c34) + 1)
  }
  # Sample unidentifiable copula parameters.
  c = sample_copula_parameters(copula_family2 = copula_family2,
                               n_sim = n_sim,
                               eq_cond_association = eq_cond_association,
                               lower = lower,
                               upper = upper
                               )
  c = cbind(rep(c12, n_sim), c[, 1], rep(c34, n_sim), c[, 2:4])
  c_list = purrr::map(.x = split(c, seq(nrow(c))), .f = as.double)
  # Sample rotation parameters of the unidentifiable copula's family does not
  # allow for negative associations.
  if (copula_family2 %in% c("gumbel", "clayton")) {
    r = sample_rotation_parameters(n_sim)
  } else {
    r = sample_rotation_parameters(n_sim, degrees = 0)
  }
  # Add rotation parameters for identifiable copulas. Rotation parameters are
  # 180 because survival copulas were fitted.
  rotation_identifiable = 180
  # The Gaussian copula is invariant to rotations and a non-zero rotation
  # parameter for the Gaussian copula will give errors. The rotation parameters
  # are therefore set to zero for the Gaussian copula.
  if (copula_family %in% c("gaussian", "frank")) rotation_identifiable = 0
  r = cbind(rep(rotation_identifiable, n_sim),
            r[, 1],
            rep(rotation_identifiable, n_sim),
            r[, 2:4])
  r_list = purrr::map(.x = split(r, seq(nrow(r))), .f = as.double)
  # For every set of sampled unidentifiable parameters, compute the
  # required quantities.

  # Put all other arguments in a list for the apply function.
  MoreArgs = list(
    copula_family1 = copula_family,
    copula_family2 = copula_family2,
    n_prec = n_prec,
    q_S0 = q_S0,
    q_S1 = q_S1
  )

  # Use multicore computing if asked.
  if(ncores > 1 & requireNamespace("parallel")){
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
  rownames(measures_df) = NULL

  colnames(c) = c("c12", "c23", "c34", "c13_2", "c24_3", "c14_23")
  colnames(r) = c("r12", "r23", "r34", "r13_2", "r24_3", "r14_23")
  # Return a data frame with the results of the sensiviity analysis.
  return(cbind(as.data.frame(measures_df), c, r))

}
