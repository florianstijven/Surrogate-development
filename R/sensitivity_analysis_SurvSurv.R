#' Sensitivity analysis for individual causal association
#'
#' The [sensitivity_analysis_SurvSurv_copula()] function performs the sensitivity analysis
#' for the individual causal association (ICA) as described by Stijven et
#' al. (2022).
#'
#' @details
#'
#' # Quantifying Surrogacy
#'
#' In the causal-inference framework to evaluate surrogate endpoints, the ICA is
#' the measure of primary interest. This measure quantifies the association
#' between the individual causal treatment effects on the surrogate (\eqn{\Delta
#' S}) and on the true endpoint (\eqn{\Delta T}). Stijven et al. (2022) proposed
#' to quantify this association through the squared informational coefficient of
#' correlation (SICC or \eqn{R^2_H}), which is based on information-theoretic
#' principles. Indeed, \eqn{R^2_H} is a transformation of the mutual information
#' between \eqn{\Delta S} and \eqn{\Delta T}, \deqn{R^2_H = 1 - e^{-2 \cdot
#' I(\Delta S; \Delta T)}.} By token of that transformation, \eqn{R^2_H} is
#' restricted to the unit interval where 0 indicates independence, and 1 a
#' functional relationship between \eqn{\Delta S} and \eqn{\Delta T}. The mutual
#' information is returned by [ica_SurvSurv_sens()] if a non-zero value is
#' specified for `minfo_prec` (see Arguments).
#'
#' The association between \eqn{\Delta S} and \eqn{\Delta T} can also be
#' quantified by Spearman's \eqn{\rho} (or Kendall's \eqn{\tau}). This quantity
#' requires appreciably less computing time than the mutual information. This
#' quantity is therefore always returned for every replication of the
#' sensitivity analysis.
#'
#' # Sensitivity Analysis
#'
#' Because \eqn{S_0} and \eqn{S_1} are never simultaneously observed in the same
#' patient, \eqn{\Delta S} is not observable, and analogously for \eqn{\Delta
#' T}. Consequently, the ICA is unidentifiable. This is solved by considering a
#' (partly identifiable) model for the full vector of potential outcomes,
#' \eqn{(T_0, S_0, S_1, T_1)'}. The identifiable parameters are estimated. The
#' unidentifiable parameters are sampled from their parameters space in each
#' replication of a sensitivity analysis. If the number of replications
#' (`n_sim`) is sufficiently large, the entire parameter space for the
#' unidentifiable parameters will be explored/sampled. In each replication, all
#' model parameters are "known" (either estimated or sampled). Consequently, the
#' ICA can be computed in each replication of the sensitivity analysis.
#'
#' The sensitivity analysis thus results in a set of values for the ICA. This
#' set can be interpreted as *all values for the ICA that are compatible with
#' the observed data*. However, the range of this set is often quite broad; this
#' means there remains too much uncertainty to make judgements regarding the
#' worth of the surrogate. To address this unwieldy uncertainty, additional
#' assumptions can be used that restrict the parameter space of the
#' unidentifiable parameters. This in turn reduces the uncertainty regarding the
#' ICA.
#'
#' # Additional Assumptions
#'
#' There are two possible types of assumptions that restrict the parameter space
#' of the unidentifiable parameters: (i) *equality* type of assumptions, and
#' (ii) *inequality* type of assumptions. These are discussed in turn in the
#' next two paragraphs.
#'
#' The equality assumptions have to be incorporated into the sensitivity
#' analysis itself. Only one type of equality assumption has been implemented;
#' this is the conditional independence assumption which can be specified
#' through the `cond_ind` argument:
#' \deqn{\tilde{S}_0 \perp \!\!\! \perp T_1 | \tilde{S}_1 \; \text{and} \;
#' \tilde{S}_1 \perp \!\!\! \perp T_0 | \tilde{S}_0 .} This can informally be
#' interpreted as ``what the control treatment does to the surrogate does not
#' provide information on the true endpoint under experimental treatment if we
#' already know what the experimental treatment does to the surrogate", and
#' analogously when control and experimental treatment are interchanged. Note
#' that \eqn{\tilde{S}_z} refers to either the actual potential surrogate
#' outcome, or a latent version. This depends on the content of `fitted_model`.
#'
#' The inequality type of assumptions have to be imposed on the data frame that
#' is returned by the current function; those assumptions are thus imposed
#' *after* running the sensitivity analysis. If `marginal_association` is set to
#' `TRUE`, the returned data frame contains additional unverifiable quantities
#' that differ across replications of the sensitivity analysis: (i) the
#' unconditional Spearman's \eqn{\rho} for all pairs of (observable/non-latent)
#' potential outcomes, and (ii) the proportions of the population strata as
#' defined by Nevo and Gorfine (2022) if semi-competing risks are present. More
#' details on the interpretation and use of these assumptions can be found in
#' Stijven et al. (2022).
#'
#'
#' @param fitted_model Returned value from [fit_model_SurvSurv()]. This object
#'   contains the estimated identifiable part of the joint distribution for the
#'   potential outcomes.
#' @param n_sim Number of replications in the *sensitivity analysis*. This value
#'   should be large enough to sufficiently explore all possible values of the
#'   ICA. The minimally sufficient number depends to a large extent on which
#'   inequality assumptions are subsequently imposed (see Additional
#'   Assumptions).
#' @param n_prec Number of Monte-Carlo samples for the *numerical approximation*
#'   of the ICA in each replication of the sensitivity analysis.
#' @param minfo_prec Number of quasi Monte-Carlo samples for the numerical
#'   integration to obtain the mutual information. If this value is 0 (default),
#'   the mutual information is not computed and `NA` is returned for that column.
#' @param restr Default value should not be modified by the user.
#' @param copula_family2 Parametric family of the unidentifiable copulas in the
#'   D-vine copula. One of the following parametric copula families:
#'   `"clayton"`, `"frank"`, `"gaussian"`, or `"gumbel"`.
#' @param ncores Number of cores used in the sensitivity analysis. The
#'   computations are computationally heavy, and this option can speed things
#'   up considerably.
#' @param get_marg_tau Boolean.
#' * `TRUE`: Return marginal association measures
#'   in each replication in terms of Spearman's rho. The proportion of harmed,
#'   protected, never diseased, and always diseased is also returned. See also
#'   Value.
#' * `FALSE` (default): No additional measures are returned.
#' @param cond_ind Boolean.
#' * `TRUE`: Assume conditional independence (see Additional Assumptions).
#' * `FALSE` (default): Conditional independence is not assumed.
#'
#' @return A data frame is returned. Each row represents one replication in the
#'   sensitivity analysis. The returned data frame always contains the following
#'   columns:
#' * `kendall`, `sp_rho`, `minfo`: ICA as quantified by Kendall's \eqn{\tau},
#' Spearman's \eqn{\rho}, and the mutual information, respectively.
#' * `c23`, `c13_2`, `c24_3`, `c14_23`: sampled copula parameters of the
#' unidentifiable copulas in the D-vine copula. The parameters correspond to the
#' parameterization of the `copula_family2` copula as in the `copula` R-package.
#' * `r23`, `r13_2`, `r24_3`, `r14_23`: sampled rotation parameters of the
#' unidentifiable copulas in the D-vine copula. These values are constant for
#' the Gaussian copula family since that copula is invariant to rotations.
#'
#' The returned data frame also contains the following columns when `get_marg_tau`
#' is `TRUE`:
#' * `sp_s0s1`, `sp_s0t0`, `sp_s0t1`, `sp_s1t0`, `sp_s1t1`, `sp_t0t1`:
#' Spearman's \eqn{\rho} between the corresponding potential outcomes. Note
#' that these associations refer to the potential time-to-composite events
#' and/or time-to-true endpoint event. In contrary, the estimated association
#' parameters from [fit_model_SurvSurv()] refer to associations between the
#' time-to-surrogate event and time-to true endpoint event.
#' * `prop_harmed`, `prop_protected`, `prop_always`, `prop_never`: proportions
#' of the corresponding population strata in each replication. These are defined
#' in Nevo and Gorfine (2022).
#' @export
#'
#' @references Stijven, F., Alonso, a., Molenberghs, G., Van Der Elst, W., Van
#'   Keilegom, I. (2022). An information-theoretic approach to the evaluation of
#'   time-to-event surrogates for time-to-event true endpoints based on causal
#'   inference.
#'
#'   Nevo, D., & Gorfine, M. (2022). Causal inference for semi-competing risks
#'   data. Biostatistics, 23 (4), 1115-1132
#'
#' @examples
#' library(Surrogate)
#' data("Ovarian")
#' # For simplicity, data is not recoded to semi-competing risks format, but the
#' # data are left in the composite event format.
#' data = data.frame(
#'   Ovarian$Pfs,
#'   Ovarian$Surv,
#'   Ovarian$Treat,
#'   Ovarian$PfsInd,
#'   Ovarian$SurvInd
#' )
#' ovarian_fitted =
#'     fit_model_SurvSurv(data = data,
#'                        copula_family = "clayton",
#'                        nknots = 1)
#' # Illustration with small number of replications and low precision
#' ica_SurvSurv_sens(ovarian_fitted,
#'                   n_sim = 5,
#'                   n_prec = 2000,
#'                   copula_family2 = "clayton")
#'
#'
sensitivity_analysis_SurvSurv_copula = function(fitted_model,
                                                composite = TRUE,
                                                n_sim,
                                                cond_ind,
                                                lower = c(-1,-1,-1,-1),
                                                upper = c(1, 1, 1, 1),
                                                marg_association = TRUE,
                                                n_prec = 5e3,
                                                minfo_prec = 0,
                                                ncores = 1

){
  # Extract relevant estimated parameters/objects for the fitted copula model.
  copula_family = fitted_model$copula_family
  copula_family2 = fitted_model$copula_family
  k = length(fitted_model$knots0)
  q_S0 = function(p) {
    flexsurv::qsurvspline(
      p = p,
      gamma = coef(fitted_model$fit_0)[1:k],
      knots = fitted_model$knots0
    )
  }
  q_S1 = function(p) {
    flexsurv::qsurvspline(
      p = p,
      gamma = coef(fitted_model$fit_1)[1:k],
      knots = fitted_model$knots1
    )
  }
  q_T0 = function(p) {
    flexsurv::qsurvspline(
      p = p,
      gamma = coef(fitted_model$fit_0)[(k + 1):(2 * k)],
      knots = fitted_model$knott0
    )
  }
  q_T1 = function(p) {
    flexsurv::qsurvspline(
      p = p,
      gamma = coef(fitted_model$fit_1)[(k + 1):(2 * k)],
      knots = fitted_model$knott1
    )
  }

  #number of parameters
  n_par = length(coef(fitted_model$fit_0))

  #initialize vectors to store measures of surrogacy in
  kendall <- sp_rho <- minfo <- 1:n_sim

  #pull association parameters from estimated parameter vectors
  #the association parameter is always the last one in the corresponding vector
  c_12 = coef(fitted_model$fit_0)[n_par]
  c_34 = coef(fitted_model$fit_1)[n_par]

  c = sample_copula_parameters(
    copula_family2 = copula_family2,
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
  # Add rotation parameters for identifiable copulas. Rotation parameters are
  # 180 because survival copulas were fitted.
  r = cbind(rep(180, n_sim), rep(180, n_sim), r)
  r_list = purrr::map(.x = split(r, seq(nrow(r))), .f = as.double)
  # For every set of sampled unidentifiable parameters, compute the
  # required quantities.

  #put all other arguments in a list for the apply function
  MoreArgs = list(
    copula_family1 = copula_family,
    copula_family2 = copula_family2,
    n_prec = n_prec,
    minfo_prec = minfo_prec,
    q_S0,
    q_T0,
    q_S1,
    q_T1,
    composite = composite,
    marginal_sp_rho = marg_association,
    seed = 1
  )
  if (ncores > 1) {
    cl  <- parallel::makeCluster(ncores)
    on.exit(expr = {
      parallel::stopCluster(cl)
      rm("cl")
    })
    #helper function
    # surrogacy_sample_sens <- surrogacy_sample_sens
    print("Starting parallel simulations")
    # surrogacy_sample_sens
    # parallel::clusterExport(cl = cl, varlist = "surrogacy_sample_sens", )
    # parallel::clusterEvalQ(cl = cl, expr = library(flexsurv))
    # parallel::clusterEvalQ(cl = cl, expr = library(rvinecopulib))
    temp = parallel::clusterMap(
      cl = cl,
      fun = compute_ICA_SurvSurv,
      c = c_list,
      r = r_list,
      MoreArgs = MoreArgs
    )
  }
  else if (ncores == 1){
    temp = mapply(FUN = compute_ICA_SurvSurv,
                  copula_par = c_list, rotation_par = r_list, MoreArgs = MoreArgs,
                  SIMPLIFY = FALSE)
  }

  measures_df = t(as.data.frame(temp))
  if(marg_association){
    colnames(measures_df) = c("ICA", "sp_rho",
                              "sp_s0s1", "sp_s0t0", "sp_s0t1",
                              "sp_s1t0", "sp_s1t1",
                              "sp_t0t1",
                              "prop_harmed", "prop_protected",
                              "prop_always", "prop_never")
  }
  else{
    colnames(measures_df) = c("ICA", "sp_rho")
  }
  rownames(measures_df) = NULL

  colnames(c) = c("c_12", "c_34", "c23", "c13_2", "c24_3", "c14_23")
  colnames(r) = c("r_12", "r_34", "r23", "r13_2", "r24_3", "r14_23")
  return(dplyr::bind_cols(as.data.frame(measures_df), c, r))

}



#
# surrogacy_sample_sens = function(fitted_model,
#                                  n_prec,
#                                  c_unid,
#                                  r_unid,
#                                  restr = TRUE,
#                                  copula_family2,
#                                  get_marg_tau = FALSE
# ){
#   #number of knots k
#   n_k = (length(coef(fitted_model$fit_0)) - 1)/2
#   if (fitted_model$copula_family == "gaussian"){
#     c12 = (exp(coef(fitted_model$fit_0)[n_k*2 + 1]) - 1)/(exp(coef(fitted_model$fit_0)[n_k*2 + 1]) + 1)
#     c34 = (exp(coef(fitted_model$fit_1)[n_k*2 + 1]) - 1)/(exp(coef(fitted_model$fit_1)[n_k*2 + 1]) + 1)
#   }
#   else{
#     c12 = coef(fitted_model$fit_0)[n_k*2 + 1]
#     c34 = coef(fitted_model$fit_1)[n_k*2 + 1]
#   }
#   c23 = c_unid[1]
#   c13_2 = c_unid[2]
#   c24_3 = c_unid[3]
#   c14_23 = c_unid[4]
#   #survival copula rotation
#   if(fitted_model$copula_family %in% c("clayton", "gumbel")){
#     rotation = 180
#   }
#   else{
#     rotation = 0
#   }
#
#   pair_copulas = list(
#     list(
#       rvinecopulib::bicop_dist(
#         family = fitted_model$copula_family,
#         rotation = rotation,
#         parameters = c12
#       ),
#       rvinecopulib::bicop_dist(
#         family = copula_family2,
#         rotation = r_unid[1],
#         parameters = c23
#       ),
#       rvinecopulib::bicop_dist(
#         family = fitted_model$copula_family,
#         rotation = rotation,
#         parameters = c34
#       )
#     ),
#     list(
#       rvinecopulib::bicop_dist(
#         family = copula_family2,
#         rotation = r_unid[2],
#         parameters = c13_2
#       ),
#       rvinecopulib::bicop_dist(
#         family = copula_family2,
#         rotation = r_unid[3],
#         parameters = c24_3
#       )
#     ),
#     list(
#       rvinecopulib::bicop_dist(
#         family = copula_family2,
#         rotation = r_unid[4],
#         parameters = c14_23
#       )
#     )
#   )
#   copula_structure = rvinecopulib::dvine_structure(order = 1:4)
#   vine_cop = rvinecopulib::vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)
#
#   u_vec = rvinecopulib::rvinecop(n = n_prec, vinecop = vine_cop, cores = 1)
#   s0 = flexsurv::qsurvspline(p = u_vec[, 2],
#                              gamma = coef(fitted_model$fit_0)[1:n_k],
#                              knots = fitted_model$knots0)
#   t0 = flexsurv::qsurvspline(p = u_vec[, 1],
#                              gamma = coef(fitted_model$fit_0)[(n_k + 1):(2 * n_k)],
#                              knots = fitted_model$knott0)
#   s1 = flexsurv::qsurvspline(p = u_vec[, 3],
#                              gamma = coef(fitted_model$fit_1)[1:n_k],
#                              knots = fitted_model$knots1)
#   t1 = flexsurv::qsurvspline(p = u_vec[, 4],
#                              gamma = coef(fitted_model$fit_1)[(n_k + 1):(2 * n_k)],
#                              knots = fitted_model$knott1)
#   if (restr) {
#     pfs0 = pmin(s0, t0)
#     pfs1 = pmin(s1, t1)
#   }
#   else{
#     pfs0 = s0
#     pfs1 = s1
#   }
#   deltaS = pfs1 - pfs0
#   deltaT = t1 - t0
#   if (get_marg_tau) {
#     return(data.frame(
#       deltaS = deltaS,
#       deltaT = deltaT,
#       s0 = pfs0,
#       s1 = pfs1,
#       t0 = t0,
#       t1 = t1
#     ))
#   }
#   else{
#     return(data.frame(deltaS = deltaS, deltaT = deltaT))
#   }
#
# }

#
# unid_cop_sample = function(copula_family2, n_sim, cond_ind){
#   #sample uniformly parameters from spearman's correlation
#   if(copula_family2 == "frank"){
#     u = runif(n = 4*n_sim, min = -1, max = 1)
#     if (cond_ind) {
#       u[(1:(4*n_sim) %% 4) == 2] = 0
#       u[(1:(4*n_sim) %% 4) == 3] = 0
#     }
#     c = sapply(X = u, FUN = copula::iRho, copula = copula::frankCopula())
#     #functions for frank copula cannot handle parameter larger than abs(35)
#     c = ifelse(abs(c) > 35, sign(c)*35, c)
#   }
#   else if(copula_family2 == "gaussian"){
#     u = runif(n = 4*n_sim, min = -1, max = 1)
#     if (cond_ind) {
#       u[(1:(4*n_sim) %% 4) == 2] = 0
#       u[(1:(4*n_sim) %% 4) == 3] = 0
#     }
#     c = sapply(X = u, FUN = copula::iRho, copula = copula::ellipCopula(family = "normal"))
#   }
#   else if(copula_family2 == "clayton"){
#     u = runif(n = 4*n_sim, min = 0, max = 1)
#     if (cond_ind) {
#       u[(1:(4*n_sim) %% 4) == 2] = 0
#       u[(1:(4*n_sim) %% 4) == 3] = 0
#     }
#     c = sapply(X = u, FUN = copula::iRho, copula = copula::claytonCopula())
#     #functions for clayton copula cannot handle parameter larger than 28
#     c = ifelse(c > 28, 28, c)
#   }
#   else if(copula_family2 == "gumbel"){
#     u = runif(n = 4*n_sim, min = 0, max = 1)
#     if (cond_ind) {
#       u[(1:(4*n_sim) %% 4) == 2] = 0
#       u[(1:(4*n_sim) %% 4) == 3] = 0
#     }
#     c = sapply(X = u, FUN = copula::iRho, copula = copula::gumbelCopula())
#     #functions for gumbel copula cannot handle parameter values larger than 50
#     c = ifelse(c > 50, 50, c)
#   }
#   return(c)
# }
