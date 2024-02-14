#' Sensitivity analysis for individual causal association
#'
#' The [sensitivity_analysis_SurvSurv_copula()] function performs the
#' sensitivity analysis for the individual causal association (ICA) as described
#' by Stijven et al. (2024).
#'
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
#' # Surrogacy in The Survival-Survival Setting
#'
#' ## General Introduction
#'
#' Stijven et al. (2024) proposed to quantify the ICA through the squared
#' informational coefficient of correlation (SICC or \eqn{R^2_H}), which is a
#' transformation of the mutaul information to the unit interval: \deqn{R^2_H =
#' 1 - e^{-2 \cdot I(\Delta S; \Delta T)}} where 0 indicates independence, and 1
#' a functional relationship between \eqn{\Delta S} and \eqn{\Delta T}. The ICA
#' (or a modified version, see next) is returned by
#' [sensitivity_analysis_SurvSurv_copula()]. Concurrently, the Spearman's
#' correlation between \eqn{\Delta S} and \eqn{\Delta T} is also returned.
#'
#' ## Issues with Composite Endpoints
#'
#' In the survival-survival setting where the surrogate is a composite endpoint,
#' care should be taken when defining the mutual information. Indeed, when
#' \eqn{S_z} is progression-free survival and \eqn{T_z} is overall survival,
#' there is a probability atom in the joint distribution of \eqn{(S_z, T_z)'}
#' because \eqn{P(S_z = T_z) > 0}. In other words, there are patient that die
#' before progressing. While this probability atom is correctly taken into
#' account in the models fitted by [fit_model_SurvSurv()], this probability atom
#' reappears when considering the distribution of \eqn{(\Delta S, \Delta T)'}
#' because \eqn{P(\Delta S = \Delta T) > 0} if we are considering PFS and OS.
#'
#' Because of the atom in the distribution of \eqn{(\Delta S, \Delta T)'}, the
#' corresponding mutual information is not defined. To solve this, the mutual
#' information is computed excluding the patients for which \eqn{\Delta S =
#' \Delta T} when `composite = TRUE`. The proportion of excluded patients is, among
#' other things, returned when `marginal_association = TRUE`. This is the proportion
#' of "never" patients following the classification of Nevo and Gorfine (2022).
#' See also Additional Assumptions.
#'
#' This modified version of the ICA quantifies the surrogacy of \eqn{S} when
#' "adjusted for the composite nature of \eqn{S}". Indeed, we exclude patients
#' where \eqn{\Delta S} perfectly predicts \eqn{\Delta T} *just because \eqn{S}
#' is a composite of \eqn{T} (and other variables).
#'
#' Other (rank-based) statistical measures of association, however,
#' remain well-defined and are thus computed without excluding any patients.
#'
#' # Sensitivity Analysis
#'
#' ## Monte Carlo Approach
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
#' ## Intervals of Ignorance and Uncertainty
#'
#' The results of the sensitivity analysis can be formalized (and summarized) in
#' intervals of ignorance and uncertainty using [sensitivity_intervals_Dvine()].
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
#' \deqn{\tilde{S}_0 \perp T_1 | \tilde{S}_1 \; \text{and} \;
#' \tilde{S}_1 \perp T_0 | \tilde{S}_0 .} This can informally be
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
#' Stijven et al. (2024).
#'
#'
#' @param fitted_model Returned value from [fit_model_SurvSurv()]. This object
#'   contains the estimated identifiable part of the joint distribution for the
#'   potential outcomes.
#' @param degrees (numeric) vector with copula rotation degrees. Defaults to
#'   `c(0, 90, 180, 270)`. This argument is not used for the Gaussian and Frank
#'   copulas since they already allow for positive and negative associations.
#' @param n_sim Number of replications in the *sensitivity analysis*. This value
#'   should be large enough to sufficiently explore all possible values of the
#'   ICA. The minimally sufficient number depends to a large extent on which
#'   inequality assumptions are subsequently imposed (see Additional
#'   Assumptions).
#' @param n_prec Number of Monte-Carlo samples for the *numerical approximation*
#'   of the ICA in each replication of the sensitivity analysis.
#' @param ncores Number of cores used in the sensitivity analysis. The
#'   computations are computationally heavy, and this option can speed things up
#'   considerably.
#' @param marg_association Boolean.
#' * `TRUE`: Return marginal association measures
#'   in each replication in terms of Spearman's rho. The proportion of harmed,
#'   protected, never diseased, and always diseased is also returned. See also
#'   Value.
#' * `FALSE` (default): No additional measures are returned.
#' @param cond_ind Boolean.
#' * `TRUE`: Assume conditional independence (see Additional Assumptions).
#' * `FALSE` (default): Conditional independence is not assumed.
#' @param sample_plots Indices for replicates in the sensitivity analysis for
#'   which the sampled individual treatment effects are plotted. Defaults to
#'   `NULL`: no plots are displayed.
#' @inheritParams estimate_mutual_information_SurvSurv
#' @inheritParams sample_copula_parameters
#' @inheritParams compute_ICA_BinCont
#' @inheritParams compute_ICA_SurvSurv
#'
#' @return A data frame is returned. Each row represents one replication in the
#'   sensitivity analysis. The returned data frame always contains the following
#'   columns:
#' * `ICA`, `sp_rho`: ICA as quantified by \eqn{R^2_h(\Delta S^*, \Delta T^*)} and
#'   \eqn{\rho_s(\Delta S, \Delta T)}.
#' * `c23`, `c13_2`, `c24_3`, `c14_23`: sampled copula parameters of the
#'   unidentifiable copulas in the D-vine copula. The parameters correspond to
#'   the parameterization of the `copula_family2` copula as in the `copula`
#'   R-package.
#' * `r23`, `r13_2`, `r24_3`, `r14_23`: sampled rotation parameters of the
#'   unidentifiable copulas in the D-vine copula. These values are constant for
#'   the Gaussian copula family since that copula is invariant to rotations.
#'
#'   The returned data frame also contains the following columns when
#'   `get_marg_tau` is `TRUE`:
#' * `sp_s0s1`, `sp_s0t0`, `sp_s0t1`, `sp_s1t0`, `sp_s1t1`, `sp_t0t1`:
#'   Spearman's \eqn{\rho} between the corresponding potential outcomes. Note
#'   that these associations refer to the potential time-to-composite events
#'   and/or time-to-true endpoint event. In contrary, the estimated association
#'   parameters from [fit_model_SurvSurv()] refer to associations between the
#'   time-to-surrogate event and time-to true endpoint event. Also note that
#'   `sp_s1t1` is constant whereas `sp_s0t0` is not. This is a particularity of
#'   the MC procedure to calculate both measures and thus not a bug.
#' * `prop_harmed`, `prop_protected`, `prop_always`, `prop_never`: proportions
#'   of the corresponding population strata in each replication. These are
#'   defined in Nevo and Gorfine (2022).
#' @export
#'
#' @references
#' Alonso, A. (2018). An information-theoretic approach for the evaluation of
#' surrogate endpoints. In Wiley StatsRef: Statistics Reference Online. John
#' Wiley & Sons, Ltd.
#'
#' Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., and Burzykowski, T.
#' (2015). On the relationship between the causal-inference and meta-analytic
#' paradigms for the validation of surrogate endpoints. Biometrics 71, 15–24.
#'
#' Stijven, F., Alonso, a., Molenberghs, G., Van Der Elst, W., Van Keilegom, I.
#' (2024). An information-theoretic approach to the evaluation of time-to-event
#' surrogates for time-to-event true endpoints based on causal inference.
#'
#' Nevo, D., & Gorfine, M. (2022). Causal inference for semi-competing risks
#' data. Biostatistics, 23 (4), 1115-1132
#'
#' @examplesIf identical(Sys.getenv("NOT_CRAN"), "true")
#' # Load Ovarian data
#' data("Ovarian")
#' # Recode the Ovarian data in the semi-competing risks format.
#' data_scr = data.frame(
#'   ttp = Ovarian$Pfs,
#'   os = Ovarian$Surv,
#'   treat = Ovarian$Treat,
#'   ttp_ind = ifelse(
#'     Ovarian$Pfs == Ovarian$Surv &
#'       Ovarian$SurvInd == 1,
#'     0,
#'     Ovarian$PfsInd
#'   ),
#'   os_ind = Ovarian$SurvInd
#' )
#' # Fit copula model.
#' fitted_model = fit_model_SurvSurv(data = data_scr,
#'                                   copula_family = "clayton",
#'                                   n_knots = 1)
#' # Illustration with small number of replications and low precision
#' sens_results = sensitivity_analysis_SurvSurv_copula(fitted_model,
#'                   n_sim = 5,
#'                   n_prec = 2000,
#'                   copula_family2 = "clayton",
#'                   cond_ind = TRUE)
#' # Compute intervals of ignorance and uncertainty. Again, the number of
#' # bootstrap replications should be larger in practice.
#' sensitivity_intervals_Dvine(fitted_model, sens_results, B = 10)
sensitivity_analysis_SurvSurv_copula = function(fitted_model,
                                                composite = TRUE,
                                                n_sim,
                                                cond_ind,
                                                lower = c(-1, -1, -1, -1),
                                                upper = c(1, 1, 1, 1),
                                                degrees = c(0, 90, 180, 270),
                                                marg_association = TRUE,
                                                copula_family2 = fitted_model$copula_family[1],
                                                n_prec = 5e3,
                                                ncores = 1,
                                                sample_plots = NULL,
                                                mutinfo_estimator = NULL,
                                                restr_time = +Inf) {
  # If copula_family2 contains only 1 element, this vector is appended to
  # the correct length.
  copula_family1 = fitted_model$copula_family
  if(length(copula_family1) == 1) copula_family1 = rep(copula_family2, 2)
  if(length(copula_family2) == 1) copula_family2 = rep(copula_family2, 4)
  # Extract relevant estimated parameters/objects for the fitted copula model.

  ks0 = length(fitted_model$knots0)
  ks1 = length(fitted_model$knots1)
  kt0 = length(fitted_model$knott0)
  kt1 = length(fitted_model$knott1)

  gammas0 = force(coef(fitted_model$fit_0)[1:ks0])
  knots0 = force(fitted_model$knots0)
  gammas1 = force(coef(fitted_model$fit_1)[1:ks1])
  knots1 = force(fitted_model$knots1)

  gammat0 = force(coef(fitted_model$fit_0)[(ks0 + 1):(ks0 + kt0)])
  knott0 = force(fitted_model$knott0)
  gammat1 = force(coef(fitted_model$fit_1)[(ks1 + 1):(ks1 + kt1)])
  knott1 = force(fitted_model$knott1)

  copula_rotations = force(fitted_model$copula_rotations)

  q_S0 = function(p) {
    flexsurv::qsurvspline(p = p,
                          gamma = gammas0,
                          knots = knots0)
  }
  q_S1 = function(p) {
    flexsurv::qsurvspline(p = p,
                          gamma = gammas1,
                          knots = knots1)
  }
  q_T0 = function(p) {
    flexsurv::qsurvspline(p = p,
                          gamma = gammat0,
                          knots = knott0)
  }
  q_T1 = function(p) {
    flexsurv::qsurvspline(p = p,
                          gamma = gammat1,
                          knots = knott1)
  }

  # number of parameters
  n_par = length(coef(fitted_model$fit_0))

  # Pull association parameters from estimated parameter vectors. The
  # association parameter is always the last one in the corresponding vector
  c12 = coef(fitted_model$fit_0)[n_par]
  c34 = coef(fitted_model$fit_1)[n_par]
  # For the Gaussian copula, fisher's Z transformation was applied. We have to
  # backtransform to the correlation scale in that case.
  if (copula_family1[1] == "gaussian") {
    c12 = (exp(2 * c12) - 1) / (exp(2 * c12) + 1)
  }
  if (copula_family1[2] == "gaussian") {
    c34 = (exp(2 * c34) - 1) / (exp(2 * c34) + 1)
  }
  c = sample_copula_parameters(
    copula_family2 = copula_family2,
    n_sim = n_sim,
    cond_ind = cond_ind,
    lower = lower,
    upper = upper
  )
  c = cbind(rep(c12, n_sim), c[, 1], rep(c34, n_sim), c[, 2:4])
  c_list = purrr::map(.x = split(c, seq(nrow(c))), .f = as.double)
  # Sample rotation parameters of the unidentifiable copula's family does not
  # allow for negative associations.
  r = sample_rotation_parameters(n_sim, degrees = degrees)
  if (copula_family2[1] %in% c("gaussian", "frank")) r[, 1] = rep(0, n_sim)
  if (copula_family2[2] %in% c("gaussian", "frank")) r[, 2] = rep(0, n_sim)
  if (copula_family2[3] %in% c("gaussian", "frank")) r[, 3] = rep(0, n_sim)
  if (copula_family2[4] %in% c("gaussian", "frank")) r[, 4] = rep(0, n_sim)
  r = cbind(rep(copula_rotations[1], n_sim),
            r[, 1],
            rep(copula_rotations[2], n_sim),
            r[, 2:4])
  r_list = purrr::map(.x = split(r, seq(nrow(r))), .f = as.double)
  # For every set of sampled unidentifiable parameters, compute the
  # required quantities.
  #put all other arguments in a list for the apply function
  MoreArgs = list(
    copula_family1 = copula_family1,
    copula_family2 = copula_family2,
    n_prec = n_prec,
    q_S0 = q_S0,
    q_T0 = q_T0,
    q_S1 = q_S1,
    q_T1 = q_T1,
    composite = composite,
    marginal_sp_rho = marg_association,
    seed = 1,
    mutinfo_estimator = mutinfo_estimator,
    restr_time = restr_time
  )
  if (ncores > 1 & requireNamespace("parallel")) {
    if (!is.null(sample_plots)){
      warning("Sample plots cannot be displayed when sensitivity analysis is executed in parallel.")
    }
    cl  <- parallel::makeCluster(ncores)
    #helper function
    # surrogacy_sample_sens <- surrogacy_sample_sens
    print("Starting parallel simulations")
    # surrogacy_sample_sens

    # Get current search path and set the same search path in the new instances
    # of R. Usually, this would not be necessary, but if the user changed the
    # search path before running this function, there could be an issue.
    search_path = .libPaths()
    force(search_path)
    parallel::clusterExport(
      cl = cl,
      varlist = c("fitted_model", "search_path"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl = cl, expr = .libPaths(new = search_path))
    temp = parallel::clusterMap(
      cl = cl,
      fun = compute_ICA_SurvSurv,
      copula_par = c_list,
      rotation_par = r_list,
      MoreArgs = MoreArgs
    )
    print("Finishing parallel simulations")
    on.exit(expr = {
      parallel::stopCluster(cl)
      rm("cl")
    })
  }
  else if (ncores == 1) {
    sample_plots_indicator = rep(FALSE, n_sim)
    sample_plots_indicator[sample_plots] = TRUE
    temp = mapply(
      FUN = compute_ICA_SurvSurv,
      copula_par = c_list,
      rotation_par = r_list,
      plot_deltas = sample_plots_indicator,
      MoreArgs = MoreArgs,
      SIMPLIFY = FALSE
    )
  }

  measures_df = t(as.data.frame(temp))
  if (marg_association) {
    colnames(measures_df) = c(
      "ICA",
      "sp_rho",
      "sp_rho_t0s0",
      "sp_rho_t0s1",
      "sp_rho_t0t1",
      "sp_rho_s0s1",
      "sp_rho_s0t1",
      "sp_rho_s1t1",
      "prop_harmed",
      "prop_protected",
      "prop_always",
      "prop_never"
    )
  }
  else{
    colnames(measures_df) = c("ICA", "sp_rho")
  }
  rownames(measures_df) = NULL

  colnames(c) = c("c12", "c23", "c34",  "c13_2", "c24_3", "c14_23")
  colnames(r) = c("r12", "r23", "r34",  "r13_2", "r24_3", "r14_23")

  sens_results = cbind(as.data.frame(measures_df), c, r)
  attr(sens_results, which = "copula_family1") = copula_family1
  attr(sens_results, which = "copula_family2") = copula_family2
  attr(sens_results, which = "composite") = composite
  return(sens_results)

}
