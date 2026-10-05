#' Assess surrogacy in the information-theoretic causal-inference framework (Individual Causal Association, ICA)
#' for a continuous surrogate and true endpoint measured repeatedly over time in a single-trial setting
#'
#'
#' This function quantifies surrogacy under Galecki's model. It is a special
#' case without subject-specific random effects. See Details below.
#'
#' @param V0,V1 Treatment-specific residual 2 x 2 covariance matrices for T,S.
#' @param p (numeric) number of time points.
#' @param grid_V Sensitivity grid for unidentified residual correlations.
#' @param tol Numerical tolerance used when checking positive definiteness and validating ICA values.
#' @param return_candidates Logical. If \code{TRUE}, the retained admissible residual covariance
#'  completion candidates are returned.
#'
#' @return
#' * Total.Num.Matrices: An object of class numeric that contains the total number of matrices that can be formed as based on the user-specified correlations in the function call.
#' * Num.Pos.Def.V: An object of class numeric that contains the number of
#'   positive definite `V` matrices.
#' * R2_Lambda: A scalar or vector that contains the individual causal association \eqn{R_{\Lambda}^2}.
#' * V_candidates: retained residual covariance candidates, if requested.
#'
#' @details
#' # Galecki's model:
#' The evolution of \eqn{y_j} over time is defined using a marginal model without subject-specific random effects.
#' \deqn{ T_{0j} = x_{T0j} \gamma_{T0} + \varepsilon_{T0j}}
#' \deqn{T_{1j} = x_{T1j} \gamma_{T1} + \varepsilon_{T1j}}
#' \deqn{S_{0j} = x_{S0j} \gamma_{S0} + \varepsilon_{S0j}}
#' \deqn{S_{1j} = x_{S1j} \gamma_{S1} + \varepsilon_{S1j}}
#'
#' It can be written as \eqn{y=x \gamma + \varepsilon}. It is assumed that \eqn{\varepsilon \sim N(0, \Sigma)},
#' or, equivalently, \eqn{y \sim N(x \gamma, \Sigma_y)}. Essentially,
#' Galecki (1994) assumed that the covariance matrix is
#' separable, i.e. it can be written as \eqn{\Sigma_y = \Sigma = R \otimes V} , where \eqn{V} is assumed to be unstructured and
#' invariant across time, and \eqn{R} reflects a general correlation matrix.
#'
#'
#' Under Galecki's model,
#' the vector of individual causal treatment effects \eqn{\Delta} is:
#' \deqn{ \Delta = \mu_{\Delta}+ \varepsilon_{\Delta} }
#' and the marginal covariance matrix of the individual causal treatment effects is
#' \eqn{\Sigma_{\Delta} = \Sigma_u \otimes R}.
#'
#' Under Galecki's model, the ICA is given by:
#' \deqn{ R_{\Lambda}^2 = 1 - (1-\rho_{u}^2)^p }
#'
#' In this model \eqn{\rho_u=\rho_{\Delta}=corr(\Delta T, \Delta S)}, which is the metric used in the cross-sectional setting.
#'
#' \deqn{
#' \rho_u =
#'  \frac{
#'    \sqrt{\sigma_{T0T0} \, \sigma_{S0S0}} \, \rho_{T0S0} +
#'      \sqrt{\sigma_{T1T1} \, \sigma_{S1S1}} \, \rho_{T1S1} -
#'      \sqrt{\sigma_{T1T1} \, \sigma_{S0S0}} \, \rho_{T1S0} -
#'      \sqrt{\sigma_{T0T0} \, \sigma_{S1S1}} \, \rho_{T0S1}
#'  }{
#'    \sqrt{
#'      (\sigma_{T0T0} + \sigma_{T1T1} - 2 \sqrt{\sigma_{T0T0} \, \sigma_{T1T1}} \, \rho_{T0T1})
#'      (\sigma_{S0S0} + \sigma_{S1S1} - 2 \sqrt{\sigma_{S0S0} \, \sigma_{S1S1}} \, \rho_{S0S1})
#'    }
#'  }
#' }
#'
#' @references
#' Deliorman, G., Pardo, M. C., Van der Elst W., Molenberghs G., & Alonso, A. (submitted in 2026).
#' Assessing Surrogate Endpoints in Longitudinal Studies: An Information-Theoretic
#' Causal Inference Approach.
#'
#'
#' @examples
#' V0 <- matrix(c(0.04513345, 0.04375355,
#'                0.04375355, 0.04866820),
#'               nrow = 2, byrow = TRUE,
#'               dimnames = list(c("T", "S"), c("T", "S")))
#'
#'
#' V1 <- matrix(c(0.05231045, 0.05177398,
#'                0.05177398, 0.05702684),
#'                nrow = 2, byrow = TRUE,
#'                dimnames = list(c("T", "S"), c("T", "S")))
#'
#' ica_galecki<- ICA_contcont_long_galecki(V0 = V0,  V1 = V1,  p=5,
#'                           grid_V = seq(-1, 1, by = 0.5),
#'                           return_candidates = FALSE,
#'                           tol = 1e-12)
#'
#' summary_longitudinal_ica(ica_galecki)
#'
#'
#'
#' @export
ICA_contcont_long_galecki <- function(V0,
                                      V1,
                                      p,
                                      grid_V = seq(-1, 1, by = 0.1),
                                      return_candidates = FALSE,
                                      tol = 1e-12) {
  p <- as.integer(p)
  if (length(p) != 1L || is.na(p) || p < 1L) stop("`p` must be a positive integer.")

  cc <- .surr_completion_components(V0, V1, grid_V, tol = tol,
                                    return_candidates = return_candidates)
  comp <- cc$components
  R2 <- 1 - (1 - comp$rho2)^p
  R2 <- .surr_validate_R2(R2)

  out <- list(
    Total.Num.Matrices = cc$total,
    Num.Pos.Def.V = nrow(comp),
    R2_Lambda = R2,
    V_candidates = if (return_candidates) comp else NULL,
    Call = match.call()
  )
  class(out) <- c("ICA_contcont_long_galecki", "surrogate_longitudinal_ica")
  out
}


