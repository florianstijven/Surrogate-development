% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensitivity_analysis_copula.R
\name{sensitivity_analysis_copula}
\alias{sensitivity_analysis_copula}
\title{Perform Sensitivity Analysis for the Individual Causal Association based on a
D-vine copula model}
\usage{
sensitivity_analysis_copula(
  fitted_model,
  n_sim,
  eq_cond_association = TRUE,
  lower = c(-1, -1, -1, -1),
  upper = c(1, 1, 1, 1),
  degrees = c(0, 90, 180, 270),
  marg_association = TRUE,
  copula_family2 = fitted_model$copula_family[1],
  n_prec = 10000,
  ncores = 1,
  ICA_estimator = NULL
)
}
\arguments{
\item{fitted_model}{Returned value from \code{\link[=fit_copula_OrdOrd]{fit_copula_OrdOrd()}},
\code{\link[=fit_copula_OrdCont]{fit_copula_OrdCont()}}, or \code{\link[=fit_copula_ContCont]{fit_copula_ContCont()}}. This object
contains the estimated identifiable part of the joint distribution for the
potential outcomes.}

\item{n_sim}{Number of replications in the \emph{sensitivity analysis}. This value
should be large enough to sufficiently explore all possible values of the
ICA. The minimally sufficient number depends to a large extent on which
inequality assumptions are subsequently imposed (see Additional
Assumptions).}

\item{eq_cond_association}{Boolean.
\itemize{
\item \code{TRUE} (default): Assume that the association in \eqn{(\tilde{S}_1, T_0)' | \tilde{S}_0}
and \eqn{(\tilde{S}_0, T_1)' | \tilde{S}_1} are the same.
\item \code{FALSE}: There is not specific a priori relationship between the
above two associations.
}}

\item{lower}{(numeric) Vector of length 4 that provides the lower limit,
\eqn{\boldsymbol{a} = (a_{23}, a_{13;2}, a_{24;3},
  a_{14;23})'}. Defaults to \code{c(-1, -1, -1, -1)}. If the provided lower limit
is smaller than what is allowed for a particular copula family, then the
copula family's lowest possible value is used instead.}

\item{upper}{(numeric) Vector of length 4 that provides the upper limit,
\eqn{\boldsymbol{b} = (b_{23}, b_{13;2}, b_{24;3},
  b_{14;23})'}. Defaults to \code{c(1, 1, 1, 1)}.}

\item{degrees}{(numeric) vector with copula rotation degrees. Defaults to
\code{c(0, 90, 180, 270)}. This argument is not used for the Gaussian and Frank
copulas since they already allow for positive and negative associations.}

\item{marg_association}{Boolean.
\itemize{
\item \code{TRUE}: Return marginal association measures
in each replication in terms of Spearman's rho. The proportion of harmed,
protected, never diseased, and always diseased is also returned. See also
Value.
\item \code{FALSE} (default): No additional measures are returned.
}}

\item{copula_family2}{Copula family of the other bivariate copulas. For the
possible options, see \code{loglik_copula_scale()}. The elements of
\code{copula_family2} correspond to \eqn{(c_{23}, c_{13;2}, c_{24;3}, c_{14;23})}.}

\item{n_prec}{Number of Monte-Carlo samples for the \emph{numerical approximation}
of the ICA in each replication of the sensitivity analysis.}

\item{ncores}{Number of cores used in the sensitivity analysis. The
computations are computationally heavy, and this option can speed things up
considerably.}

\item{ICA_estimator}{Function that estimates the ICA between the first two
arguments which are numeric vectors. See also \code{\link[=compute_ICA_OrdOrd]{compute_ICA_OrdOrd()}},
\code{\link[=compute_ICA_OrdCont]{compute_ICA_OrdCont()}}, and \code{\link[=compute_ICA_ContCont]{compute_ICA_ContCont()}}.}
}
\value{
A data frame is returned. Each row represents one replication in the
sensitivity analysis. The returned data frame always contains the following
columns:
\itemize{
\item \code{R2H}, \code{sp_rho}: ICA as quantified by \eqn{R^2_H} and Spearman's rho, respectively.
\item \code{c12}, \code{c34}: estimated copula parameters.
\item \code{c23}, \code{c13_2}, \code{c24_3}, \code{c14_23}: sampled copula parameters of the
unidentifiable copulas in the D-vine copula. The parameters correspond to
the parameterization of the \code{copula_family2} copula as in the \code{copula}
R-package.
\item \code{r12}, \code{r34}: Fixed rotation parameters for the two identifiable copulas.
\item \code{r23}, \code{r13_2}, \code{r24_3}, \code{r14_23}: Sampled rotation parameters of the
unidentifiable copulas in the D-vine copula. These values are constant for
the Gaussian copula family since that copula is invariant to rotations.
}

The returned data frame also contains the following columns when
\code{marg_association} is \code{TRUE}:
\itemize{
\item \code{sp_s0s1}, \code{sp_s0t0}, \code{sp_s0t1}, \code{sp_s1t0}, \code{sp_s1t1}, \code{sp_t0t1}:
Spearman's rho between the corresponding potential outcomes. Note that these
associations refer to the observable potential outcomes. In contrast, the
estimated association parameters from \code{\link[=fit_copula_OrdOrd]{fit_copula_OrdOrd()}} and
\link{fit_copula_OrdCont} refer to associations on a latent scale.
}
}
\description{
Perform Sensitivity Analysis for the Individual Causal Association based on a
D-vine copula model
}
\section{Information-Theoretic Causal Inference Framework}{
The information-theoretic causal inference (ITCI) is a general framework to
evaluate surrogate endpoints in the single-trial setting (Alonso et al.,
2015). In this framework, we focus on the individual causal effects,
\eqn{\Delta S = S_1 - S_0} and \eqn{\Delta T = T_1 - T_0} where \eqn{S_z}
and \eqn{T_z} are the potential surrogate end true endpoint under treatment
\eqn{Z = z}.

In the ITCI framework, we say that \eqn{S} is a good surrogate for \eqn{T}
if
\emph{\eqn{\Delta S} conveys a substantial amount of information on \eqn{\Delta T}}
(Alonso, 2018). This amount of shared information can generally be quantified
by the mutual information between \eqn{\Delta S} and \eqn{\Delta T},
denoted by \eqn{I(\Delta S; \Delta T)}. However, the mutual information lies
in \eqn{[0, + \infty]} which complicates the interpretation. In addition,
the mutual information may not be defined in specific scenarios where
absolute continuity of certain probability measures fails. Therefore, the
mutual information is transformed, and possibly modified, to enable a simple
interpretation in light of the definition of surrogacy. The resulting measure
is termed the individual causal association (ICA). This is explained in
the next sections.

While the definition of surrogacy in the ITCI framework rests on information
theory, shared information is closely related to statistical association. Hence,
we can also define the ICA in terms of statistical association measures, like
Spearman's rho and Kendall's tau. The advantage of the latter are that they
are well-known, simple and rank-based measures of association.
}

\section{Individual Causal Association}{
Many association measures can operationalize the ICA. For each setting, we
consider one default definition for the ICA which follows from the mutual
information.
\subsection{Continuous-Continuous}{

The ICA is defined as the squared informational coefficient of correlation
(SICC or \eqn{R^2_H}), which is a transformation of the mutual information
to the unit interval: \deqn{R^2_h = 1 - e^{-2 \cdot I(\Delta S; \Delta T)}}
where 0 indicates independence, and 1 a functional relationship between
\eqn{\Delta S} and \eqn{\Delta T}. If \eqn{(\Delta S, \Delta T)'} is bivariate
normal, the ICA equals the Pearson correlation between \eqn{\Delta S} and
\eqn{\Delta T}.
}

\subsection{Ordinal-Continuous}{

The ICA is defined as the following transformation of the mutual information:
\deqn{R^2_H = \frac{I(\Delta S; \Delta T)}{H(\Delta T)},}
where \eqn{I(\Delta S; \Delta T)} is the mutual information and \eqn{H(\Delta T)}
the entropy.
}

\subsection{Ordinal-Ordinal}{

The ICA is defined as the following transformation of the mutual information:
\deqn{R^2_H = \frac{I(\Delta S; \Delta T)}{\min \{H(\Delta S), H(\Delta T) \}},}
where \eqn{I(\Delta S; \Delta T)} is the mutual information, and \eqn{H(\Delta S)}
and \eqn{H(\Delta T)} the entropy of \eqn{\Delta S} and \eqn{\Delta T},
respectively.
}
}

\section{Sensitivity Analysis}{
\subsection{Monte Carlo Approach}{

Because \eqn{S_0} and \eqn{S_1} are never simultaneously observed in the same
patient, \eqn{\Delta S} is not observable, and analogously for \eqn{\Delta
T}. Consequently, the ICA is unidentifiable. This is solved by considering a
(partly identifiable) model for the full vector of potential outcomes,
\eqn{(T_0, S_0, S_1, T_1)'}. The identifiable parameters are estimated. The
unidentifiable parameters are sampled from their parameters space in each
replication of a sensitivity analysis. If the number of replications
(\code{n_sim}) is sufficiently large, the entire parameter space for the
unidentifiable parameters will be explored/sampled. In each replication, all
model parameters are "known" (either estimated or sampled). Consequently, the
ICA can be computed in each replication of the sensitivity analysis.

The sensitivity analysis thus results in a set of values for the ICA. This
set can be interpreted as \emph{all values for the ICA that are compatible with
the observed data}. However, the range of this set is often quite broad; this
means there remains too much uncertainty to make judgements regarding the
worth of the surrogate. To address this unwieldy uncertainty, additional
assumptions can be used that restrict the parameter space of the
unidentifiable parameters. This in turn reduces the uncertainty regarding the
ICA.
}

\subsection{Intervals of Ignorance and Uncertainty}{

The results of the sensitivity analysis can be formalized (and summarized) in
intervals of ignorance and uncertainty using \code{\link[=sensitivity_intervals_Dvine]{sensitivity_intervals_Dvine()}}.
}
}

\section{Additional Assumptions}{
There are two possible types of assumptions that restrict the parameter space
of the unidentifiable parameters: (i) \emph{equality} type of assumptions, and
(ii) \emph{inequality} type of assumptions. These are discussed in turn in the
next two paragraphs.

The equality assumptions have to be incorporated into the sensitivity
analysis itself. Only one type of equality assumption has been implemented;
this is the conditional independence assumption:
\deqn{\tilde{S}_0 \perp T_1 | \tilde{S}_1 \; \text{and} \;
\tilde{S}_1 \perp T_0 | \tilde{S}_0 .} This can informally be
interpreted as ``what the control treatment does to the surrogate does not
provide information on the true endpoint under experimental treatment if we
already know what the experimental treatment does to the surrogate", and
analogously when control and experimental treatment are interchanged. Note
that \eqn{\tilde{S}_z} refers to either the actual potential surrogate
outcome, or a latent version. This depends on the content of \code{fitted_model}.

The inequality type of assumptions have to be imposed on the data frame that
is returned by the current function; those assumptions are thus imposed
\emph{after} running the sensitivity analysis. If \code{marginal_association} is set to
\code{TRUE}, the returned data frame contains additional unverifiable quantities
that differ across replications of the sensitivity analysis: (i) the
unconditional Spearman's \eqn{\rho} for all pairs of (observable/non-latent)
potential outcomes, and (ii) the proportions of the population strata as
defined by Nevo and Gorfine (2022) if semi-competing risks are present. More
details on the interpretation and use of these assumptions can be found in
Stijven et al. (2024).
}

\references{
Alonso, A. (2018). An information-theoretic approach for the evaluation of
surrogate endpoints. In Wiley StatsRef: Statistics Reference Online. John
Wiley & Sons, Ltd.

Alonso, A., Van der Elst, W., Molenberghs, G., Buyse, M., and Burzykowski, T.
(2015). On the relationship between the causal-inference and meta-analytic
paradigms for the validation of surrogate endpoints. Biometrics 71, 15–24.
}
