# Surrogate 3.3

## Surrogate 3.3.0

* Added and updated functions for the survival-categorical, survival-binary, and
survival-continuous settings in the meta-analytic framework.
* Added functions to perform the two-stage federated surrogacy evaluation.

## Surrogate 3.3.1

* General implementation of the information-theoretic causal inference framework
based on D-vine copula models for the ordinal-ordinal, ordinal-continuous, and
continuous-continuous settings.

## Surrogate 3.3.3

* Fix tests that fail on CRAN. 

# Surrogate 3.2

## Surrogate 3.2.6

* Updated functions for the surrogate predictive function in the binary-continuous
setting (information-theoretic causal inference framework).

## Surrogate 3.2.5

* Added functions `survbin()` and `survcat()` to evaluate categorical surrogates 
for time-to-event true endpoints in the meta-analytic framework.

## Surrogate 3.2.4

* In the sensitivity analyses, we now have to option to set conditional 
  association parameters equal to each other through `eq_cond_association` 
* Bug fixed in `sensitivity_analysis_SurvSurv()`. When the unidentifiable copula
  family was `"clayton"`, then, if the copula parameter on Spearman's rho scale
  was very close to zero, the conversion to the original copula parameter scale
  was wrong.

## Surrogate 3.2.2

* Added uncertainty intervals to summarize sensitivity analyses
* Extended flexibility for survival-survival vine copula models and associated sensitivity analysis
* Add functions to evaluate multivariate surrogates in the meta-analytic framework

## Surrogate 3.2.1

* Added `ICA.BinCont.BS()` function. This function allows for evaluating a 
continuous for a binary true endpoint in the information-theoretic causal 
inference framework. This function serves the same purpose as `ICA.BinCont()`,
but additionally allows for taking sampling variability into account.

## Surrogate 3.2.0

* Added a `NEWS.md` file to track changes to the package.
