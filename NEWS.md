# Surrogate (development version)

* In the sensitivity analyses, we now have to option to set conditional 
  association parameters equal to each other through `eq_cond_association` 
* Bug fixed in `sensitivity_analysis_SurvSurv()`. When the unidentifiable copula
  family was `"clayton"`, then, if the copula parameter on Spearman's rho scale
  was very close to zero, the conversion to the original copula parameter scale
  was wrong.

# Surrogate 3.2

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
