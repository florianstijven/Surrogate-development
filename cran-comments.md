── R CMD check results ───────────────────────────────────────────────────────────────────────────── Surrogate 3.4.2 ────
Duration: 8m 25.5s

❯ checking package dependencies ... NOTE
  Imports includes 21 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking examples ... [108s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                                        user system elapsed
  sensitivity_intervals_Dvine          15.79   1.46   18.26
  sensitivity_analysis_SurvSurv_copula 15.98   0.95   18.47
  ICA_contcont_long_galecki             5.28   0.07    5.65
  ICA_alpha_ContCont                    5.19   0.08    5.61

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

* The timestamp note is related to the local machine in which R CMD CHECK was run and should
  not be produced on other machines. 
* The note package dependencies is a consequence of this package bundling 
  methods from many fields in statistics. 
* Only a few examples take longer than 5s to run. These examples are skipped on CRAN.
* In a previous version of the package, CRAN raised issues with URLs in the documentation.
  These issues were false positives and can be ignored. 
