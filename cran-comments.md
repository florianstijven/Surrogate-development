── R CMD check results ─────────────────────────────────────────────── Surrogate 3.4.1.9000 ────
Duration: 8m 34.4s

❯ checking CRAN incoming feasibility ... [21s] NOTE
  Maintainer: 'Wim Van Der Elst <wim.vanderelst@gmail.com>'
  
  Version contains large components (3.4.1.9000)

❯ checking package dependencies ... NOTE
  Imports includes 21 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking examples ... [109s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                                        user system elapsed
  sensitivity_intervals_Dvine          15.86   1.45   17.73
  sensitivity_analysis_SurvSurv_copula 15.20   0.86   18.83
  ICA_contcont_long_galecki             5.42   0.08    5.79
  ICA_alpha_ContCont                    5.19   0.06    5.83
  marginal_gof_scr_S_plot               3.52   0.25    5.09

0 errors ✔ | 0 warnings ✔ | 4 notes ✖

* The timestamp note is related to the local machine in which R CMD CHECK was run and should
not be produced on other machines. 
* The note package dependencies is a consequence of this package bundling 
  methods from many fields in statistics. 
* Only a few examples take longer than 5s to run. On more powerful systems, this
  note may not be produced.
