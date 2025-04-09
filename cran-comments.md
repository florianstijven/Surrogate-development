── R CMD check results ──────────────────────────────────────────────────── Surrogate 3.3.3.9000 ────
Duration: 7m 45.5s

❯ checking CRAN incoming feasibility ... [31s] NOTE
  Maintainer: 'Wim Van Der Elst <wim.vanderelst@gmail.com>'
  
  Version contains large components (3.3.3.9000)

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking examples ... [82s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                                        user system elapsed
  sensitivity_analysis_SurvSurv_copula 15.31   1.27   16.90
  sensitivity_intervals_Dvine          15.09   0.95   16.25

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

* The note is related to the local machine in which R CMD CHECK was run and should
not be produced on other machines. 
