## Resubmission
This is a resubmission. In this version I have:

* Skipped an additional test that caused issues on the OPENBLAS check of CRAN. 
  This test failed because of small numerical differences in the OPENBLAS check.
  On all other platforms and R versions, this test returns exactly the expected
  values.
