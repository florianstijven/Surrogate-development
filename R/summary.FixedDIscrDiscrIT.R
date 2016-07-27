summary.FixedDiscrDiscrIT<-function(object,...,Object){
  
  if (missing(Object)){Object<-object}
  cat("\nFunction call:\n\n")
  print(Object$Call)
  
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials:",nrow(Object$Trial.Spec.Results))
  cat("\nTotal number of patients:",sum(Object$Trial.Spec.Results$Obs.per.trial))
  
  cat("\n\n\n## Summary of data retained in analysis")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials/patients retained in analysis for R2h:")
  cat("\nTrials:",sum(Object$Trial.Spec.Results$R2h.Included=="Y"))
  cat("\nPatients:",sum(Object$Trial.Spec.Results$Obs.per.trial[Object$Trial.Spec.Results$R2h.Included=="Y"]))
  
  cat("\n\nTotal number of trials/patients retained in analysis for R2ht:")
  cat("\nTrials:",sum(Object$Trial.Spec.Results$R2ht.Included=="Y"))
  cat("\nPatients:",sum(Object$Trial.Spec.Results$Obs.per.trial[Object$Trial.Spec.Results$R2ht.Included=="Y"]))
  
  cat("\n\n\n## Summary of separation in R2ht")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
  cat("\n\nNumber of trials where the penalized likelihood approach was applied to deal with separation in the:")
      
  cat("\nSurrogate:",sum(Object$Trial.Spec.Results$Separation.S=="YES"&Object$Trial.Spec.Results$R2ht.Included=="Y"))
  cat("\nTrue outcome:",sum(Object$Trial.Spec.Results$Separation.T=="YES"&Object$Trial.Spec.Results$R2ht.Included=="Y"))

  cat("\n\n\n# Information-theoretic surrogacy estimates summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  cat("Trial-level surrogacy (R2ht): \n")
  print(format(round(Object$R2ht, 4), nsmall = 4))
  cat("\nIndividual-level surrogacy (R2h): \n")
  print(format(round(Object$R2h, 4), nsmall = 4))
  
  }