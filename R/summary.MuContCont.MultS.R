
summary.MufixedContCont.MultS <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1])
  cat("\nTotal number of patients in experimental treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat==1]), 
      "\nTotal number of patients in control treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat!=1])) 
  
  cat("\n\n\n# Results (CI's based on Lee's approach; Lee, 1971)")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  print(format(round(Object$Trial.R2.Lee, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Indiv.R2.Lee, 4), nsmall = 4))
  cat("\n")
  
  if (is.na(Object$Trial.R2.Boot$`CI lower limit`)==FALSE){
    cat("\n\n")  
    cat("\n\n\n# CI's based on a non-parametric bootstrap")
    cat("\n\n")
    print(format(round(Object$Trial.R2.Boot, 4), nsmall = 4))
    cat("\n")
    print(format(round(Object$Indiv.R2.Boot, 4), nsmall = 4))
    cat("\n")
    
  }
  
}


summary.MumixedContCont.MultS <- function(object, ..., Object){
  
  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Data summary and descriptives")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\nTotal number of trials: ", nrow(Object$Obs.Per.Trial))
  cat("\nTotal number of patients: ", dim(Object$Data.Analyze)[1])
  cat("\nTotal number of patients in experimental treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat==1]), 
      "\nTotal number of patients in control treatment group: ", length(Object$Data.Analyze$Treat[Object$Data.Analyze$Treat!=1])) 
  
  cat("\n\n\n# Results (CI's based on Lee's approach; Lee, 1971)")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cat("\n\n")
  print(format(round(Object$Trial.R2.Lee, 4), nsmall = 4))
  cat("\n")
  print(format(round(Object$Indiv.R2.Lee, 4), nsmall = 4))
  cat("\n")
  

}