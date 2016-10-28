summary.MaxEntContCont <- function(object, ..., Object){
  
  options(digits = 4)
  
  if (missing(Object)){Object <- object} 
  
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n# Maximum entropy ICA:")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~\n")  
  cat(Object$ICA.Max.Ent) #, " (with entropy = ", Object$Max.Ent, ")\n", sep = "")
  
  cat("\n\n# Obtained under assumed correlations:")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")  
  print(Object$Table.ICA.Entropy[1,3:8], row.names = "")
  cat("\n")
}
