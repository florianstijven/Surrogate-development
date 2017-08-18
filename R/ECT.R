ECT <- function(Perc=.95, H_Max, N){
Delta_H <- qchisq(c(1-Perc), df=8, lower.tail=FALSE) / (2 * N)

Lower_H <- H_Max - Delta_H
Upper_H <- H_Max

fit <-   
  list(Lower_H=as.numeric(Lower_H), Upper_H=as.numeric(Upper_H), Perc=as.numeric(Perc), Call=match.call())  

class(fit) <- "ECT"
fit

}


summary.ECT <- function(object, ..., Object){

  if (missing(Object)){Object <- object} 

  cat("Results application Entropy Concentration Theorem:\n")
  cat("--------------------------------------------------\n")
  cat((Object$Perc)*100, "% of all vectors p satisfy ", 
      round(Object$Lower_H, digits = 4), " <= H(p) <= ", round(Object$Upper_H, digits=4), sep="")
}


