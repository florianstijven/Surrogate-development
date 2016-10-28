MaxEntContCont <- function(x, T0T0, T1T1, S0S0, S1S1){
entropy <- NULL

for (i in 1: dim(x$Pos.Def)[1]){

T0T1 <- x$Pos.Def$T0T1[i]  
T0S0 <- x$Pos.Def$T0S0[i]
T0S1 <- x$Pos.Def$T0S1[i]
T1S0 <- x$Pos.Def$T1S0[i]
T1S1 <- x$Pos.Def$T1S1[i]
S0S1 <- x$Pos.Def$S0S1[i]
Sigma_c <- diag(4)         
Sigma_c[2,1] <- Sigma_c[1,2] <- T0T1 * (sqrt(T0T0)*sqrt(T1T1))
Sigma_c[3,1] <- Sigma_c[1,3] <- T0S0 * (sqrt(T0T0)*sqrt(S0S0))
Sigma_c[4,1] <- Sigma_c[1,4] <- T0S1 * (sqrt(T0T0)*sqrt(S1S1))
Sigma_c[3,2] <- Sigma_c[2,3] <- T1S0 * (sqrt(T1T1)*sqrt(S0S0))
Sigma_c[4,2] <- Sigma_c[2,4] <- T1S1 * (sqrt(T1T1)*sqrt(S1S1))
Sigma_c[4,3] <- Sigma_c[3,4] <- S0S1 * (sqrt(S0S0)*sqrt(S1S1))
Sigma_c[1,1] <- T0T0
Sigma_c[2,2] <- T1T1
Sigma_c[3,3] <- S0S0
Sigma_c[4,4] <- S1S1

entropy_here <- 2 * log2(2 * pi) + (.5 * log2(det(Sigma_c)))
entropy <- cbind(entropy, entropy_here)
}

results <- data.frame(cbind(x$ICA, as.numeric(entropy), x$Pos.Def))
names(results) <- c("ICA", "Entropy", "T0T1", "T0S0", "T0S1", "T1S0", "T1S1", "S0S1")
results_all <- results[order(results$Entropy, decreasing = TRUE),]
results <- results[order(results$Entropy, decreasing = TRUE),][1,]


fit <- 
  list(ICA.Max.Ent=as.numeric(results[1]), Max.Ent=as.numeric(results[2]), Entropy=as.numeric(entropy), 
       Table.ICA.Entropy=results_all, ICA.Fit=x, Call=match.call())

class(fit) <- "MaxEntContCont"
fit

}