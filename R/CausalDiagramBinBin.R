

#x = ICA; Min = 0.1; Max = 1; Monotonicity = "True.Endp"



CausalDiagramBinBin <- function(x, Values="Corrs", Theta_T0S0, Theta_T1S1, Min=0, Max=1, Cex.Letters=3, 
                                Cex.Corrs=2, 
                                Lines.Rel.Width=TRUE, Col.Pos.Neg=TRUE,
                                Monotonicity, Histograms.Correlations=FALSE, Densities.Correlations=FALSE) {
  
  if (class(x)!="ICA.BinBin") {stop("The function CausalDiagramBinBin should be applied to an object of class ICA.BinBin.")}
  
  if (missing(Theta_T0S0)) {Theta_T0S0 <- 1}
  if (missing(Theta_T1S1)) {Theta_T1S1 <- 1} 
  
  dat <- data.frame(cbind(x$Pi.Vectors, x$R2_H, x$Theta_T, x$Theta_S), stringsAsFactors = TRUE)
  sub <- dat[dat$x.R2_H >= Min & dat$x.R2_H <= Max,] 
  
  if ((Monotonicity=="No" | Monotonicity=="True.Endp" | Monotonicity=="Surr.Endp" | Monotonicity=="Surr.True.Endp")==FALSE){
    stop("The Monotonicity=... argument is not correctly specified. \n")
  }
  
  # Modifications if only one monotonicity setting was specified in ICA and causal diagram calls
  if (is.null(sub$Monotonicity)==TRUE){
    sub$Monotonicity <- x$Monotonicity
    if (all(sub$Monotonicity=="Surr.Endp")){sub$Monotonicity <- "Surr"; 
    sub$Pi_0010 <- sub$Pi_1010 <- sub$Pi_1110 <- sub$Pi_0110 <- 0}
    if (all(sub$Monotonicity=="True.Endp")){sub$Monotonicity <- "True"; 
    sub$Pi_1000 <- sub$Pi_1010 <- sub$Pi_1001 <- sub$Pi_1011 <- 0}
    if (all(sub$Monotonicity=="Surr.True.Endp")){sub$Monotonicity <- "SurrTrue"; 
    sub$Pi_0010 <- sub$Pi_1010 <- sub$Pi_1110 <- sub$Pi_0110 <- sub$Pi_1000 <- sub$Pi_1010 <- sub$Pi_1001 <- sub$Pi_1011 <- 0}
    }
   
  
  if (Monotonicity=="No"){sub <- sub[sub$Monotonicity=="No",]}
  if (Monotonicity=="True.Endp"){sub <- sub[c(sub$Monotonicity=="True"),]}
  if (Monotonicity=="Surr.Endp"){sub <- sub[c(sub$Monotonicity=="Surr"),]}
  if (Monotonicity=="Surr.True.Endp"){sub <- sub[c(sub$Monotonicity=="SurrTrue"),]}
  
  
  
  
  if (dim(sub)[1]==0) {stop("There are no valid observations within the specified R2_H range [min, max].")}
  
  cat("Note. The figure is based on ", dim(sub)[1], " observations. \n", sep="")
  

  if (Values == "Corrs"){
  
  #S0 T0
  cell_00 <- sub$Pi_0000 + sub$Pi_0100 + sub$Pi_0001 + sub$Pi_0101    
  cell_10 <- sub$Pi_0010 + sub$Pi_0110 + sub$Pi_0011 + sub$Pi_0111   
  cell_01 <- sub$Pi_1000 + sub$Pi_1001 + sub$Pi_1101 + sub$Pi_1100  
  cell_11 <- sub$Pi_1010 + sub$Pi_1110 + sub$Pi_1011 + sub$Pi_1111  
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  
  I <- (cell_00 * log2(cell_00/(sum0_ * sum_0)))+
    (cell_10 * log2(cell_10/(sum1_ * sum_0)))+
    (cell_01 * log2(cell_01/(sum0_ * sum_1)))+
    (cell_11 * log2(cell_11/(sum1_ * sum_1)))
  H_S <- - ((sum_0 * log2(sum_0)) + (sum_1 * log2(sum_1)))
  H_T <- - ((sum0_ * log2(sum0_)) + (sum1_ * log2(sum1_))) 
  R2_H_S0T0 <- 
    I/(min(H_S, H_T))
  
  
  #S1 T1
  cell_00 <- sub$Pi_0000 + sub$Pi_0010 + sub$Pi_1000 + sub$Pi_1010    
  cell_10 <- sub$Pi_0001 + sub$Pi_1001 + sub$Pi_1011 + sub$Pi_0011     
  cell_01 <- sub$Pi_0100 + sub$Pi_1110 + sub$Pi_0110 + sub$Pi_1100 
  cell_11 <- sub$Pi_0101 + sub$Pi_1101 + sub$Pi_1111 + sub$Pi_0111 
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  I <- (cell_00 * log2(cell_00/(sum0_ * sum_0)))+
    (cell_10 * log2(cell_10/(sum1_ * sum_0)))+
    (cell_01 * log2(cell_01/(sum0_ * sum_1)))+
    (cell_11 * log2(cell_11/(sum1_ * sum_1)))
  H_S <- - ((sum_0 * log2(sum_0)) + (sum_1 * log2(sum_1)))
  H_T <- - ((sum0_ * log2(sum0_)) + (sum1_ * log2(sum1_))) 
  R2_H_S1T1 <- 
    I/(min(H_S, H_T)) 
  
  #S0 S1
  cell_00 <- sub$Pi_0000 + sub$Pi_0100 + sub$Pi_1000 + sub$Pi_1100     
  cell_10 <- sub$Pi_0010 + sub$Pi_1010 + sub$Pi_1110 + sub$Pi_0110      
  cell_01 <- sub$Pi_0001 + sub$Pi_0101 + sub$Pi_1001 + sub$Pi_1101      
  cell_11 <- sub$Pi_1011 + sub$Pi_1111 + sub$Pi_0011 + sub$Pi_0111 
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  R2_H_S0S1 <- matrix(NA, ncol = length(sub$Monotonicity)) 
  for (i in 1: length(sub$Monotonicity)){
    if (sub$Monotonicity[i]=="No" | sub$Monotonicity[i]=="True"){
      I <- (cell_00[i] * log2(cell_00[i]/(sum0_[i] * sum_0[i])))+
        (cell_10[i] * log2(cell_10[i]/(sum1_[i] * sum_0[i])))+  #
        (cell_01[i] * log2(cell_01[i]/(sum0_[i] * sum_1[i])))+
        (cell_11[i] * log2(cell_11[i]/(sum1_[i] * sum_1[i])))
      H_S <- - ((sum_0[i] * log2(sum_0[i])) + (sum_1[i] * log2(sum_1[i])))
      H_T <- - ((sum0_[i] * log2(sum0_[i])) + (sum1_[i] * log2(sum1_[i]))) 
      R2_H_S0S1[i] <- 
        as.numeric(I/(min(H_S, H_T))) 
    }
    
    if (sub$Monotonicity[i]=="Surr" | sub$Monotonicity[i]=="SurrTrue"){
      I <- (cell_00[i] * log2(cell_00[i]/(sum0_[i] * sum_0[i])))+
        (cell_01[i] * log2(cell_01[i]/(sum0_[i] * sum_1[i])))+
        (cell_11[i] * log2(cell_11[i]/(sum1_[i] * sum_1[i])))
      H_S <- - ((sum_0[i] * log2(sum_0[i])) + (sum_1[i] * log2(sum_1[i])))
      H_T <- - ((sum0_[i] * log2(sum0_[i])) + (sum1_[i] * log2(sum1_[i]))) 
      R2_H_S0S1[i] <- 
        as.numeric(I/(min(H_S, H_T))) 
    }
  }  
  
  #T0 T1
  cell_00 <- sub$Pi_0000 + sub$Pi_0010 + sub$Pi_0001 + sub$Pi_0011          
  cell_10 <- sub$Pi_1000 + sub$Pi_1010 + sub$Pi_1001 + sub$Pi_1011     
  cell_01 <- sub$Pi_0100 + sub$Pi_0101 + sub$Pi_0110 + sub$Pi_0111         
  cell_11 <- sub$Pi_1110 + sub$Pi_1101 + sub$Pi_1111 + sub$Pi_1100  
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  R2_H_T0T1 <- matrix(NA, ncol = length(sub$Monotonicity))  
  
  for (i in 1: length(sub$Monotonicity)){
    
    if (sub$Monotonicity[i]=="No" | sub$Monotonicity[i]=="Surr"){
      
      I <- (cell_00[i] * log2(cell_00[i]/(sum0_[i] * sum_0[i])))+
        (cell_10[i] * log2(cell_10[i]/(sum1_[i] * sum_0[i])))+
        (cell_01[i] * log2(cell_01[i]/(sum0_[i] * sum_1[i])))+
        (cell_11[i] * log2(cell_11[i]/(sum1_[i] * sum_1[i])))  
      H_S <- - ((sum_0[i] * log2(sum_0[i])) + (sum_1[i] * log2(sum_1[i])))
      H_T <- - ((sum0_[i] * log2(sum0_[i])) + (sum1_[i] * log2(sum1_[i]))) 
      R2_H_T0T1[i] <- 
        I/(min(H_S, H_T)) 
    }
    
    if (sub$Monotonicity[i]=="True" | sub$Monotonicity[i]=="SurrTrue"){
      
      I <- (cell_00[i] * log2(cell_00[i]/(sum0_[i] * sum_0[i])))+
        (cell_01[i] * log2(cell_01[i]/(sum0_[i] * sum_1[i])))+
        (cell_11[i] * log2(cell_11[i]/(sum1_[i] * sum_1[i])))  
      H_S <- - ((sum_0[i] * log2(sum_0[i])) + (sum_1[i] * log2(sum_1[i])))
      H_T <- - ((sum0_[i] * log2(sum0_[i])) + (sum1_[i] * log2(sum1_[i]))) 
      R2_H_T0T1[i] <- 
        I/(min(H_S, H_T)) 
    }
    
  }
  
  
  #S0 T1
  cell_00 <- sub$Pi_0000 + sub$Pi_0001 + sub$Pi_1000 + sub$Pi_1001         
  cell_10 <- sub$Pi_0010 + sub$Pi_1010 + sub$Pi_1011 + sub$Pi_0011     
  cell_01 <- sub$Pi_0100 + sub$Pi_0101 + sub$Pi_1101 + sub$Pi_1100          
  cell_11 <- sub$Pi_1110 + sub$Pi_1111 + sub$Pi_0110 + sub$Pi_0111   
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  I <- (cell_00 * log2(cell_00/(sum0_ * sum_0)))+
    (cell_10 * log2(cell_10/(sum1_ * sum_0)))+
    (cell_01 * log2(cell_01/(sum0_ * sum_1)))+
    (cell_11 * log2(cell_11/(sum1_ * sum_1)))
  H_S <- - ((sum_0 * log2(sum_0)) + (sum_1 * log2(sum_1)))
  H_T <- - ((sum0_ * log2(sum0_)) + (sum1_ * log2(sum1_))) 
  R2_H_S0T1 <- 
    I/(min(H_S, H_T)) 
  
  
  #S1 T0
  cell_00 <- sub$Pi_0000 + sub$Pi_0100 + sub$Pi_0010 + sub$Pi_0110          
  cell_10 <- sub$Pi_0001 + sub$Pi_0101 + sub$Pi_0011 + sub$Pi_0111      
  cell_01 <- sub$Pi_1000 + sub$Pi_1010 + sub$Pi_1110 + sub$Pi_1100         
  cell_11 <- sub$Pi_1001 + sub$Pi_1101 + sub$Pi_1011 + sub$Pi_1111
  sum_0 <- cell_00 + cell_10  
  sum_1 <- cell_01 + cell_11
  sum0_ <- cell_00 + cell_01 
  sum1_ <- cell_10 + cell_11
  
  I <- (cell_00 * log2(cell_00/(sum0_ * sum_0)))+
    (cell_10 * log2(cell_10/(sum1_ * sum_0)))+
    (cell_01 * log2(cell_01/(sum0_ * sum_1)))+
    (cell_11 * log2(cell_11/(sum1_ * sum_1)))
  H_S <- - ((sum_0 * log2(sum_0)) + (sum_1 * log2(sum_1)))
  H_T <- - ((sum0_ * log2(sum0_)) + (sum1_ * log2(sum1_))) 
  R2_H_S1T0 <- 
    I/(min(H_S, H_T)) 
  
  med_T0T1 <- round(median(R2_H_T0T1, na.rm = T), digits=2) 
  med_T0S0 <- round(median(R2_H_S0T0), digits=2)
  med_T0S1 <- round(median(R2_H_S1T0), digits=2)
  med_T1S0 <- round(median(R2_H_S0T1), digits=2)
  med_T1S1 <- round(median(R2_H_S1T1), digits=2)
  med_S0S1 <- round(median(R2_H_S0S1, na.rm = T), digits=2)
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  plot(0:10, 0:10, axes=FALSE, xlab="", ylab="", type="n")  
  par(oma=c(0, 0, 0, 0))
  text(1, 9, expression(S[0]), cex=Cex.Letters)
  text(1, 1, expression(S[1]), cex=Cex.Letters)
  text(0.4, 5, med_S0S1, cex=Cex.Corrs)
  text(9, 9, expression(T[0]), cex=Cex.Letters)
  text(9, 1, expression(T[1]), cex=Cex.Letters)
  text(4, 6.9, med_T1S0, cex=Cex.Corrs)
  text(4, 3.1, med_T0S1, cex=Cex.Corrs)
  text(5, 9.5, med_T0S0, cex=Cex.Corrs)
  text(5, 0.5, med_T1S1, cex=Cex.Corrs)
  text(9.6, 5, med_S0S1, cex=Cex.Corrs)
  
  col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1
    segments(x0=1, y0=8, x1=1, y1=2, lwd=1, col=col_S0S1)
    segments(x0=1.5, y0=8, x1=8.5, y1=2, lwd=1, col=col_T1S0)
    segments(x0=1.5, y0=2, x1=8.5, y1=8, lwd=1, col=col_T0S1)
    segments(x0=1.5, y0=9, x1=8.5, y1=9, lwd=1, col=col_T0S0)
    segments(x0=1.5, y0=1, x1=8.5, y1=1, lwd=1, col=col_T1S1)
    segments(x0=9, y0=8, x1=9, y1=2, lwd=1,  col=col_T0T1)
       
  }
  
if  (Histograms.Correlations==TRUE){
  dev.new(width = 10, height = 10)
  par(mar = c(4.4, 2.5, 2, 1))
  par(mfrow=c(2, 2))
  hist(R2_H_T0T1, main=expression(paste(R[h]^2, "(", T[0], ",", T[1], ")")), xlab="", col="grey")
  hist(R2_H_S1T0, main=expression(paste(R[h]^2, "(", T[0], ",", S[1], ")")), xlab="", col="grey")
  hist(R2_H_S0T1, main=expression(paste(R[h]^2, "(", T[1], ",", S[0], ")")), xlab="", col="grey")
  hist(R2_H_S0S1, main=expression(paste(R[h]^2, "(", S[0], ",", S[1], ")")), xlab="", col="grey")
  par(mfrow=c(1, 1))
}  
 
  
if  (Densities.Correlations==TRUE){
    dev.new(width = 10, height = 10)
    par(mar = c(4.4, 2.5, 2, 1))
    par(mfrow=c(2, 2))
    plot(density(R2_H_T0T1), main=expression(paste(R[h]^2, "(", T[0], ",", T[1], ")")), xlab="", col=1)
    plot(density(R2_H_S1T0), main=expression(paste(R[h]^2, "(", T[0], ",", S[1], ")")), xlab="", col=1)
    plot(density(R2_H_S0T1), main=expression(paste(R[h]^2, "(", T[1], ",", S[0], ")")), xlab="", col=1)
    plot(density(R2_H_S0S1), main=expression(paste(R[h]^2, "(", S[0], ",", S[1], ")")), xlab="", col=1)
    par(mfrow=c(1, 1))
  }  
  
  
   
if (Values=="ORs"){
  
  # ORs
  pi_T_00 <- sub$Pi_0000 + sub$Pi_0010 + sub$Pi_0001 + sub$Pi_0011 
  pi_T_01 <- sub$Pi_0100 + sub$Pi_0101 + sub$Pi_0110 + sub$Pi_0111
  pi_T_10 <- sub$Pi_1000 + sub$Pi_1010 + sub$Pi_1001 + sub$Pi_1011
  pi_T_11 <- sub$Pi_1110 + sub$Pi_1101 + sub$Pi_1111 + sub$Pi_1100
  pi_S_00 <- sub$Pi_0000 + sub$Pi_0100 + sub$Pi_1000 + sub$Pi_1100  
  pi_S_01 <- sub$Pi_0001 + sub$Pi_0101 + sub$Pi_1001 + sub$Pi_1101 
  pi_S_10 <- sub$Pi_0010 + sub$Pi_1010 + sub$Pi_1110 + sub$Pi_0110   
  pi_S_11 <- sub$Pi_1011 + sub$Pi_1111 + sub$Pi_0011 + sub$Pi_0111
  
  Theta_T0T1 <- (pi_T_00 * pi_T_11)/(pi_T_10 * pi_T_01)
  Theta_S0S1 <- (pi_S_00 * pi_S_11)/(pi_S_10 * pi_S_01)
  
  p_00 <- sub$Pi_0000 + sub$Pi_0100 + sub$Pi_0010 + sub$Pi_0110
  p_01 <- sub$Pi_1000 + sub$Pi_1010 + sub$Pi_1110 + sub$Pi_1100
  p_10 <- sub$Pi_0001 + sub$Pi_0101 + sub$Pi_0011 + sub$Pi_0111
  p_11 <- sub$Pi_1001 + sub$Pi_1101 + sub$Pi_1011 + sub$Pi_1111 
  
  q_00 <- sub$Pi_0000 + sub$Pi_0001 + sub$Pi_1000 + sub$Pi_1001
  q_01 <- sub$Pi_0100 + sub$Pi_0101 + sub$Pi_1101 + sub$Pi_1100
  q_10 <- sub$Pi_0010 + sub$Pi_1010 + sub$Pi_1011 + sub$Pi_0011
  q_11 <- sub$Pi_1110 + sub$Pi_1111 + sub$Pi_0110 + sub$Pi_0111
  
  Theta_T0S1 <-(p_00 * p_11) / (p_01 * p_10)
  Theta_T1S0 <-(q_00 * q_11) / (q_01 * q_10)
  
  med_T0T1 <- round(median(Theta_T0T1), digits=2)
  med_T0S0 <- round(median(Theta_T0S0), digits=2)
  med_T0S1 <- round(median(Theta_T0S1), digits=2)
  med_T1S0 <- round(median(Theta_T1S0), digits=2)
  med_T1S1 <- round(median(Theta_T1S1), digits=2)
  med_S0S1 <- round(median(Theta_S0S1), digits=2)
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  plot(0:10, 0:10, axes=FALSE, xlab="", ylab="", type="n")  
  par(oma=c(0, 0, 0, 0))
  text(1, 9, expression(S[0]), cex=Cex.Letters)
  text(1, 1, expression(S[1]), cex=Cex.Letters)
  text(0.4, 5, med_S0S1, cex=Cex.Corrs)
  text(9, 9, expression(T[0]), cex=Cex.Letters)
  text(9, 1, expression(T[1]), cex=Cex.Letters)
  text(4, 6.9, med_T1S0, cex=Cex.Corrs)
  text(4, 3.1, med_T0S1, cex=Cex.Corrs)
  text(5, 9.5, med_T0S0, cex=Cex.Corrs)
  text(5, 0.5, med_T1S1, cex=Cex.Corrs)
  text(9.6, 5, med_T0T1, cex=Cex.Corrs)
  
  
  if (Lines.Rel.Width==TRUE){
    
    max_med <- max(med_S0S1[is.finite(med_S0S1)==TRUE], med_T1S0[is.finite(med_T1S0)==TRUE], med_T0S1[is.finite(med_T0S1)==TRUE],
                   med_T0S0[is.finite(med_T0S0)==TRUE], med_T1S1[is.finite(med_T1S1)==TRUE], med_T0T1[is.finite(med_T0T1)==TRUE])
    
    if (Col.Pos.Neg==FALSE) {col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1}
    
    if (Col.Pos.Neg==TRUE) {
      col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1
      if (med_S0S1<1) {col_S0S1 <- "red"} 
      if (med_T1S0<1) {col_T1S0 <- "red"}
      if (med_T0S1<1) {col_T0S1 <- "red"}
      if (med_T0S0<1) {col_T0S0 <- "red"} 
      if (med_T1S1<1) {col_T1S1 <- "red"} 
      if (med_T0T1<1) {col_T0T1 <- "red"}
    }
    
    if (med_S0S1=="Inf") {segments(x0=1, y0=8, x1=1, y1=2, lwd=1+2, col=col_S0S1)}
    if (med_S0S1!="Inf") {segments(x0=1, y0=8, x1=1, y1=2, lwd=1+med_S0S1/max_med, col=col_S0S1)}
    segments(x0=1.5, y0=8, x1=8.5, y1=2, lwd=1+med_T1S0/max_med, col=col_T1S0)
    segments(x0=1.5, y0=2, x1=8.5, y1=8, lwd=1+med_T0S1/max_med, col=col_T0S1)
    segments(x0=1.5, y0=9, x1=8.5, y1=9, lwd=1+med_T0S0/max_med, col=col_T0S0)
    segments(x0=1.5, y0=1, x1=8.5, y1=1, lwd=1+med_T1S1/max_med, col=col_T1S1)
    if (med_T0T1=="Inf") {segments(x0=9, y0=8, x1=9, y1=2, lwd=1+2,  col=col_T0T1)}
    if (med_T0T1!="Inf") {segments(x0=9, y0=8, x1=9, y1=2, lwd=1+med_T0T1/max_med,  col=col_T0T1)}
  }
  
  if (Lines.Rel.Width==FALSE){
    
    if (Col.Pos.Neg==FALSE) {col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1}
    
    if (Col.Pos.Neg==TRUE) {
      col_S0S1 <- col_T1S0 <- col_T0S1 <- col_T0S0 <- col_T1S1 <- col_T0T1 <- 1
      if (med_S0S1<0) {col_S0S1 <- "red"} 
      if (med_T1S0<0) {col_T1S0 <- "red"}
      if (med_T0S1<0) {col_T0S1 <- "red"}
      if (med_T0S0<0) {col_T0S0 <- "red"} 
      if (med_T1S1<0) {col_T1S1 <- "red"} 
      if (med_T0T1<0) {col_T0T1 <- "red"}
    }
    
    if (med_S0S1=="Inf") {segments(x0=1, y0=8, x1=1, y1=2, lwd=1, col=col_S0S1)}
    if (med_S0S1!="Inf") {segments(x0=1, y0=8, x1=1, y1=2, lwd=1, col=col_S0S1)}
    segments(x0=1.5, y0=8, x1=8.5, y1=2, lwd=1, col=col_T1S0)
    segments(x0=1.5, y0=2, x1=8.5, y1=8, lwd=1, col=col_T0S1)
    segments(x0=1.5, y0=9, x1=8.5, y1=9, lwd=1, col=col_T0S0)
    segments(x0=1.5, y0=1, x1=8.5, y1=1, lwd=1, col=col_T1S1)
    if (med_T0T1=="Inf") {segments(x0=9, y0=8, x1=9, y1=2, lwd=1,  col=col_T0T1)}
    if (med_T0T1!="Inf") {segments(x0=9, y0=8, x1=9, y1=2, lwd=1,  col=col_T0T1)}
    
  }
  
  }
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  
  if  (Histograms.Correlations==TRUE & Values == "Corrs"){
  fit <- 
    list(R2_H_T0T1=as.numeric(R2_H_T0T1), R2_H_S1T0=as.numeric(R2_H_S1T0), R2_H_S0T1=as.numeric(R2_H_S0T1),
         R2_H_S0S1=as.numeric(R2_H_S0S1), R2_H_S0T0=as.numeric(R2_H_S0T0), R2_H_S1T1=as.numeric(R2_H_S1T1),
         Call=match.call())   
  }
  } 


