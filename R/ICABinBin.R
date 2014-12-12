ICA.BinBin <- function(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1, Monotonicity=c("General"), 
                       Sum_Pi_f = seq(from=0.01, to=0.99, by=.01), M=10000, Seed=sample(1:100000, size=1)) {   
  
  Seed.orig <- Seed
  aantal <- length(Sum_Pi_f)*M 
  set.seed(Seed.orig)
  Seed <- round(runif(min=(1+aantal), max=2147483647, n=1), digits=0)
  Seed.generalCase <- Seed   #for gen case LS
  
  vector_b <- matrix(data=c(1, pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1), ncol=1)  
  
  sum_pi_f <- Pi.Vectors <- pi_f_all <- pi_r_all <- C3_all <- R2_H_all <- theta_T_all <- theta_S_all <- H_Delta_T_all <- pi_all <- NULL
  
  if (Monotonicity=="General"){
    
    T_No <- ICA.BinBin(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1, Seed=Seed.orig, Monotonicity=c("No"), Sum_Pi_f, M)
    sum_pi_f <- Pi.Vectors <- pi_f_all <- pi_r_all <- C3_all <- R2_H_all <- theta_T_all <- theta_S_all <- H_Delta_T_all <- pi_all <- NULL
    T_True <- ICA.BinBin(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1, Seed=Seed.orig, Monotonicity=c("True.Endp"), Sum_Pi_f, M)
    sum_pi_f <- Pi.Vectors <- pi_f_all <- pi_r_all <- C3_all <- R2_H_all <- theta_T_all <- theta_S_all <- H_Delta_T_all <- pi_all <- NULL
    T_Surr <- ICA.BinBin(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1, Seed=Seed.orig, Monotonicity=c("Surr.Endp"), Sum_Pi_f, M)
    sum_pi_f <- Pi.Vectors <- pi_f_all <- pi_r_all <- C3_all <- R2_H_all <- theta_T_all <- theta_S_all <- H_Delta_T_all <- pi_all <- NULL
    T_SurrTrue <- ICA.BinBin(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1, Seed=Seed.orig, Monotonicity=c("Surr.True.Endp"), Sum_Pi_f, M)
    sum_pi_f <- Pi.Vectors <- pi_f_all <- pi_r_all <- C3_all <- R2_H_all <- theta_T_all <- theta_S_all <- H_Delta_T_all <- pi_all <- NULL
    
    try(Pi.Vectors.No <- cbind(T_No$Pi.Vectors, "No"), silent=TRUE)
    try(Pi.Vectors.T <- cbind(T_True$Pi.Vectors[,1:5], 0, 0, 0, T_True$Pi.Vectors[,6:7], 0, T_True$Pi.Vectors[,8:13], "True"), silent=TRUE)
    try(Pi.Vectors.S <- cbind(T_Surr$Pi.Vectors[,1:2], 0, T_Surr$Pi.Vectors[,3:5], 0, T_Surr$Pi.Vectors[,8], 0, T_Surr$Pi.Vectors[,9],
                              T_Surr$Pi.Vectors[,6], T_Surr$Pi.Vectors[,10], 0, T_Surr$Pi.Vectors[,7], 
                              T_Surr$Pi.Vectors[,11:13], "Surr"), silent=TRUE)
    try(Pi.Vectors.ST <- cbind(T_SurrTrue$Pi.Vectors[,1:2], 0, T_SurrTrue$Pi.Vectors[,3:4], 0, 0, 0, 0, T_SurrTrue$Pi.Vectors[,5], 0,
                               T_SurrTrue$Pi.Vectors[,6], 0, T_SurrTrue$Pi.Vectors[,7:10], "SurrTrue"), silent=TRUE)
    
    try(colnames(Pi.Vectors.No) <- colnames(Pi.Vectors.T) <-  colnames(Pi.Vectors.S) <- colnames(Pi.Vectors.ST) <- 
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity"), silent=TRUE) 
    
    if (exists("Pi.Vectors.No")==FALSE){Pi.Vectors.No <- NULL}
    if (exists("Pi.Vectors.T")==FALSE){Pi.Vectors.T <- NULL}
    if (exists("Pi.Vectors.S")==FALSE){Pi.Vectors.S <- NULL}
    if (exists("Pi.Vectors.ST")==FALSE){Pi.Vectors.ST <- NULL}
    
    try(Pi.Vectors <- rbind(Pi.Vectors.No, Pi.Vectors.T, Pi.Vectors.S, Pi.Vectors.ST), silent=TRUE)
    C3_all <- c(T_No$C3, T_True$C3, T_Surr$C3, T_SurrTrue$C3)
    R2_H_all  <- c(T_No$R2_H, T_True$R2_H, T_Surr$R2_H, T_SurrTrue$R2_H)
    theta_T_all <- c(T_No$Theta_T, T_True$Theta_T, T_Surr$Theta_T, T_SurrTrue$Theta_T)
    theta_S_all <- c(T_No$Theta_S, T_True$Theta_S, T_Surr$Theta_S, T_SurrTrue$Theta_S)
    H_Delta_T_all <- c(T_No$H_Delta_T, T_True$H_Delta_T, T_Surr$H_Delta_T, T_SurrTrue$H_Delta_T)  
    
  }
  
  if (Monotonicity=="No"){
    
    A_r <- matrix(data=c(1, 0, 0, 0, 0, 0, 0,
                         1, 0, 0, 0, 1, 0, 0,
                         1, 0, 0, 0, 0, 1, 0,
                         1, 0, 0, 0, 0, 0, 1,
                         1, 0, 0, 1, 0, 0, 0,
                         1, 0, 1, 0, 0, 0, 0,
                         1, 1, 0, 0, 0, 0, 0), ncol=7)
    A_f <- matrix(data=c(1, 0, 1, 0, 0, 0, 1,
                         1, 1, 0, 0, 1, 0, 0,
                         1, 0, 1, 1, 0, 0, 0,
                         1, 1, 0, 0, 0, 0, 1,
                         1, 1, 0, 1, 0, 0, 0,
                         1, 0, 0, 0, 1, 1, 0,
                         1, 0, 0, 0, 0, 1, 1,
                         1, 0, 0, 1, 0, 1, 0,
                         1, 0, 1, 0, 1, 0, 0), ncol=9)
    A <- cbind(A_r, A_f)
    
    
    #restrictions
    min_pi_1001 <- min(pi1_0_, pi_0_1)
    min_pi_1110 <- min(pi1_1_, pi_1_0)
    min_pi_1101 <- min(pi1_0_, pi_1_1)
    min_pi_1011 <- min(pi1_1_, pi_0_1)
    min_pi_1111 <- min(pi1_1_, pi_1_1)
    min_pi_0110 <- min(pi0_1_, pi_1_0)
    min_pi_0011 <- min(pi0_1_, pi_0_1)
    min_pi_0111 <- min(pi0_1_, pi_1_1)
    min_pi_1100 <- min(pi1_0_, pi_1_0)
    
    max_val <- min((min_pi_1001 + min_pi_1110 + min_pi_1101 + min_pi_1011 + min_pi_1111 + 
      min_pi_0110 + min_pi_0011 + min_pi_0111 + min_pi_1100), 1)
    min_val <- min(Sum_Pi_f)
    by_val = (max(Sum_Pi_f)-min(Sum_Pi_f))/(length(Sum_Pi_f)-1)
        Sum_Pi_f <- seq(min_val, max_val, by=by_val)
    
    Sum_Pi_f_number <- length(Sum_Pi_f)
    
    for (k in 1: Sum_Pi_f_number){
      
      for (i in 1: M){

        Seed <- Seed-1   
        set.seed(Seed)
        
        pi_f <- (RandVec(a=0, b=Sum_Pi_f[k], s=Sum_Pi_f[k], n=9, 1, Seed=Seed))$RandVecOutput  
        
        pi_r <- solve(A_r) %*% (vector_b - (A_f %*% pi_f))
        
        if ((sum(pi_r >= -0.00000001 & pi_r <= 1.00000001) == 7)==TRUE) {
          
          for (l in 1: length(pi_r)){
            if (pi_r[l] < 0) {pi_r[l] <- c(0)}
            if (pi_r[l] > 1) {pi_r[l] <- c(1)}
          }
          
          pi <- rbind(pi_r, pi_f)
          
          mat1 <- pi[7] 
          mat2 <- pi[3] + pi[9] 
          mat3 <- pi[13] 
          mat4 <- pi[6] + pi[11] 
          mat5 <- pi[1] + pi[14] + pi[16] + pi[12] 
          mat6 <- pi[2] + pi[15]
          mat7 <- pi[8]
          mat8 <- pi[4] + pi[10] 
          mat9 <- pi[5]
          
          Delta_c_mat <-
            matrix(data=c(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9), nrow=3)
          
          C3 <- 1
          C3_all <- cbind(C3_all, C3)
          
          #R2_H
          sum_S_min1 <- mat1+mat2+mat3
          sum_S_0 <- mat4+mat5+mat6
          sum_S_1 <- mat7+mat8+mat9
          
          sum_T_min1 <- mat1+mat4+mat7
          sum_T_0 <- mat2+mat5+mat8
          sum_T_1 <- mat3+mat6+mat9
          
          I_Delta_T_Delta_S <- 
            (mat1*log2(mat1/(sum_S_min1*sum_T_min1)))+(mat2*log2(mat2/(sum_S_min1*sum_T_0)))+(mat3*log2(mat3/(sum_S_min1*sum_T_1)))+
            (mat4*log2(mat4/(sum_S_0*sum_T_min1)))+(mat5*log2(mat5/(sum_S_0*sum_T_0)))+(mat6*log2(mat6/(sum_S_0*sum_T_1)))+
            (mat7*log2(mat7/(sum_S_1*sum_T_min1)))+(mat8*log2(mat8/(sum_S_1*sum_T_0)))+(mat9*log2(mat9/(sum_S_1*sum_T_1)))
          
          H_Delta_T <-  
            -(((mat1+mat4+mat7)*log2(mat1+mat4+mat7))+ 
                ((mat2+mat5+mat8)*log2(mat2+mat5+mat8))+
                ((mat3+mat6+mat9)*log2(mat3+mat6+mat9)))
          H_Delta_T_all <- cbind(H_Delta_T_all, H_Delta_T)
          H_Delta_S <-   
            -(((mat1+mat2+mat3)*log2(mat1+mat2+mat3))+ 
                ((mat4+mat5+mat6)*log2(mat4+mat5+mat6))+
                ((mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
          R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
          R2_H_all <- cbind(R2_H_all, R2_H)
          
          
          # Association potential outcomes true (theta_T) and surrogate (theta_S)
          pi_T_00 <- pi[1] + pi[3] + pi[4] + pi[14] 
          pi_T_01 <- pi[2] + pi[5] + pi[13] + pi[15]
          pi_T_10 <- pi[6] + pi[7] + pi[8] + pi[11]
          pi_T_11 <- pi[9] + pi[10] + pi[12] + pi[16]
          
          pi_S_00 <- pi[1] + pi[2] + pi[6] + pi[16]  
          pi_S_01 <- pi[4] + pi[5] + pi[8] + pi[10] 
          pi_S_10 <- pi[3] + pi[7] + pi[9] + pi[13]   
          pi_S_11 <- pi[11] + pi[12] + pi[14] + pi[15]
          
          theta_T <- (pi_T_00 * pi_T_11)/(pi_T_10 * pi_T_01)
          theta_S <- (pi_S_00 * pi_S_11)/(pi_S_10 * pi_S_01)
          theta_T_all <- cbind(theta_T_all, theta_T)  
          theta_S_all <- cbind(theta_S_all, theta_S)
          
          pi_all <- cbind(pi_all, pi) 
        }
      }
    }
    
    num <- dim(R2_H_all)[2]
   
    if (is.null(num)==TRUE){stop(c("No pi vectors found that are in line with the restrictions based on the data."))}
    
    sum_pi_f <- colSums(pi_all[8:16,]) 
    Pi.Vectors <- 
      data.frame(t(rbind(pi_all, sum_pi_f))) 
    colnames(Pi.Vectors) <- c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", 
                              "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f") 

    
    
  } # end no mono
  
   
  if (Monotonicity=="True.Endp"){
    
    A_star_r <- matrix(data=c(1, 0, 0, 0, 0, 0, 0,
                              1, 0, 0, 0, 1, 0, 0,
                              1, 0, 0, 0, 0, 1, 0,
                              1, 0, 0, 0, 0, 0, 1,
                              1, 0, 0, 1, 0, 0, 0,
                              1, 1, 0, 0, 1, 0, 0,
                              1, 0, 1, 1, 0, 0, 0), ncol=7)
    A_star_f <- matrix(data=c(1, 1, 0, 1, 0, 0, 0,
                              1, 0, 0, 0, 1, 1, 0,
                              1, 0, 0, 0, 0, 1, 1,
                              1, 0, 0, 1, 0, 1, 0,
                              1, 0, 1, 0, 1, 0, 0), ncol=5)
    A_star <- cbind(A_star_r, A_star_f)
    
    pi_f_all <- pi_r_all <- C3_all <- R2_H_all <- theta_T_all <- theta_S_all <- H_Delta_T_all <- NULL
    
    #restrictions
    min_pi_1111 <- min(pi1_1_, pi_1_1)
    min_pi_0110 <- min(pi0_1_, pi_1_0)
    min_pi_0011 <- min(pi0_1_, pi_0_1)
    min_pi_0111 <- min(pi0_1_, pi_1_1)
    min_pi_1100 <- min(pi1_0_, pi_1_0)
    
    max_val <- min((min_pi_1111 + min_pi_0110 + min_pi_0011 + min_pi_0111 + min_pi_1100), 1)
    min_val <- min(Sum_Pi_f)
    by_val = (max(Sum_Pi_f)-min(Sum_Pi_f))/(length(Sum_Pi_f)-1)
    Sum_Pi_f <- seq(min_val, max_val, by=by_val)
    
    Sum_Pi_f_number <- length(Sum_Pi_f)
    
    
    for (k in 1: Sum_Pi_f_number){
      
      for (i in 1: M){
        
        Seed <- Seed-1   
        set.seed(Seed)
        
        pi_f <- (RandVec(a=0, b=Sum_Pi_f[k], s=Sum_Pi_f[k], n=5, 1, Seed=Seed))$RandVecOutput  
        
        pi_r <-   
          solve(A_star_r) %*% (vector_b - (A_star_f %*% pi_f))
        
        if ((sum(pi_r >= -0.00000001 & pi_r <= 1.00000001) == 7)==TRUE) {
          
          for (l in 1: length(pi_r)){
            if (pi_r[l] < 0) {pi_r[l] <- c(0)}
            if (pi_r[l] > 1) {pi_r[l] <- c(1)}
          }
          
          pi <- rbind(pi_r, pi_f)
          
          mat1 <- c(0) 
          mat2 <- pi[3] + pi[6] 
          mat3 <- pi[9] 
          mat4 <- c(0) + c(0)  
          mat5 <- pi[1] + pi[10] + pi[12] + pi[8] 
          mat6 <- pi[2] + pi[11]
          mat7 <- c(0)
          mat8 <- pi[4] + pi[7] 
          mat9 <- pi[5]
          
          Delta_c_mat <-
            matrix(data=c(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9), nrow=3)
          
          C3 <- 1
          C3_all <- cbind(C3_all, C3)
          
          #R2_H
          sum_S_min1 <- mat1+mat2+mat3
          sum_S_0 <- mat4+mat5+mat6
          sum_S_1 <- mat7+mat8+mat9
          
          sum_T_min1 <- mat1+mat4+mat7
          sum_T_0 <- mat2+mat5+mat8
          sum_T_1 <- mat3+mat6+mat9
          
          I_Delta_T_Delta_S <- 
            0+(mat2*log2(mat2/(sum_S_min1*sum_T_0)))+(mat3*log2(mat3/(sum_S_min1*sum_T_1)))+
            0+(mat5*log2(mat5/(sum_S_0*sum_T_0)))+(mat6*log2(mat6/(sum_S_0*sum_T_1)))+
            0+(mat8*log2(mat8/(sum_S_1*sum_T_0)))+(mat9*log2(mat9/(sum_S_1*sum_T_1)))
          
          H_Delta_T <-  
            -(0+((mat2+mat5+mat8)*log2(mat2+mat5+mat8))+((mat3+mat6+mat9)*log2(mat3+mat6+mat9)))
          H_Delta_T_all <- cbind(H_Delta_T_all, H_Delta_T)
          H_Delta_S <-  -(((mat1+mat2+mat3)*log2(mat1+mat2+mat3))+((mat4+mat5+mat6)*log2(mat4+mat5+mat6))+((mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
          R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
          R2_H_all <- cbind(R2_H_all, R2_H)
          
          
          # Association potential outcomes true (theta_T) and surrogate (theta_S)
          
          pi_T_00 <- pi[1] + pi[3] + pi[4] + pi[10] 
          pi_T_01 <- pi[2] + pi[5] + pi[9] + pi[11]
          pi_T_10 <- 0 + 0 + 0 + 0
          pi_T_11 <- pi[6] + pi[7] + pi[8] + pi[12]
          
          
          pi_S_00 <- pi[1] + pi[2] + 0 + pi[12]  
          pi_S_01 <- pi[4] + pi[5] + 0 + pi[7] 
          pi_S_10 <- pi[3] + 0 + pi[6] + pi[9]   
          pi_S_11 <- 0 + pi[8] + pi[10] + pi[11]
          
          theta_T <- (pi_T_00 * pi_T_11)/(pi_T_10 * pi_T_01)
          theta_S <- (pi_S_00 * pi_S_11)/(pi_S_10 * pi_S_01)
          theta_T_all <- cbind(theta_T_all, theta_T)  
          theta_S_all <- cbind(theta_S_all, theta_S)
          
          pi_all <- cbind(pi_all, pi) 
          
        }
      }
    } 
    
    num <- dim(R2_H_all)[2]
    
    
    if (is.null(num)==TRUE){stop(c("No pi vectors found that are in line with the restrictions based on the data."))}
    
    sum_pi_f <- colSums(pi_all[8:12,]) 
    Pi.Vectors <- 
      data.frame(t(rbind(pi_all, sum_pi_f))) 
    colnames(Pi.Vectors) <- c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", 
                              "Pi_1110", "Pi_1101", "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f")  
    
    
  } # end mono T
  
  
  
  if (Monotonicity=="Surr.Endp"){
    
    A_star_r <- matrix(data=c(1, 0, 0, 0, 0, 0, 0,
                              1, 0, 0, 0, 1, 0, 0,
                              1, 0, 0, 0, 0, 0, 1,
                              1, 0, 0, 1, 0, 0, 0,
                              1, 0, 1, 0, 0, 0, 0,
                              1, 1, 0, 0, 0, 0, 1,
                              1, 0, 0, 0, 0, 1, 1), ncol=7)
    A_star_f <- matrix(data=c(1, 0, 1, 0, 0, 0, 1,
                              1, 0, 1, 1, 0, 0, 0,
                              1, 1, 0, 1, 0, 0, 0,
                              1, 0, 0, 1, 0, 1, 0,
                              1, 0, 1, 0, 1, 0, 0), ncol=5)
    A_star <- cbind(A_star_r, A_star_f)
    
    #restrictions
    min_pi_1001 <- min(pi1_0_, pi_0_1)
    min_pi_1101 <- min(pi1_0_, pi_1_1)
    min_pi_1111 <- min(pi1_1_, pi_1_1)
    min_pi_0111 <- min(pi0_1_, pi_1_1)
    min_pi_1100 <- min(pi1_0_, pi_1_0)
    
    max_val <- min((min_pi_1001 + min_pi_1101 + min_pi_1111 + min_pi_0111 + min_pi_1100), 1)
    min_val <- min(Sum_Pi_f)
    by_val = (max(Sum_Pi_f)-min(Sum_Pi_f))/(length(Sum_Pi_f)-1)
    Sum_Pi_f <- seq(min_val, max_val, by=by_val)
    
    Sum_Pi_f_number <- length(Sum_Pi_f)
    
    
    for (k in 1: Sum_Pi_f_number){
      
      for (i in 1: M){
        
        Seed <- Seed-1  
        set.seed(Seed)
        
        pi_f <- (RandVec(a=0, b=Sum_Pi_f[k], s=Sum_Pi_f[k], n=5, 1, Seed=Seed))$RandVecOutput  
        
        pi_r <- solve(A_star_r) %*% (vector_b - (A_star_f %*% pi_f))
        
        if ((sum(pi_r >= -0.00000001 & pi_r <= 1.00000001) == 7)==TRUE) {
          
          for (l in 1: length(pi_r)){
            if (pi_r[l] < 0) {pi_r[l] <- c(0)}
            if (pi_r[l] > 1) {pi_r[l] <- c(1)}
          }      
          
          pi <- rbind(pi_r, pi_f)
          
          mat1 <- 0 
          mat2 <- 0 
          mat3 <- 0 
          mat4 <- pi[5] + pi[6] 
          mat5 <- pi[1] + pi[7] + pi[12] + pi[10] 
          mat6 <- pi[2] + pi[11]
          mat7 <- pi[8]
          mat8 <- pi[3] + pi[9] 
          mat9 <- pi[4]
          
          
          Delta_c_mat <-
            matrix(data=c(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9), nrow=3)
          
          C3 <- 1
          C3_all <- cbind(C3_all, C3)
          
          #R2_H
          sum_S_min1 <- mat1+mat2+mat3
          sum_S_0 <- mat4+mat5+mat6
          sum_S_1 <- mat7+mat8+mat9
          
          sum_T_min1 <- mat1+mat4+mat7
          sum_T_0 <- mat2+mat5+mat8
          sum_T_1 <- mat3+mat6+mat9
          
          I_Delta_T_Delta_S <- 
            0+(mat4*log2(mat4/(sum_S_0*sum_T_min1)))+(mat5*log2(mat5/(sum_S_0*sum_T_0)))+(mat6*log2(mat6/(sum_S_0*sum_T_1)))+
            (mat7*log2(mat7/(sum_S_1*sum_T_min1)))+(mat8*log2(mat8/(sum_S_1*sum_T_0)))+(mat9*log2(mat9/(sum_S_1*sum_T_1)))
          
          H_Delta_T <-  
            -(((mat1+mat4+mat7)*log2(mat1+mat4+mat7))+ 
                ((mat2+mat5+mat8)*log2(mat2+mat5+mat8))+
                ((mat3+mat6+mat9)*log2(mat3+mat6+mat9)))
          H_Delta_T_all <- cbind(H_Delta_T_all, H_Delta_T)
          H_Delta_S <-   
            -(0+((mat4+mat5+mat6)*log2(mat4+mat5+mat6))+((mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
          R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
          R2_H_all <- cbind(R2_H_all, R2_H)
          
          
          # Association potential outcomes true (theta_T) and surrogate (theta_S)
          pi_T_00 <- pi[1] + 0 + pi[3] + pi[7] 
          pi_T_01 <- pi[2] + pi[4] + 0 + pi[11]
          pi_T_10 <- pi[5] + 0 + pi[8] + pi[6]
          pi_T_11 <- 0 + pi[9] + pi[10] + pi[12]
          
          pi_S_00 <- pi[1] + pi[2] + pi[5] + pi[12]  
          pi_S_01 <- pi[3] + pi[4] + pi[8] + pi[10] 
          pi_S_10 <- 0 + 0 + 0 + 0   
          pi_S_11 <- pi[6] + pi[10] + pi[7] + pi[11]
          
          theta_T <- (pi_T_00 * pi_T_11)/(pi_T_10 * pi_T_01)
          theta_S <- (pi_S_00 * pi_S_11)/(pi_S_10 * pi_S_01)
          theta_T_all <- cbind(theta_T_all, theta_T)  
          theta_S_all <- cbind(theta_S_all, theta_S)
          
          pi_all <- cbind(pi_all, pi) 
          
          
        }
      }
    }
    
    num <- dim(R2_H_all)[2]
    
   
      if (is.null(num)==TRUE){stop(c("No pi vectors found that are in line with the restrictions based on the data."))}
    
    sum_pi_f <- colSums(pi_all[8:12,]) 
    
    Pi.Vectors <- 
      data.frame(t(rbind(pi_all, sum_pi_f))) 
    colnames(Pi.Vectors) <- c("Pi_0000", "Pi_0100", "Pi_0001", "Pi_0101", "Pi_1000",  
                              "Pi_1011", "Pi_0011", "Pi_1001", "Pi_1101", "Pi_1111", "Pi_0111", "Pi_1100", "Sum.Pi.f")   
    
   
  } # end mono S  
  
  if (Monotonicity=="Surr.True.Endp"){  
    
    A_star_r <- matrix(data=c(1, 0, 0, 0, 0, 0, 0,
                              1, 0, 0, 0, 1, 0, 0,
                              1, 0, 0, 0, 0, 0, 1,
                              1, 0, 0, 1, 0, 0, 0,
                              1, 0, 1, 1, 0, 0, 0,
                              1, 1, 0, 1, 0, 0, 0, 
                              1, 0, 0, 0, 0, 1, 1), ncol=7)
    A_star_f <- matrix(data=c(1, 0, 0, 1, 0, 1, 0,
                              1, 0, 1, 0, 1, 0, 0), ncol=2)
    
    A_star <- cbind(A_star_r, A_star_f)
    
    
    #restrictions
    min_pi_0111 <- min(pi0_1_, pi_1_1)
    min_pi_1100 <- min(pi1_0_, pi_1_0)
    
    max_val <- min((min_pi_0111 + min_pi_1100), 1)
    min_val <- min(Sum_Pi_f)
    by_val = (max(Sum_Pi_f)-min(Sum_Pi_f))/(length(Sum_Pi_f)-1)
    Sum_Pi_f <- seq(min_val, max_val, by=by_val)
    
    Sum_Pi_f_number <- length(Sum_Pi_f)
    
    for (k in 1: Sum_Pi_f_number){
      
      for (i in 1: M){
        
        Seed <- Seed-1   
        set.seed(Seed)
        
        pi_f <- (RandVec(a=0, b=Sum_Pi_f[k], s=Sum_Pi_f[k], n=2, 1, Seed=Seed))$RandVecOutput  
        
        pi_r <-   
          solve(A_star_r) %*% (vector_b - (A_star_f %*% pi_f))
        
        if ((sum(pi_r >= -0.00000001 & pi_r <= 1.00000001) == 7)==TRUE) {
          
          for (l in 1: length(pi_r)){
            if (pi_r[l] < 0) {pi_r[l] <- c(0)}
            if (pi_r[l] > 1) {pi_r[l] <- c(1)}
          }
          
          pi <- rbind(pi_r, pi_f)
          
          mat1 <- 0 
          mat2 <- 0 + 0 
          mat3 <- 0 
          mat4 <- 0 + 0 
          mat5 <- pi[1] + pi[7] + pi[9] + pi[6] 
          mat6 <- pi[2] + pi[8]
          mat7 <- 0
          mat8 <- pi[3] + pi[5] 
          mat9 <- pi[4]
          
          Delta_c_mat <-
            matrix(data=c(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9), nrow=3)
          
          C3 <- 1
          C3_all <- cbind(C3_all, C3)
          
          #R2_H
          sum_S_min1 <- mat1+mat2+mat3
          sum_S_0 <- mat4+mat5+mat6
          sum_S_1 <- mat7+mat8+mat9
          
          sum_T_min1 <- mat1+mat4+mat7
          sum_T_0 <- mat2+mat5+mat8
          sum_T_1 <- mat3+mat6+mat9
          
          I_Delta_T_Delta_S <- 
            0+0+(mat5*log2(mat5/(sum_S_0*sum_T_0)))+(mat6*log2(mat6/(sum_S_0*sum_T_1)))+
            0+(mat8*log2(mat8/(sum_S_1*sum_T_0)))+(mat9*log2(mat9/(sum_S_1*sum_T_1)))
          
          H_Delta_T <-  
            -(0+((mat2+mat5+mat8)*log2(mat2+mat5+mat8))+((mat3+mat6+mat9)*log2(mat3+mat6+mat9)))
          H_Delta_T_all <- cbind(H_Delta_T_all, H_Delta_T)
          H_Delta_S <-   
            -(0+((mat4+mat5+mat6)*log2(mat4+mat5+mat6))+
                ((mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
          R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
          R2_H_all <- cbind(R2_H_all, R2_H)
          
          
          # Association potential outcomes true (theta_T) and surrogate (theta_S)
          pi_T_00 <- pi[1] + 0 + pi[3] + pi[7] 
          pi_T_01 <- pi[2] + pi[4] + 0 + pi[8]
          pi_T_10 <- 0 + 0 + 0 + 0
          pi_T_11 <- 0 + pi[5] + pi[6] + pi[9]
          
          pi_S_00 <- pi[1] + pi[2] + 0 + pi[9]  
          pi_S_01 <- pi[3] + pi[4] + 0 + pi[5] 
          pi_S_10 <- 0 + 0 + 0 + 0   
          pi_S_11 <- 0 + pi[6] + pi[7] + pi[8]
          
          theta_T <- (pi_T_00 * pi_T_11)/(pi_T_10 * pi_T_01)
          theta_S <- (pi_S_00 * pi_S_11)/(pi_S_10 * pi_S_01)
          theta_T_all <- cbind(theta_T_all, theta_T)  
          theta_S_all <- cbind(theta_S_all, theta_S)
          
          pi_all <- cbind(pi_all, pi) 
          
        }
      }
    }
    
    num <- dim(R2_H_all)[2]
    
    
      if (is.null(num)==TRUE){stop(c("No pi vectors found that are in line with the restrictions based on the data."))}
  
        
    sum_pi_f <- colSums(pi_all[8:9,]) 
    Pi.Vectors <- 
      data.frame(t(rbind(pi_all, sum_pi_f))) 
    colnames(Pi.Vectors) <- c("Pi_0000", "Pi_0100", "Pi_0001", "Pi_0101",  
                              "Pi_1101", "Pi_1111", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f") 
        
  } # end mono S and T
  
  fit <- 
    list(Pi.Vectors= Pi.Vectors, R2_H=as.numeric(R2_H_all), Theta_T = as.numeric(theta_T_all), 
         Theta_S = as.numeric(theta_S_all), H_Delta_T = as.numeric(H_Delta_T_all), 
         Monotonicity=Monotonicity, Call=match.call())      
  
  class(fit) <- "ICA.BinBin"
  fit
  
} 