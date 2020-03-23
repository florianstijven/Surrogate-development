comb27.BinBin <-function(pi1_1_, pi1_0_, pi_1_1, pi_1_0, 
  pi0_1_, pi_0_1, Monotonicity=c("No"),M=1000, Seed=1){

#############################################################
# Run diagnostic on Monotonicity and Prediction             #
#############################################################


mon=length(Monotonicity)
for (i in 1:mon){if (!(Monotonicity[i]=="No" | Monotonicity[i]=="True.Endp" | Monotonicity[i]=="Surr.Endp" | Monotonicity[i]=="Surr.True.Endp"))
  { cat("Monotonicity='",Monotonicity[i],"'is not a valid option",sep=" ")
}
}


if (mon != 1) {
  cat("\nNote. Only 1 monotonicity assumption can be used as input \n")
} 



     
#############################################################
# Initialize                                                #
#############################################################


R2_H_all<-Pi_f_all<-Grid<-I<-seed<-pi_f_all<- result<-POS<-H_Delta_T_all <-result<-monotonicity_all<-Pe_all<-Prediction_all<-NULL
pi0000<-pi0100<-pi0001<-pi0101<-pi1101<-pi1111<-pi0011<-pi0111<-pi1100<-pi0010<-pi1110<-pi0110<-pi1010<-pi1000<-pi1001<-pi1011<-NULL
Pe_lowerbound_all<-I_Delta_T_Delta_S_all<-H_Delta_S_all<-support_delta_T_all<-combo_all<-pi_all<-index_all<-NULL


pi_0_0=1-pi_1_0 - pi_1_1 - pi_0_1
pi0_0_=1-pi1_0_-pi1_1_-pi0_1_

     
for (monotonicity in Monotonicity) 
{

#############################################################
#1) Define A_r and A_f                                      #
#############################################################


    A_r <- matrix(data=c(1, 0, 0, 0, 0, 0, 0,
                         1, 0, 0, 0, 1, 0, 0,
                         1, 0, 0, 0, 0, 0, 1,
                         1, 0, 0, 1, 0, 0, 0,
                         1, 0, 1, 1, 0, 0, 0,
                         1, 1, 0, 1, 0, 0, 0,
                         1, 0, 0, 0, 0, 1, 1), ncol=7)

    A_f <- matrix(data=c(1, 0, 0, 1, 0, 1, 0,
                         1, 0, 1, 0, 1, 0, 0,
                         1, 0, 0, 0, 0, 1, 0,
                         1, 1, 0, 0, 1, 0, 0,
                         1, 0, 0, 0, 1, 1, 0,
                         1, 1, 0, 0, 0, 0, 0,
                         1, 0, 1, 0, 0, 0, 0,
                         1, 0, 1, 0, 0, 0, 1,
                         1, 1, 0, 0, 0, 0, 1), ncol=9)
A=cbind(A_r,A_f)

#############################################################
#  Apply Restrictions in sampling                           #
#############################################################

#restrictions based on observed and monotonicity assumptions
    
    min_pi_0111 <- min(pi0_1_, pi_1_1)
    min_pi_1100 <- min(pi1_0_, pi_1_0)
    min_pi_0010 <- min(pi0_1_, pi_0_0)
    min_pi_1110 <- min(pi1_1_, pi_1_0)
    min_pi_0110 <- min(pi0_1_, pi_1_0)
    min_pi_1010 <- min(pi1_1_, pi_0_0)
    min_pi_1000 <- min(pi1_0_, pi_0_0)
    min_pi_1001 <- min(pi1_0_, pi_0_1)
    min_pi_1011 <- min(pi1_1_, pi_0_1)
 

##################################################################################################################
# Based on monotonicity assumptions: delete columns of matrix A, determine A_f and override sampling minima      #
##################################################################################################################

if (monotonicity=="True.Endp") {A<-A[,-16:-13]; min_pi_1010=min_pi_1000=min_pi_1001=min_pi_1011=0}
if (monotonicity=="Surr.Endp") {A<-A[,-10:-13]; min_pi_0010=min_pi_1110=min_pi_0110=min_pi_1010=0}
if (monotonicity=="Surr.True.Endp" ) {A<-A[,-16:-10]; min_pi_1010=min_pi_1000=min_pi_1001=min_pi_1011=min_pi_0010=min_pi_1110=min_pi_0110=0}


A_f=A[,8:ncol(A)]
free=ncol(A_f)
A=cbind(A_r,A_f)
invAr=solve(A_r)

     
#############################################################
# Define vector b                                           #
#############################################################


vector_b <- matrix(data=c(1, sample(size = 1, pi1_1_), sample(size = 1, pi1_0_), sample(size = 1, pi_1_1), 
                        sample(size = 1, pi_1_0), sample(size = 1, pi0_1_), sample(size = 1, pi_0_1)), ncol=1) 

    
#############################################################
# Specify B matrix                                          #
#############################################################

 B<-matrix(data=c(0,0,0,0,1,0,0,0,0,  
                  0,0,0,0,0,1,0,0,0,  
                  0,0,0,0,0,0,0,1,0, 
                  0,0,0,0,0,0,0,0,1,
                  0,0,0,0,0,0,0,1,0,
                  0,0,0,0,1,0,0,0,0,
                  0,0,0,0,1,0,0,0,0,
                  0,0,0,0,0,1,0,0,0,
                  0,0,0,0,1,0,0,0,0,
                  0,1,0,0,0,0,0,0,0,
                  0,1,0,0,0,0,0,0,0,
                  0,0,1,0,0,0,0,0,0,
                  1,0,0,0,0,0,0,0,0,
                  0,0,0,1,0,0,0,0,0,
                  0,0,0,0,0,0,1,0,0,
                  0,0,0,1,0,0,0,0,0), ncol=16)


#############################################################
#  Determine size of the sampling space                     #
#############################################################

    grid<-0.01
      
    pi_1 <- seq(0, min_pi_0111, by=grid)
    pi_2 <- seq(0, min_pi_1100, by=grid) 
    pi_3 <- seq(0, min_pi_0010, by=grid) 
    pi_4 <- seq(0, min_pi_1110, by=grid) 
    pi_5 <- seq(0, min_pi_0110, by=grid) 
    pi_6 <- seq(0, min_pi_1010, by=grid) 
    pi_7 <- seq(0, min_pi_1000, by=grid) 
    pi_8 <- seq(0, min_pi_1001, by=grid) 
    pi_9 <- seq(0, min_pi_1011, by=grid) 
    
    

#############################################################
#  Apply Monte Carlo Simulation                             #
#############################################################

for (index in 1:M){
Seed=Seed+1
set.seed(Seed)
 
        pi_0111 <- runif(n=1, min = 0, max=min_pi_0111)
        pi_1100 <- runif(n=1, min = 0, max=min_pi_1100)
        pi_0010 <- runif(n=1, min = 0, max=min_pi_0010)
        pi_1110 <- runif(n=1, min = 0, max=min_pi_1110)
        pi_0110 <- runif(n=1, min = 0, max=min_pi_0110)
        pi_1010 <- runif(n=1, min = 0, max=min_pi_1010)
        pi_1000 <- runif(n=1, min = 0, max=min_pi_1000)
        pi_1001 <- runif(n=1, min = 0, max=min_pi_1001)
        pi_1011 <- runif(n=1, min = 0, max=min_pi_1011)
   
        pi_f <- matrix(data = c(pi_0111, pi_1100, pi_0010, pi_1110, pi_0110, pi_1010, pi_1000, pi_1001, pi_1011), ncol = 1)
        pi_f_reduced<-pi_f
        if (monotonicity=="True.Endp") { pi_f_reduced<-pi_f[-6:-9]}
        if (monotonicity=="Surr.Endp") {pi_f_reduced<-pi_f[-3:-6]}
        if (monotonicity=="Surr.True.Endp" ) {pi_f_reduced<-pi_f[-3:-9]}

       
        pi_r <- solve(A_r) %*% (vector_b - (A_f %*% pi_f_reduced))
        
        if ((sum(pi_r >= 0 & pi_r <= 1) == 7)==TRUE) {
          
          for (l in 1: length(pi_r)){
            if (pi_r[l] < 0) {pi_r[l] <- c(0)}
            if (pi_r[l] > 1) {pi_r[l] <- c(1)}
          }
          
          pi <- rbind(pi_r, pi_f)
          
    
#############################################################
# Specify Prediction function matrix                        #
#############################################################

counter=ifelse(monotonicity=="Surr.Endp" | monotonicity=="Surr.True.Endp",1,3)

for (i in 1:counter){
 for (j in 4:6){
  for (k in 7:9){
  reeks=rep(0,81)
   combo=c(i,j,k)
   combinatie=ifelse(counter==3,i*100+j*10+k,j*10+k)
   combinatie2=ifelse(counter==3,paste(i-2,"/",j-5,"/",k-8),paste(j-5,"/",k-8))
   if (combo[1]==1){reeks[1]=1; reeks[11]=1;reeks[21]=1;} 
   if (combo[1]==2){reeks[28]=1; reeks[38]=1;reeks[48]=1;} 
   if (combo[1]==3){reeks[55]=1; reeks[65]=1;reeks[75]=1;} 
   if (combo[2]==4){reeks[4]=1; reeks[14]=1;reeks[24]=1;} 
   if (combo[2]==5){reeks[31]=1; reeks[41]=1;reeks[51]=1;} 
   if (combo[2]==6){reeks[58]=1; reeks[68]=1;reeks[78]=1;} 
   if (combo[3]==7){reeks[7]=1; reeks[17]=1;reeks[27]=1;} 
   if (combo[3]==8){reeks[34]=1; reeks[44]=1;reeks[54]=1;} 
   if (combo[3]==9){reeks[61]=1; reeks[71]=1;reeks[81]=1;} 

P<-t(matrix(data=reeks, ncol=9))

    
         mat=B%*%pi
         mat1=mat[1];mat2=mat[2];mat3=mat[3];mat4=mat[4];mat5=mat[5];mat6=mat[6];mat7=mat[7];mat8=mat[8];mat9=mat[9];

 
          sum_S_min1 <- mat1+mat2+mat3
          sum_S_0 <- mat4+mat5+mat6
          sum_S_1 <- mat7+mat8+mat9
          
          sum_T_min1 <- mat1+mat4+mat7
          sum_T_0 <- mat2+mat5+mat8
          sum_T_1 <- mat3+mat6+mat9

         fout=P%*%B%*%pi
         fout1=fout[1];fout2=fout[2];fout3=fout[3];fout4=fout[4];fout5=fout[5];fout6=fout[6];fout7=fout[7];fout8=fout[8];fout9=fout[9];
          
          Pe<-fout2+fout3+fout4+fout6+fout7+fout8
         
                                 
          I_Delta_T_Delta_S <- 
            ifelse(mat1==0,0,mat1*log2(mat1/(sum_S_min1*sum_T_min1)))+ifelse(mat2==0,0,mat2*log2(mat2/(sum_S_min1*sum_T_0)))+
            ifelse(mat3==0,0,mat3*log2(mat3/(sum_S_min1*sum_T_1)))  + ifelse(mat4==0,0,mat4*log2(mat4/(sum_S_0*sum_T_min1)))+
            ifelse(mat5==0,0,mat5*log2(mat5/(sum_S_0*sum_T_0))) +     ifelse(mat6==0,0,mat6*log2(mat6/(sum_S_0*sum_T_1)))+
            ifelse(mat7==0,0,mat7*log2(mat7/(sum_S_1*sum_T_min1)))+   ifelse(mat8==0,0,mat8*log2(mat8/(sum_S_1*sum_T_0)))+
            ifelse(mat9==0,0,mat9*log2(mat9/(sum_S_1*sum_T_1)))

          
          H_Delta_T <-  
            -(ifelse(mat1+mat4+mat7==0,0,(mat1+mat4+mat7)*log2(mat1+mat4+mat7))+ 
              ifelse(mat2+mat5+mat8==0,0,(mat2+mat5+mat8)*log2(mat2+mat5+mat8))+
              ifelse(mat3+mat6+mat9==0,0,(mat3+mat6+mat9)*log2(mat3+mat6+mat9)))

        
          H_Delta_S <-  
            -(ifelse(mat1+mat2+mat3==0,0,(mat1+mat2+mat3)*log2(mat1+mat2+mat3))+ 
              ifelse(mat4+mat5+mat6==0,0,(mat4+mat5+mat6)*log2(mat4+mat5+mat6))+
              ifelse(mat7+mat8+mat9==0,0,(mat7+mat8+mat9)*log2(mat7+mat8+mat9)))
 
          R2_H <- I_Delta_T_Delta_S / min(H_Delta_T, H_Delta_S)
       
         
          I_Delta_T_Delta_S_all=rbind(I_Delta_T_Delta_S_all,I_Delta_T_Delta_S)
          H_Delta_T_all <-rbind(H_Delta_T_all, H_Delta_T)
          H_Delta_S_all <-rbind(H_Delta_S_all, H_Delta_S)
          R2_H_all <- rbind(R2_H_all, R2_H)
          monotonicity_all<-rbind(monotonicity_all,monotonicity) 
          Pe_all<-rbind(Pe_all,Pe)
          combo_all<-rbind(combo_all,combinatie2)
          pi_all <- cbind(pi_all, pi) 
          index_all<-cbind(index_all,index)

 
                    
   }}}
     
  }
 }
}
fit<-data.frame(index=as.numeric(index_all),Monotonicity=as.character(monotonicity_all),Pe=as.numeric(Pe_all,Pe),
                combo=as.character(combo_all),
                R2_H=as.numeric(R2_H_all), H_Delta_T=as.numeric(H_Delta_T_all),
                H_Delta_S=as.numeric(H_Delta_S_all),
                I_Delta_T_Delta_S=as.numeric(I_Delta_T_Delta_S_all), stringsAsFactors = TRUE)
class(fit)<-"comb27.BinBin"
fit

}


