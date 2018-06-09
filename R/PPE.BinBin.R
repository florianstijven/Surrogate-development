#
# matrices Paul
# M=number of valid vectors
# Monotonicity always 'No'
#


PPE.BinBin <-function(pi1_1_, pi1_0_, pi_1_1, pi_1_0, pi0_1_, pi_0_1, M=1000, Seed=1)
{

     
#############################################################
# Initialize                                                #
#############################################################


   R2_H_all<-Pi_f_all<-Grid<-I<-seed<-pi_f_all<- result<-POS<-H_Delta_T_all <-result<-Pe_all<-monotonicity_all<-Prediction_all<-NULL
   pi0000<-pi0100<-pi0001<-pi0101<-pi1101<-pi1111<-pi0011<-pi0111<-pi1100<-pi0010<-pi1110<-pi0110<-pi1010<-pi1000<-pi1001<-pi1011<-RPE_all<-NULL
   I_Delta_T_Delta_S_all<-H_Delta_S_all<-support_delta_T_all<-combo_all<-pi_all<-index_all<-PPE_T_all<-PPE_all<-NULL


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
   
   invAr=solve(A_r)

 
#############################################################
#  Apply Monte Carlo Simulation                             #
#############################################################

   
set.seed(Seed)
   
index<-0
runs<-0
while (index < M){
  
runs<-runs+1


######################################################################################
# Define vector b : this is a sampling step that allows to account for uncertainty.  #
# Constraint added to have valid vectors                                             #
######################################################################################


   vector_b <- matrix(data=c(1, sample(size = 1, pi1_1_), sample(size = 1, pi1_0_), sample(size = 1, pi_1_1), 
                        sample(size = 1, pi_1_0), sample(size = 1, pi0_1_), sample(size = 1, pi_0_1)), ncol=1) 

   valid=((vector_b[4]+vector_b[6]+vector_b[7]<=1) & (vector_b[2]+vector_b[3]+vector_b[5]<=1))

   if (valid==TRUE){

   PI1_1_=vector_b[2]
   PI1_0_=vector_b[3]
   PI_1_1=vector_b[4]
   PI_1_0=vector_b[5]
   PI0_1_=vector_b[6]
   PI_0_1=vector_b[7]
#this is the corrected version
   PI_0_0=1-vector_b[4]-vector_b[5]-vector_b[7]
   PI0_0_=1-vector_b[2]-vector_b[3]-vector_b[6]



#############################################################
#  Apply Restrictions in sampling based on observed data    #
#############################################################


#restrictions based on observed 
    
   min_pi_0111 <- min(PI0_1_, PI_1_1)
   min_pi_1100 <- min(PI1_0_, PI_1_0)
   min_pi_0010 <- min(PI0_1_, PI_0_0)
   min_pi_1110 <- min(PI1_1_, PI_1_0)
   min_pi_0110 <- min(PI0_1_, PI_1_0)
   min_pi_1010 <- min(PI1_1_, PI_0_0)
   min_pi_1000 <- min(PI1_0_, PI_0_0)
   min_pi_1001 <- min(PI1_0_, PI_0_1)
   min_pi_1011 <- min(PI1_1_, PI_0_1)


 
   #randomly generate freely varying parameters
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
  
   pi_r <- invAr %*% (vector_b - (A_f %*% pi_f))
   
   monotonicity=c('No')
        
    if ((sum(pi_r >= 0 & pi_r <= 1) == 7)==TRUE) {
          
          for (l in 1: length(pi_r)){
            if (pi_r[l] < 0) {pi_r[l] <- c(0)}
            if (pi_r[l] > 1) {pi_r[l] <- c(1)}
          }
          
    pi <- rbind(pi_r, pi_f)
   
    
    
######################################################################
# Compute PP_T, PPE and RPE from the individual causal effects table.#
######################################################################
    
     mat=B%*%pi
     mat1=mat[1];mat2=mat[2];mat3=mat[3];mat4=mat[4];mat5=mat[5];mat6=mat[6];mat7=mat[7];mat8=mat[8];mat9=mat[9];
 
     sum_S_min1 <- mat1+mat2+mat3
     sum_S_0 <- mat4+mat5+mat6
     sum_S_1 <- mat7+mat8+mat9
          
     sum_T_min1 <- mat1+mat4+mat7
     sum_T_0 <- mat2+mat5+mat8
     sum_T_1 <- mat3+mat6+mat9

     PPE_T<-1-max(sum_T_min1,sum_T_0,sum_T_1)

     max_min1=max(mat1,mat2,mat3)
     max_0=max(mat4,mat5,mat6)
     max_1=max(mat7,mat8,mat9)
         
          
     PPE<-1-(max_min1+max_0+max_1)
     RPE<-1-PPE/PPE_T
       
              
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
     monotonicity_all <- rbind(monotonicity_all, monotonicity)
   
         
     pi_all <- cbind(pi_all, pi) 
     index<-index+1
     index_all<-cbind(index_all,index)
     PPE_T_all<-rbind(PPE_T_all,PPE_T)
     PPE_all<-rbind(PPE_all,PPE)
     RPE_all<-rbind(RPE_all,RPE)
   
     num <- dim(R2_H_all)[2]
    
     Pi.Vectors <- data.frame(t(pi_all)) 
      # colnames(Pi.Vectors) <- c("Pi_0000", "Pi_0100", "Pi_0001", "Pi_0101", "Pi_1101", "Pi_1111", "Pi_0011", 
      #                            "Pi_0111", "Pi_1100", "Pi_0010", "Pi_1110", "Pi_0110", "Pi_1010", "Pi_1000", "Pi_1001", "Pi_1011") 
    
    }#valid=true
   }#pi_r>0
   }#while index<M

print(runs)

#Pi.Vectors= Pi.Vectors terugzetten
fit<-data.frame(index=as.numeric(index_all),PPE=as.numeric(PPE_all),RPE=as.numeric(RPE_all), Monotonicity=as.character(monotonicity_all),
                PPE_T=as.numeric(PPE_T_all), R2_H=as.numeric(R2_H_all), H_Delta_T=as.numeric(H_Delta_T_all),
                H_Delta_S=as.numeric(H_Delta_S_all), I_Delta_T_Delta_S=as.numeric(I_Delta_T_Delta_S_all))
class(fit)<-"PPE.BinBin"
fit
} #end function


