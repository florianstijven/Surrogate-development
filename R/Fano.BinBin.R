
Fano.BinBin <-function(pi1_ ,  pi_1 , rangepi10=c(0,min(pi1_,1-pi_1)), fano_delta=c(0.1), M=100, Seed=1)
                           
{



#############################################################
# Initialize                                                #
#############################################################


seed<-pi_f_all<-H_Delta_T_all<-R2_HL_all <-delta_all<-uncertainty_all<-NULL
pi00<-pi01<-pi11<-pi10<-maxpi10<-minpi10<-samplepi10<-PPE_T_all<-NULL


     
#############################################################
#1) Define A_r and A_f                                      #
#############################################################


A_r <- matrix(data=c(1, 0, 0, 
                     1, 1, 1,
                     1, 0, 1), ncol=3)

A_f <- matrix(data=c(1, 1, 0), ncol=1)

A=cbind(A_r,A_f)
invAr=solve(A_r)

   

#############################################################
# Specify B matrix                                          #
#############################################################

 B<-matrix(data=c(0,1,0,  
                  0,1,0, 
                  0,0,1,
                  1,0,0), ncol=4)


 
#############################################################
#  Apply Monte Carlo Simulation                             #
#############################################################

if (rangepi10[1]>rangepi10[2]) {
cat("\nNote. range for pi10 is not properly defined \n")
}



##########################################################################
# Define vector b : this is a sampling step to account for uncertainty   #
##########################################################################

N=length(pi1_)

for (uncertainty in 1:N){

Seed=Seed+uncertainty
set.seed(Seed)
vector_b <- matrix(data=c(1, sample(size = 1, pi1_), sample(size = 1, pi_1)), ncol=1) 

for (i in 1:M){



#############################################################
#  Apply Restrictions in sampling                           #
#############################################################

     
    
#we only use this one : take the minimum of the theoretical minimum and the invoked maximum

    #max_pi_10 <- min(rangepi10[2],min(vector_b[2], 1-vector_b[3]))
     max_pi_10 <- min(rangepi10[2],min(pi1_ , 1-pi_1))

 
 
   if (M>1) {samplepi_10 <- rangepi10[1]+ (i-1)*(max_pi_10 - rangepi10[1])/(M-1) }
   if (M==1){samplepi_10 <- rangepi10[1]}
        
           
        pi_f <- matrix(data = c(samplepi_10), ncol = 1)   
       
        pi_r <- solve(A_r) %*% (vector_b - (A_f %*% pi_f))
    
        
        if ((sum(pi_r >= 0 & pi_r <= 1) == 3)==TRUE) {
          
          for (l in 1: length(pi_r)){
            if (pi_r[l] < 0) {pi_r[l] <- c(0)}
            if (pi_r[l] > 1) {pi_r[l] <- c(1)}
          }
          
          pi <- rbind(pi_r, pi_f)
          
    
         mat=B%*%pi
         mat1=mat[1];mat2=mat[2];mat3=mat[3];
       
         sum_T_min1 <- mat1
         sum_T_0 <- mat2
         sum_T_1 <- mat3
         PPE_T=1-max(mat1,mat2,mat3)


        for (delta in fano_delta) 
        {
          support_delta_T=sum(ifelse(sum_T_min1>0,1,0),ifelse(sum_T_0>0,1,0),ifelse(sum_T_1>0,1,0)) 
          f_delta=-delta*log2(delta) - (1-delta)*log2(1-delta) + delta*log2(support_delta_T-1)
      
          H_Delta_T <-  
            -(ifelse(mat1==0,0,(mat1)*log2(mat1))+ 
              ifelse(mat2==0,0,(mat2)*log2(mat2))+
              ifelse(mat3==0,0,(mat3)*log2(mat3)))
          
          R2_HL=max(0,1-f_delta/H_Delta_T)

          R2_HL_all <- rbind(R2_HL_all, R2_HL)
          delta_all<-rbind(delta_all,delta)
          H_Delta_T_all<-rbind(H_Delta_T_all, H_Delta_T)
          PPE_T_all<-rbind(PPE_T_all,PPE_T)
          maxpi10<-rbind(maxpi10,rangepi10[2])
          minpi10<-rbind(minpi10,rangepi10[1])
          samplepi10<-rbind(samplepi10,samplepi_10)
          uncertainty_all<-rbind(uncertainty_all,uncertainty)

          
          pi00<-rbind(pi00,pi[1])
          pi11<-rbind(pi11,pi[2])
          pi01<-rbind(pi01,pi[3])
          pi10<-rbind(pi10,pi[4])
                   
  }
 }
}
}
fit<-list(R2_HL=as.numeric(R2_HL_all), H_Delta_T=as.numeric(H_Delta_T_all), PPE_T=as.numeric(PPE_T_all),minpi10=as.numeric(minpi10),maxpi10=as.numeric(maxpi10),samplepi10=as.numeric(samplepi10),
      delta=as.numeric(delta_all),uncertainty=as.numeric(uncertainty_all), pi_00=as.numeric(pi00),pi_11=as.numeric(pi11),pi_01=as.numeric(pi01),pi_10=as.numeric(pi10))
class(fit)<-"Fano.BinBin"
fit
}




