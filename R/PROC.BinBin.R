

PROC.BinBin <-function(Dataset=Dataset, Surr=Surr, True=True, Treat=Treat, BS=FALSE, seqs=250, MC_samples=1000, Seed=1){

  
  Surr <- Dataset[,paste(substitute(Surr))]
  True <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  
  Dataset <- na.exclude(data.frame(cbind(Surr, True, Treat)))
  if (length(unique(Dataset$Treat))!=2) stop("Please make sure that the treatment variable has only 2 levels.")
  if ((sort(unique(Dataset$Treat))[1]==c(-1) & sort(unique(Dataset$Treat))[2]==c(1))==FALSE)
    stop("Please make sure that the treatment is coded as control = -1 and experimental = 1.")
  
  if (length(unique(Dataset$Surr)) > 2) {stop("\nWarning: S should be coded with only two values!\n\n")}
  if ((min(unique(Dataset$Surr)) != 0) && (length(unique(Dataset$Surr)) == 2)){stop("\nWarning: S should be coded as 0/1!\n\n")}
  if ((max(unique(Dataset$Surr)) != 1) && (length(unique(Dataset$Surr)) == 2)){stop("\nWarning: S should be coded as 0/1!\n\n")}
  
  if (length(unique(Dataset$True)) > 2) {stop("\nWarning: T should be coded with only two values!\n\n")}
  if ((min(unique(Dataset$True)) != 0) && (length(unique(Dataset$True)) == 2)){stop("\nWarning: T should be coded as 0/1!\n\n")}
  if ((max(unique(Dataset$True)) != 1) && (length(unique(Dataset$True)) == 2)){stop("\nWarning: T should be coded as 0/1!\n\n")}

    
rpe_all<-ppe_all<-R2_H_all<-ppe_T_all<-NULL



for (i in 1:seqs){
#bootstrap
if (BS==TRUE) {
A1<-subset(Dataset,Dataset$Treat==-1)
B1<-subset(Dataset,Dataset$Treat== 1)
A2=A1[sample(nrow(A1),replace=T),  ]
B2=B1[sample(nrow(B1),replace=T),  ]
bsdb=as.data.frame(rbind(A2,B2))                        
}
  
if (BS==FALSE) {
bsdb=as.data.frame(Dataset)                        
  }
#determine marginal probabilities
pi1_1_ <- (dim(subset(x=bsdb, subset=(bsdb$True==1 & bsdb$Surr==1 & bsdb$Treat==-1)))[1])/dim(bsdb[bsdb$Treat==-1,])[1] 
pi0_1_ <- (dim(subset(x=bsdb, subset=(bsdb$True==0 & bsdb$Surr==1 & bsdb$Treat==-1)))[1])/dim(bsdb[bsdb$Treat==-1,])[1]
pi1_0_ <- (dim(subset(x=bsdb, subset=(bsdb$True==1 & bsdb$Surr==0 & bsdb$Treat==-1)))[1])/dim(bsdb[bsdb$Treat==-1,])[1]
pi0_0_ <- (dim(subset(x=bsdb, subset=(bsdb$True==0 & bsdb$Surr==0 & bsdb$Treat==-1)))[1])/dim(bsdb[bsdb$Treat==-1,])[1] 

pi_1_1 <- (dim(subset(x=bsdb, subset=(bsdb$True==1 & bsdb$Surr==1 & bsdb$Treat==1)))[1])/dim(bsdb[bsdb$Treat==1,])[1] 
pi_1_0 <- (dim(subset(x=bsdb, subset=(bsdb$True==1 & bsdb$Surr==0 & bsdb$Treat==1)))[1])/dim(bsdb[bsdb$Treat==1,])[1]
pi_0_1 <- (dim(subset(x=bsdb, subset=(bsdb$True==0 & bsdb$Surr==1 & bsdb$Treat==1)))[1])/dim(bsdb[bsdb$Treat==1,])[1]
pi_0_0 <- (dim(subset(x=bsdb, subset=(bsdb$True==0 & bsdb$Surr==0 & bsdb$Treat==1)))[1])/dim(bsdb[bsdb$Treat==1,])[1]

Seed=Seed+i
set.seed(Seed)
#run of PPE.BinBin
ppe <- PPE.BinBin(pi1_1_=pi1_1_, pi0_1_=pi0_1_, pi1_0_=pi1_0_,
                    pi_1_1=pi_1_1, pi_1_0=pi_1_0, pi_0_1=pi_0_1, 
                    M=MC_samples,Seed=Seed) 

  
  
ppe_T_all=rbind(ppe_T_all,ppe$PPE_T)  
rpe_all=rbind(rpe_all,ppe$RPE)
ppe_all=rbind(ppe_all,ppe$PPE)
R2_H_all=rbind(R2_H_all,ppe$R2_H)
}

fit=list(RPE=as.numeric(rpe_all),PPE=as.numeric(ppe_all),R2_H=as.numeric(R2_H_all),PPE_T=as.numeric(ppe_T_all))
class(fit)<-"PPE.BinBin"
fit
}









