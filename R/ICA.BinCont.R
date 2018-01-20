
#M = 100; Dataset = CIGTS_BinCont
#Surr = IOP_12; True = IOP_96_Bin; Treat=Treat; Diff.Sigma=FALSE

#library(Surrogate); data(Schizo_BinCont); attach(Schizo_BinCont)
#M = 1; Dataset = Schizo_BinCont; attach(Schizo_BinCont)
#Surr = PANSS; True = CGI_Bin; Treat=Treat; Diff.Sigma=F; Min.Pval.S_0 = 0.01 
#Min.Pval.S_1 = 0.01; Diff.Sigma=FALSE; G_pi_00 = seq(0, 1, by=.01); 
#G_rho_01_00=seq(-1, 1, by=.01); G_rho_01_01=seq(-1, 1, by=.01); 
#G_rho_01_10=seq(-1, 1, by=.01); G_rho_01_11=seq(-1, 1, by=.01); 
#M=1000; Seed=123
#Save.Plots="/Users/WimVDE/Google Drive/Werk/_Recent/Ariel - BinCont/Appendix/Plots Mixtures/"
#Test.Fit.Details=TRUE; Keep.All=FALSE


ICA.BinCont <- function(Dataset, Surr, True, Treat, Diff.Sigma=FALSE, 
         G_pi_00 = seq(0, 1, by=.01), 
         G_rho_01_00=seq(-1, 1, by=.01), G_rho_01_01=seq(-1, 1, by=.01), 
         G_rho_01_10=seq(-1, 1, by=.01), G_rho_01_11=seq(-1, 1, by=.01), 
         M=1000, Seed=123, 
         Plots=TRUE, 
         Save.Plots="No", Test.Fit.Mixture=FALSE, Test.Fit.Mixture.Alpha=0.01, 
         Test.Fit.Details=FALSE, Keep.All=FALSE){          
  
  if (Test.Fit.Details==TRUE & Test.Fit.Mixture==FALSE){
    stop("Model fit details were requested (Test.Fit.Details=TRUE), but this requires that Test.Fit.Mixture=TRUE")
  }
  
  Nboot <- 1000
  
  set.seed(Seed)
  
  totaal <- M*length(G_pi_00)*length(G_rho_01_00)*length(G_rho_01_10)*length(G_rho_01_01)*length(G_rho_01_11)

  Surr <- (Dataset[,paste(substitute(Surr))])
  True <- (Dataset[,paste(substitute(True))])
  Treat <- Dataset[,paste(substitute(Treat))]
  
if (min(na.exclude(Treat))!=c(-1)) {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
if (na.exclude(max(Treat))!=c(1))  {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
if (length(unique(na.exclude(True)))>2) {stop("\nThe true endpoint should be binary.")}
if (min(na.exclude(True))!=c(0)) {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
if (max(na.exclude(True))!=c(1))  {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  
data_no_miss <- data.frame(na.exclude(cbind(Surr, True, Treat)))
data_conttreat <- subset(data_no_miss, Treat=="-1")      
data_exptreat <- subset(data_no_miss, Treat=="1")      
S_0 <- data_conttreat$Surr
S_1 <- data_exptreat$Surr

Surr <- na.exclude(Surr)   #LS
True <- na.exclude(True)   #LS
Treat <- na.exclude(Treat) #LS

pi_punt1 <- mean(data_exptreat$True) # E(T|Z=1)
pi_punt0 <- 1 - pi_punt1 # 1-E(T|Z=1)
pi_1punt <- mean(data_conttreat$True) # E(T|Z=0)
pi_0punt <- 1 - pi_1punt # 1 - E(T|Z=0)

count <- 0
R2_H_all <- PD_OK_all <- pi_00_all <- pi_10_all <- pi_01_all <- pi_11_all <- NULL
G_rho_01_00_all <- G_rho_01_01_all <- G_rho_01_10_all <- G_rho_01_11_all <- NULL
pi_Delta_T_min1_all <- pi_Delta_T_0_all <- pi_Delta_T_1_all <- NULL
pi_0_00_e_all <- pi_0_01_e_all <- pi_0_10_e_all <- pi_0_11_e_all <- NULL 
mu_0_00_all <- mu_0_01_all <- mu_0_10_all <- mu_0_11_all <- NULL
sigma_00_00_all <- sigma_00_01_all <- sigma_00_10_all <- sigma_00_11_all <- NULL
pi_1_00_e_all <- pi_1_01_e_all <- pi_1_10_e_all <- pi_1_11_e_all <- NULL 
mu_1_00_all <- mu_1_01_all <- mu_1_10_all <- mu_1_11_all <- NULL
sigma_11_00_all <- sigma_11_01_all <- sigma_11_10_all <- sigma_11_11_all <- NULL
aantal <- 0
voor_seed <- Seed
Fit.Mixture_S_0_OK <- Fit.Mixture_S_1_OK <- Fit.Mixture_S_0_OK_all <- Fit.Mixture_S_1_OK_all <- NULL
Test.Fit.Details_all <- NULL

while (aantal < M){
#flush.console(); 
#cat("\nAantal:", aantal)
  
#set.seed(i+Seed)
count <- count+1
voor_seed <- voor_seed + 1

G_rho_01_00_hier <- sample(G_rho_01_00, size = 1)
G_rho_01_01_hier <- sample(G_rho_01_01, size = 1)
G_rho_01_10_hier <- sample(G_rho_01_10, size = 1)
G_rho_01_11_hier <- sample(G_rho_01_11, size = 1)

set.seed(count)
pi_00_hier <- sample(G_pi_00, size=1)   
pi_10_hier <- pi_punt0 - pi_00_hier #UD
pi_01_hier <- pi_0punt -  pi_00_hier #UD
pi_11_hier <- pi_1punt -  pi_10_hier #UD 

#((pi_00_hier >= 0 & pi_10_hier >= 0 & pi_01_hier >= 0 & pi_11_hier >= 0) & 
#   (pi_00_hier <= 1 & pi_10_hier <= 1 & pi_01_hier <= 1 & pi_11_hier <= 1))

if ((pi_00_hier >= 0 & pi_10_hier >= 0 & pi_01_hier >= 0 & pi_11_hier >= 0) & 
  (pi_00_hier <= 1 & pi_10_hier <= 1 & pi_01_hier <= 1 & pi_11_hier <= 1)){

#flush.console()
#cat("\n\nPi's hier:", cbind(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), "\n")
  
pi_Delta_T_min1 <- pi_10_hier
pi_Delta_T_0 <- pi_00_hier + pi_11_hier
pi_Delta_T_1 <- pi_01_hier 

mix1 <- mix2 <- mix.object <- Sum_M_S_0 <- Sum_M_S_1 <- M1 <- M2 <- M3 <- M4 <- NULL
set.seed(voor_seed)
try(mix1 <- invisible(mixtools::normalmixEM(arbvar = Diff.Sigma, x = na.exclude(S_0), #maxrestarts = 100000,
                              lambda=c(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), k=4, 
                              maxit = 10000000)), silent=TRUE) 
set.seed(voor_seed*43)
try(mix2 <- invisible(mixtools::normalmixEM(arbvar = Diff.Sigma, x = na.exclude(S_1), 
                              lambda=c(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), k=4, 
                              maxit = 10000000)), silent=TRUE)

if (exists("mix1")==TRUE & exists("mix2")==TRUE & is.null(mix1)==FALSE & is.null(mix2)==FALSE){
  
# S 0
mix.object <- mix1
k <- ncol(mix.object$posterior)
x_S_0 <- x <- sort(mix.object$x)
a <- hist(x, plot = FALSE) 
maxy <- max(max(a$density), .3989*mix.object$lambda/mix.object$sigma)
if (length(mix.object$mu) == k && length(mix.object$sigma) == 1) {
 arbmean <- TRUE
 arbvar <- FALSE
}
if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {  
 arbmean <- TRUE
 arbvar <- TRUE
}
M1 <- mix.object$lambda[1] * dnorm(x, mean = mix.object$mu[1 * arbmean + (1 - arbmean)], 
                                     sd = mix.object$sigma[1 * arbvar + (1 - arbvar)])
M2 <- mix.object$lambda[2] * dnorm(x, mean = mix.object$mu[2 * arbmean + (1 - arbmean)], 
                                     sd = mix.object$sigma[2 * arbvar + (1 - arbvar)])
M3 <- mix.object$lambda[3] * dnorm(x, mean = mix.object$mu[3 * arbmean + (1 - arbmean)], 
                                     sd = mix.object$sigma[3 * arbvar + (1 - arbvar)])
M4 <- mix.object$lambda[4] * dnorm(x, mean = mix.object$mu[4 * arbmean + (1 - arbmean)], 
                                     sd = mix.object$sigma[4 * arbvar + (1 - arbvar)])
Sum_M_S_0 <- M1 + M2 + M3 + M4


if (Test.Fit.Mixture == TRUE){
  
  # The functions  below were programmed by Rand Wilcox
  # see http://dornsife.usc.edu/labs/rwilcox/software/
  qcomhdMC<-function(x,y,q=c(.25,.5,.75),
                     nboot=2000,plotit=FALSE,SEED=TRUE,xlab="Group 1",ylab="Est.1-Est.2",
                     alpha=.01,ADJ.CI=TRUE){
    
    #
    # Compare quantiles using pb2gen
    # via hd estimator. Tied values are allowed.
    #
    #  ADJ.CI=TRUE means that the confidence intervals are adjusted based on the level used by the corresponding 
    #  test statistic. If a test is performed with at the .05/3 level, for example, the confidence returned has 
    #  1-.05/3 probability coverage.
    # 
    # When comparing lower or upper quartiles, both power and the probability of Type I error
    # compare well to other methods that have been derived.
    # q: can be used to specify the quantiles to be compared
    # q defaults to comparing the .1,.25,.5,.75, and .9 quantiles
    #
    #   Function returns p-values and critical p-values based on Hochberg's method.
    #
    
    if(SEED)set.seed(2)
    pv=NULL
    output=matrix(NA,nrow=length(q),ncol=10)
    dimnames(output)<-list(NULL,c("q","n1","n2","est.1","est.2","est.1_minus_est.2","ci.low","ci.up","p_crit","p-value"))
    for(i in 1:length(q)){
      output[i,1]=q[i]
      output[i,2]=length(elimna(x))
      output[i,3]=length(elimna(y))
      output[i,4]=hd(x,q=q[i])
      output[i,5]=hd(y,q=q[i])
      output[i,6]=output[i,4]-output[i,5]
      temp=qcom.sub(x,y,nboot=nboot,q=q[i],SEED=FALSE,alpha=alpha)
      output[i,7]=temp$ci[1]
      output[i,8]=temp$ci[2]
      output[i,10]=temp$p.value                                                                                                      
    }                                                                                                                              
    temp=order(output[,10],decreasing=TRUE)                                                                                        
    zvec=alpha/c(1:length(q))                                                                                                      
    output[temp,9]=zvec     
    if(ADJ.CI){
      for(i in 1:length(q)){
        temp=pb2gen(x,y,nboot=nboot,est=hd,q=q[i],SEED=FALSE,alpha=output[i,9],pr=FALSE)
        output[i,7]=temp$ci[1]
        output[i,8]=temp$ci[2]
        output[i,10]=temp$p.value
      }
      temp=order(output[,10],decreasing=TRUE) 
    }                                                                                                       
    output <- data.frame(output)                                                                                                   
    output$signif=rep("YES",nrow(output))
    for(i in 1:nrow(output)){
      if(output[temp[i],10]>output[temp[i],9])output$signif[temp[i]]="NO"                                                            
      if(output[temp[i],10]<=output[temp[i],9])break                                                                                 
    }                                                                                                                              
    if(plotit){                                                                                                                    
      xax=rep(output[,4],3)                                                                                                          
      yax=c(output[,6],output[,7],output[,8])                                                                                        
      plot(xax,yax,xlab=xlab,ylab=ylab,type="n")                                                                                     
      points(output[,4],output[,6],pch="*")                                                                                          
      lines(output[,4],output[,6])                                                                                                   
      points(output[,4],output[,7],pch="+")                                                                                          
      points(output[,4],output[,8],pch="+")                                                                                          
    }                                                                                                                              
    output}                                                                                                                         
  
    pool.a.list<-function(x){
      #
      # x has list mode. Pool all of the data into a single R variable.
      #
      if(!is.list(x))stop("x should have list mode")
      pts=NULL
      for(j in 1:length(x))pts=c(pts,x[[j]])
      pts
    }
    
    onestep<-function(x,bend=1.28,na.rm=FALSE,MED=TRUE){
      #
      #  Compute one-step M-estimator of location using Huber's Psi.
      #  The default bending constant is 1.28
      #  
      #  MED=TRUE: initial estimate is the median
      #  Otherwise use modified one-step M-estimator
      #
      if(na.rm)x<-x[!is.na(x)]
      if(MED)init.loc=median(x)
      if(!MED)init.loc=mom(x,bend=bend)
      y<-(x-init.loc)/mad(x)  #mad in splus is madn in the book.
      A<-sum(hpsi(y,bend))
      B<-length(x[abs(y)<=bend])
      onestep<-median(x)+mad(x)*A/B
      onestep
    }
    
    hpsi<-function(x,bend=1.28){
      #
      #   Evaluate Huber`s Psi function for each value in the vector x
      #   The bending constant defaults to 1.28.
      #
      hpsi<-ifelse(abs(x)<=bend,x,bend*sign(x))
      hpsi
    }
    mom<-function(x,bend=2.24,na.rm=TRUE){
      #
      #  Compute MOM-estimator of location.
      #  The default bending constant is 2.24
      #
      if(na.rm)x<-x[!is.na(x)] #Remove missing values
      flag1<-(x>median(x)+bend*mad(x))
      flag2<-(x<median(x)-bend*mad(x))
      flag<-rep(T,length(x))
      flag[flag1]<-F
      flag[flag2]<-F
      mom<-mean(x[flag])
      mom
    }
    elimna<-function(m){
      #
      # remove any rows of data having missing values
      #
      DONE=FALSE
      if(is.list(m) && is.matrix(m)){
        z=pool.a.list(m)
        m=matrix(z,ncol=ncol(m))
        DONE=TRUE
      }
      if(!DONE){
        if(is.list(m) && is.matrix(m[[1]])){
          for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
          e=m
          DONE=TRUE
        }}
      if(!DONE){
        if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
          for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
          e=m
          DONE=TRUE
        }}
      if(!DONE){
        #if(!is.list(m)){
        #if(is.null(dim(m)))
        m<-as.matrix(m)
        ikeep<-c(1:nrow(m))
        for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
        e<-m[ikeep[ikeep>=1],]
        #}
      }
      e
    }
    hd<-function(x,q=.5,na.rm=TRUE,STAND=NULL){
      #
      #  Compute the Harrell-Davis estimate of the qth quantile
      #
      #  The vector x contains the data,
      #  and the desired quantile is q
      #  The default value for q is .5.
      #
      if(na.rm)x=elimna(x)
      n<-length(x)
      m1<-(n+1)*q
      m2<-(n+1)*(1-q)
      vec<-seq(along=x)
      w<-pbeta(vec/n,m1,m2)-pbeta((vec-1)/n,m1,m2)  # W sub i values
      y<-sort(x)
      hd<-sum(w*y)
      hd
    }
    qcom.sub<-function(x,y,q,alpha=.05,nboot=2000,SEED=TRUE){
      #
      x<-x[!is.na(x)] # Remove any missing values in x
      y<-y[!is.na(y)] # Remove any missing values in y
      if(SEED)set.seed(2) # set seed of random number generator so that
      #             results can be duplicated.
      datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
      datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
      datax=listm(t(datax))
      datay=listm(t(datay))
      bvecx<-parallel::mclapply(datax,hd,q,mc.preschedule=TRUE)
      bvecy<-parallel::mclapply(datay,hd,q,mc.preschedule=TRUE)
      bvecx=as.vector(matl(bvecx))
      bvecy=as.vector(matl(bvecy))
      bvec<-sort(bvecx-bvecy)
      low<-round((alpha/2)*nboot)+1
      up<-nboot-low
      temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
      sig.level<-2*(min(temp,1-temp))
      se<-var(bvec)
      list(est.1=hd(x,q),est.2=hd(y,q),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
    }
    listm<-function(x){
      #
      # Store the data in a matrix or data frame in a new
      # R variable having list mode.
      # Col 1 will be stored in y[[1]], col 2 in y[[2]], and so on.
      #
      if(is.null(dim(x)))stop("The argument x must be a matrix or data frame")
      y<-list()
      for(j in 1:ncol(x))y[[j]]<-x[,j]
      y
    }
    matl<-function(x){
      #
      # take data in list mode and store it in a matrix
      #
      J=length(x)
      nval=NA
      for(j in 1:J)nval[j]=length(x[[j]])
      temp<-matrix(NA,ncol=J,nrow=max(nval))
      for(j in 1:J)temp[1:nval[j],j]<-x[[j]]
      temp
    }
    pb2gen<-function(x,y,alpha=.05,nboot=100,est=onestep,SEED=TRUE,pr=FALSE,...){
      #
      #   Compute a bootstrap confidence interval for the
      #   the difference between any two parameters corresponding to
      #   independent groups.
      #   By default, M-estimators are compared.
      #   Setting est=mean, for example, will result in a percentile
      #   bootstrap confidence interval for the difference between means.
      #   Setting est=onestep will compare M-estimators of location.
      #   The default number of bootstrap samples is nboot=2000
      #
      x<-x[!is.na(x)] # Remove any missing values in x
      y<-y[!is.na(y)] # Remove any missing values in y
      if(SEED)set.seed(2) # set seed of random number generator so that
      #             results can be duplicated.
      datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
      datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
      bvecx<-apply(datax,1,est,...)
      bvecy<-apply(datay,1,est,...)
      bvec<-sort(bvecx-bvecy)
      low<-round((alpha/2)*nboot)+1
      up<-nboot-low
      temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
      sig.level<-2*(min(temp,1-temp))
      se<-var(bvec)
      list(est.1=est(x,...),est.2=est(y,...),est.dif=est(x,...)-est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
    }
  ############################################################

result <- result_S_0 <- result_S_1 <- NULL
# Test fit mixture
options(warn=-1)
A <- NULL
set.seed(voor_seed)
A <- sample(x = x_S_0, replace = T, prob = Sum_M_S_0, size = 5000)
result <- result_S_0 <- qcomhdMC(x = S_0, y = A, nboot = Nboot, alpha = Test.Fit.Mixture.Alpha)
Fit.Mixture_S_0_OK <- all(result$signif=="No")
#if (Test.Fit.Mixture.Alpha==0){Fit.Mixture_S_0_OK = TRUE}
options(warn=0)
}
if (Test.Fit.Mixture == FALSE){Fit.Mixture_S_0_OK = TRUE}


# S 1
mix.object <- mix2
k <- ncol(mix.object$posterior)
x_S_1 <- x <- sort(mix.object$x)
a <- hist(x, plot = FALSE)
maxy <- max(max(a$density), .3989*mix.object$lambda/mix.object$sigma)
if (length(mix.object$mu) == k && length(mix.object$sigma) == 1) {
  arbmean <- TRUE
  arbvar <- FALSE
}
if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {  
  arbmean <- TRUE
  arbvar <- TRUE
}
M1 <- mix.object$lambda[1] * dnorm(x, mean = mix.object$mu[1 * arbmean + (1 - arbmean)], 
                                   sd = mix.object$sigma[1 * arbvar + (1 - arbvar)])
M2 <- mix.object$lambda[2] * dnorm(x, mean = mix.object$mu[2 * arbmean + (1 - arbmean)], 
                                   sd = mix.object$sigma[2 * arbvar + (1 - arbvar)])
M3 <- mix.object$lambda[3] * dnorm(x, mean = mix.object$mu[3 * arbmean + (1 - arbmean)], 
                                   sd = mix.object$sigma[3 * arbvar + (1 - arbvar)])
M4 <- mix.object$lambda[4] * dnorm(x, mean = mix.object$mu[4 * arbmean + (1 - arbmean)], 
                                   sd = mix.object$sigma[4 * arbvar + (1 - arbvar)])
Sum_M_S_1 <- M1 + M2 + M3 + M4


if (Test.Fit.Mixture == TRUE){
result <- result_S_1 <- NULL
options(warn=-1)
B <- NULL
set.seed(voor_seed*25)
B <- sample(x = x_S_1, replace = T, prob = Sum_M_S_1, size = 5000)
result <- result_S_1 <- qcomhdMC(S_1, B, nboot = Nboot, alpha = Test.Fit.Mixture.Alpha)
Fit.Mixture_S_1_OK <- all(result$signif=="No")
#if (Test.Fit.Mixture.Alpha==0){Fit.Mixture_S_0_OK = TRUE}
options(warn=0)
}
if (Test.Fit.Mixture == FALSE){Fit.Mixture_S_1_OK = TRUE}


if (Test.Fit.Details == TRUE){
  temp_S_0 <- cbind("S_0", result_S_0); temp_S_1 <- cbind("S_1", result_S_1)
  names(temp_S_0)[1] <- names(temp_S_1)[1] <- "Outcome"
  temp <- rbind(temp_S_0, temp_S_1)
  temp <- temp[-c(3,4)]
  temp <- cbind((aantal+1), temp)
  names(temp)[1] <- c("Run") 
  Test.Fit.Details_all <- rbind(Test.Fit.Details_all, temp)
  }



if (Plots == TRUE){
  
  if (Save.Plots!="No"){  
  
  options(warn=-1)
  dir.create(Save.Plots)
  options(warn=0)
  
  pdf(paste(Save.Plots, "/Run", (aantal+1), ".pdf", sep=""),width=6,height=4,paper='special') 
  par(mfrow=c(1,2))
  max_val <- max(max(density(S_0)$y), max(Sum_M_S_0))
  
  hist(S_0, freq = F, ylim=c(0, max_val), xlab=expression(S[0]), 
       main=bquote(paste(S[0], " (run "~.(aantal+1), ")")), col="grey")
  lines(x_S_0, Sum_M_S_0, col="red", lwd=2, lty=2)
  
  #mixtools::plot.mixEM(mix1, whichplots = 2)
  #lines(x_S_0, Sum_M_S_0, lwd=3, col="grey", lty=3)
    
  #plot(density(S_0), ylim=c(0, max_val))
  #legend("topright", legend=c("observed", "Mixture"), col=c("black", "red"), lty=c(1,1), cex=.5)
    
  #mixtools::plot.mixEM(mix2, whichplots = 2)
  #lines(x_S_1, Sum_M_S_1, lwd=3, col="grey", lty=3)
    
  max_val <- max(max(density(S_0)$y), max(Sum_M_S_1))
  hist(S_1, freq = F, ylim=c(0, max_val), xlab=expression(S[1]), 
       main=bquote(paste(S[1], " (run "~.(aantal+1), ")")), col="grey")
  lines(x_S_1, Sum_M_S_1, col="red", lwd=2, lty=2)
  dev.off()
    }
 
    # Show plots on screen   
  par(mfrow=c(1,2))
  max_val <- max(max(density(S_0)$y), max(Sum_M_S_0))
  
  hist(S_0, freq = F, ylim=c(0, max_val), xlab=expression(S[0]), 
       main=bquote(paste(S[0], " (run "~.(aantal+1), ")")), col="grey")
  lines(x_S_0, Sum_M_S_0, col="red", lwd=2, lty=2)
  
  #mixtools::plot.mixEM(mix1, whichplots = 2)
  #lines(x_S_0, Sum_M_S_0, lwd=3, col="grey", lty=3)
  
  #plot(density(S_0), ylim=c(0, max_val))
  #legend("topright", legend=c("observed", "Mixture"), col=c("black", "red"), lty=c(1,1), cex=.5)
  
  #mixtools::plot.mixEM(mix2, whichplots = 2)
  #lines(x_S_1, Sum_M_S_1, lwd=3, col="grey", lty=3)
  
  max_val <- max(max(density(S_0)$y), max(Sum_M_S_1))
  hist(S_1, freq = F, ylim=c(0, max_val), xlab=expression(S[1]), 
       main=bquote(paste(S[1], " (run "~.(aantal+1), ")")), col="grey")
  lines(x_S_1, Sum_M_S_1, col="red", lwd=2, lty=2)
    par(mfrow=c(1,1))
}

# mixture components f(S_0) for use later
pi_0_00_e <- mix1$lambda[1]
pi_0_01_e <- mix1$lambda[2] 
pi_0_10_e <- mix1$lambda[3]
pi_0_11_e <- mix1$lambda[4]
mu_0_00 <- mix1$mu[1] 
mu_0_01 <- mix1$mu[2] 
mu_0_10 <- mix1$mu[3] 
mu_0_11 <- mix1$mu[4] 
sigma_00_00 <- mix1$sigma[1]**2  # notatie sigma_00_00 refereert naar var ipv SD (mix$sigma is SD)
sigma_00_01 <- mix1$sigma[2]**2 
sigma_00_10 <- mix1$sigma[3]**2 
sigma_00_11 <- mix1$sigma[4]**2 

# mixture components f(S_1)
pi_1_00_e <- mix2$lambda[1]
pi_1_01_e <- mix2$lambda[2] 
pi_1_10_e <- mix2$lambda[3]
pi_1_11_e <- mix2$lambda[4]
mu_1_00 <- mix2$mu[1] 
mu_1_01 <- mix2$mu[2] 
mu_1_10 <- mix2$mu[3] 
mu_1_11 <- mix2$mu[4] 
sigma_11_00 <- mix2$sigma[1]**2 
sigma_11_01 <- mix2$sigma[2]**2 
sigma_11_10 <- mix2$sigma[3]**2 
sigma_11_11 <- mix2$sigma[4]**2

# Check PD covariance matrices 
mat_a <- matrix(c(sigma_00_00, (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))), 
                  (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))), 
                  sigma_11_00), nrow = 2, byrow = TRUE)
eigen_mat_a <- min(eigen((mat_a))$values)

mat_b <- matrix(c(sigma_00_01, (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))), 
                  (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))), sigma_11_01), nrow = 2, byrow = TRUE)
eigen_mat_b <- min(eigen((mat_b))$values)

mat_c <- matrix(c(sigma_00_10, (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))), 
                  (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))), sigma_11_10), nrow = 2, byrow = TRUE)
eigen_mat_c <- min(eigen((mat_c))$values)

mat_d <- matrix(c(sigma_00_11, (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))), 
                  (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))), sigma_11_11), nrow = 2, byrow = TRUE)
eigen_mat_d <- min(eigen((mat_d))$values)



if (((((Fit.Mixture_S_0_OK == TRUE) & (Fit.Mixture_S_1_OK == TRUE))) & 
  ((eigen_mat_a > 0) & (eigen_mat_b > 0) & (eigen_mat_c > 0) & (eigen_mat_d > 0))) |
  (((Keep.All==TRUE)) & 
     ((eigen_mat_a > 0) & (eigen_mat_b > 0) & (eigen_mat_c > 0) & (eigen_mat_d > 0))))
{

mu_s00 <- mu_1_00 - mu_0_00
mu_s01 <- mu_1_01 - mu_0_01
mu_s10 <- mu_1_10 - mu_0_10
mu_s11 <- mu_1_11 - mu_0_11
omega_1 <- pi_00_hier / (pi_00_hier + pi_11_hier); omega_2 <- 1 - omega_1

sigma_s00 <- sigma_00_00 + sigma_11_00 - (2 * (sqrt(sigma_00_00 * sigma_11_00) * G_rho_01_00_hier)) #var
sigma_s01 <- sigma_00_01 + sigma_11_01 - (2 * (sqrt(sigma_00_01 * sigma_11_01) * G_rho_01_01_hier))
sigma_s10 <- sigma_00_10 + sigma_11_10 - (2 * (sqrt(sigma_00_10 * sigma_11_10) * G_rho_01_10_hier))
sigma_s11 <- sigma_00_11 + sigma_11_11 - (2 * (sqrt(sigma_00_11 * sigma_11_11) * G_rho_01_11_hier))


f_Delta_S <- function(val){
  v <- NA
  v <- (pi_00_hier * dnorm(val, mu_s00, sd = sqrt(sigma_s00))) + 
       (pi_01_hier * dnorm(val, mu_s01, sd = sqrt(sigma_s01))) +
       (pi_10_hier * dnorm(val, mu_s10, sd = sqrt(sigma_s10))) +
       (pi_11_hier * dnorm(val, mu_s11, sd = sqrt(sigma_s11)))
  return(v)
}

# I_10
#Seed <- Seed+1; set.seed(Seed)
S1 <- rnorm(n=10000, mean = mu_s10, sd = sqrt(sigma_s10))
I_10 <- mean(log(dnorm(S1, mu_s10, sd = sqrt(sigma_s10)) / f_Delta_S(S1)))

# I_01
#Seed <- Seed+1; set.seed(Seed)
S2 <- rnorm(n=10000, mean = mu_s01, sd = sqrt(sigma_s01))
I_01 <- mean(log(dnorm(S2, mu_s01, sd = sqrt(sigma_s01)) / f_Delta_S(S2)))

# I_00
#Seed <- Seed+1; set.seed(Seed)
S3 <- rnorm(n=10000, mean = mu_s00, sd = sqrt(sigma_s00))
I_00 <- mean(log(
  ((omega_1 * dnorm(S3, mu_s00, sd = sqrt(sigma_s00))) + 
  ((omega_2 * dnorm(S3, mu_s11, sd = sqrt(sigma_s11))))) / 
    f_Delta_S(S3)))

# I_11
#Seed <- Seed+1; set.seed(Seed)
S4 <- rnorm(n=10000, mean = mu_s11, sd = sqrt(sigma_s11))
I_11 <-mean(log(
  ((omega_1 * dnorm(S4, mu_s00, sd = sqrt(sigma_s00))) + 
     ((omega_2 * dnorm(S4, mu_s11, sd = sqrt(sigma_s11))))) / 
    f_Delta_S(S4)))

# R^2_H
I_Delta_T_Delta_S <- (pi_00_hier * I_00) + (pi_01_hier * I_01) + (pi_10_hier * I_10) + (pi_11_hier * I_11)
H_Delta_S <- 
  - ((pi_Delta_T_min1 * log(pi_Delta_T_min1)) + 
     (pi_Delta_T_0 * log(pi_Delta_T_0)) +
     (pi_Delta_T_1 * log(pi_Delta_T_1)))
R2_H <- I_Delta_T_Delta_S / H_Delta_S

# aantal
aantal <- aantal + 1
cat("\n", (aantal/M)*100, "% done. \n", sep="") 

# save results for output

R2_H_all <- c(R2_H_all, R2_H)
G_rho_01_00_all <- cbind(G_rho_01_00_all, G_rho_01_00_hier)
G_rho_01_01_all <- cbind(G_rho_01_01_all, G_rho_01_01_hier)
G_rho_01_10_all <- cbind(G_rho_01_10_all, G_rho_01_10_hier)
G_rho_01_11_all <- cbind(G_rho_01_11_all, G_rho_01_11_hier)
pi_00_all <- cbind(pi_00_all, pi_00_hier); pi_10_all <- cbind(pi_10_all, pi_10_hier)
pi_11_all <- cbind(pi_11_all, pi_11_hier); pi_01_all <- cbind(pi_01_all, pi_01_hier)
pi_Delta_T_min1_all <- cbind(pi_Delta_T_min1_all, pi_Delta_T_min1)
pi_Delta_T_0_all <- cbind(pi_Delta_T_0_all, pi_Delta_T_0)
pi_Delta_T_1_all <- cbind(pi_Delta_T_1_all, pi_Delta_T_1)
pi_0_00_e_all <- cbind(pi_0_00_e_all, pi_0_00_e)
pi_0_01_e_all <- cbind(pi_0_01_e_all, pi_0_01_e)
pi_0_10_e_all <- cbind(pi_0_10_e_all, pi_0_10_e)
pi_0_11_e_all <- cbind(pi_0_11_e_all, pi_0_11_e)
mu_0_00_all <- cbind(mu_0_00_all, mu_0_00)
mu_0_01_all <- cbind(mu_0_01_all, mu_0_01)
mu_0_10_all <- cbind(mu_0_10_all, mu_0_10)
mu_0_11_all <- cbind(mu_0_11_all, mu_0_11)
sigma_00_00_all <- cbind(sigma_00_00_all, sigma_00_00)
sigma_00_01_all <- cbind(sigma_00_01_all, sigma_00_01)
sigma_00_10_all <- cbind(sigma_00_10_all, sigma_00_10)
sigma_00_11_all <- cbind(sigma_00_11_all, sigma_00_11)
pi_1_00_e_all <- cbind(pi_1_00_e_all, pi_1_00_e)
pi_1_01_e_all <- cbind(pi_1_01_e_all, pi_1_01_e)
pi_1_10_e_all <- cbind(pi_1_10_e_all, pi_1_10_e)
pi_1_11_e_all <- cbind(pi_1_11_e_all, pi_1_11_e)
mu_1_00_all <- cbind(mu_1_00_all, mu_1_00)
mu_1_01_all <- cbind(mu_1_01_all, mu_1_01)
mu_1_10_all <- cbind(mu_1_10_all, mu_1_10)
mu_1_11_all <- cbind(mu_1_11_all, mu_1_11)
sigma_11_00_all <- cbind(sigma_11_00_all, sigma_11_00)
sigma_11_01_all <- cbind(sigma_11_01_all, sigma_11_01)
sigma_11_10_all <- cbind(sigma_11_10_all, sigma_11_10)
sigma_11_11_all <- cbind(sigma_11_11_all, sigma_11_11)
Fit.Mixture_S_0_OK_all <- c(Fit.Mixture_S_0_OK_all, Fit.Mixture_S_0_OK)
Fit.Mixture_S_1_OK_all <- cbind(Fit.Mixture_S_1_OK_all, Fit.Mixture_S_1_OK)


} # einde if eigen_mat_a en eigen_mat_b pos def

}

} # einde for (a in 1: M) loop
}
       
fit <- 
  list(R2_H=R2_H_all, pi_00=as.numeric(pi_00_all), pi_01=as.numeric(pi_01_all), 
       pi_10=as.numeric(pi_10_all), pi_11=as.numeric(pi_11_all), 
       G_rho_01_00=as.numeric(G_rho_01_00_all), 
       G_rho_01_01=as.numeric(G_rho_01_01_all),
       G_rho_01_10=as.numeric(G_rho_01_10_all), 
       G_rho_01_11=as.numeric(G_rho_01_11_all),
       pi_Delta_T_min1=as.numeric(pi_Delta_T_min1_all), 
       pi_Delta_T_0=as.numeric(pi_Delta_T_0_all),
       pi_Delta_T_1=as.numeric(pi_Delta_T_1_all),
       pi_0_00=as.numeric(pi_0_00_e_all),
       pi_0_01=as.numeric(pi_0_01_e_all),
       pi_0_10=as.numeric(pi_0_10_e_all),
       pi_0_11=as.numeric(pi_0_11_e_all), 
       mu_0_00=as.numeric(mu_0_00_all),
       mu_0_01=as.numeric(mu_0_01_all),
       mu_0_10=as.numeric(mu_0_10_all),
       mu_0_11=as.numeric(mu_0_11_all),
       sigma2_00_00=as.numeric(sigma_00_00_all),
       sigma2_00_01=as.numeric(sigma_00_01_all),
       sigma2_00_10=as.numeric(sigma_00_10_all),
       sigma2_00_11=as.numeric(sigma_00_11_all),
       pi_1_00=as.numeric(pi_1_00_e_all),
       pi_1_01=as.numeric(pi_1_01_e_all),
       pi_1_10=as.numeric(pi_1_10_e_all),
       pi_1_11=as.numeric(pi_1_11_e_all), 
       mu_1_00=as.numeric(mu_1_00_all),
       mu_1_01=as.numeric(mu_1_01_all),
       mu_1_10=as.numeric(mu_1_10_all),
       mu_1_11=as.numeric(mu_1_11_all),
       sigma2_11_00=as.numeric(sigma_11_00_all),
       sigma2_11_01=as.numeric(sigma_11_01_all),
       sigma2_11_10=as.numeric(sigma_11_10_all),
       sigma2_11_11=as.numeric(sigma_11_11_all),
       Fit.Mixture_S_0_OK=Fit.Mixture_S_0_OK_all,
       Fit.Mixture_S_1_OK=Fit.Mixture_S_0_OK_all,
       Test.Fit.Details=Test.Fit.Details_all,
       Call=match.call())   

class(fit) <- "ICA.BinCont"
fit

}  

