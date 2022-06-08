plot.ICA.BinCont <- function(x, Histogram.ICA=TRUE, Mixmean=TRUE, Mixvar=TRUE, Deviance=TRUE,
                             Type="Percent", Labels=FALSE, ...){

  Object <- x 
  
  
  if (Histogram.ICA==TRUE){
    Xlab <- expression(R[H]^2)
    Main = " "
    
    dev.new()
  
    if (Type=="Density"){    
    plot(density(x$R2_H, na.rm = T), xlab=Xlab, ylab="Density", main=Main, lwd=2, ...)
    }
    if (Type=="Freq"){
      h <- hist(Object$R2_H, plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab, ylab="Frequency", main=Main, col="grey", ...)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab, ylab="Frequency", main=Main, labels=labs, col="grey", ...)
      }
    }
    
    if (Type=="Percent"){
      
      h <- hist(Object$R2_H, plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h, freq=F, xlab=Xlab, ylab="Percentage", main=Main, col="grey", ...)
      }
      if (Labels==TRUE){
        plot(h, freq=F, xlab=Xlab, ylab="Percentage", main=Main, labels=labs, col="grey", ...)
      }
      }
      
    if (Type=="CumPerc"){
      
      h <- hist(Object$R2_H, breaks=length(Object$R2_H), plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab, ylab="Cumulative percentage", col=0, main=Main, ...)
      lines(x=h$mids, y=cumulative)
    }
  } # end of histogram
  
  
if (Mixmean==TRUE){
  
Xlab <- expression(Run)
Ylab <- expression(Mean)
Main_S0 <- expression(S[0])
Main_S1 <- expression(S[1])
  
dev.new()
  par(mfrow=c(1,2))
  plot(Object$mean_Y_S0, type = "p", ylim = c(round(min(Object$mean_Y_S0, Object$mean.S0))-1, 
                                              round(max(Object$mean_Y_S0, Object$mean.S0))+1), 
       main = Main_S0, xlab = Xlab, ylab = Ylab, ...)
  abline(h = Object$mean.S0, col = "red", lty = 2, lwd = 4)
  
  plot(Object$mean_Y_S1, type = "p", ylim = c(round(min(Object$mean_Y_S1, Object$mean.S1))-1, round(max(Object$mean_Y_S1,Object$mean.S1))+1), 
       main = Main_S1, xlab = Xlab, ylab = Ylab, ...)
  abline(h = Object$mean.S1, col = "red", lty = 2, lwd = 4)
  par(mfrow=c(1,1))
  
}
  

  if (Mixvar==TRUE){
    
    Xlab <- expression(Run)
    Ylab <- expression(Variance)
    Main_S0 <- expression(S[0])
    Main_S1 <- expression(S[1])
    
    dev.new()
    par(mfrow=c(1,2))
    plot(Object$var_Y_S0, type = "p", ylim = c(round(min(Object$var_Y_S0, Object$var.S0))-1, 
                                                round(max(Object$var_Y_S0, Object$var.S0))+1), 
         main = Main_S0, xlab = Xlab, ylab = Ylab, ...)
    abline(h = Object$var.S0, col = "red", lty = 2, lwd = 4)
    
    plot(Object$var_Y_S1, type = "p", ylim = c(round(min(Object$var_Y_S1, Object$var.S1))-1, round(max(Object$var_Y_S1,Object$var.S1))+1), 
         main = Main_S1, xlab = Xlab, ylab = Ylab, ...)
    abline(h = Object$var.S1, col = "red", lty = 2, lwd = 4)
    par(mfrow=c(1,1))
    
  }
  
  
  
  if (Deviance==TRUE){
  dev.new()
  par(mfrow=c(1,2))
  boxplot(Object$dev_S0, main = expression(S[0]), ylab = "deviance", ylim = NULL)
  boxplot(Object$dev_S1, main = expression(S[1]), ylab = "deviance", ylim = NULL)
  par(mfrow=c(1,2))
  }
  
  
  }

 