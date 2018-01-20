plot.SPF.BinCont <- function(x, Type="Frequency", Col="grey", Main, 
          Xlab=TRUE, ...){
                            
  Object <- x 

  if (missing(Col)) {Col = "grey"}
  
  if (Type=="Most.Likely.DeltaT"){
    if (missing(Main)){Main <- bquote(paste(Delta, "S in interval ["~.(Object$a), ",", ~.(Object$b), "]"))}
    probs <- na.exclude(data.frame(cbind(Object$P_Delta_T_min1, Object$P_Delta_T_0, Object$P_Delta_T_1)))
    names(probs) <- c(-1, 0, 1)
    z<-apply(probs,1,which.max) 
    Most.Likely.Delta.T <- as.numeric(as.character(names(probs)[z])) 
  
    if (Xlab==TRUE){
    hist(Most.Likely.Delta.T, col=Col, xlim = c(-1.5, 1.5), breaks = seq(-1.5, 1.5), xaxt="no", freq = F,
         main=Main, ylab="Percentage", 
         xlab=expression(paste("Most likely ", Delta, T, " given interval ", Delta, S)))
    } else {
      hist(Most.Likely.Delta.T, col=Col, xlim = c(-1.5, 1.5), breaks = seq(-1.5, 1.5), xaxt="no", freq = F,
           main=Main, ylab="Percentage", 
           xlab="")  
    }
  
    axis(side=1, at=c(-1, 0, 1)); box()  
  }
  
  
  if (Type=="Frequency"){
    if (missing(Main)){Main <- bquote(paste("P(", Delta, "T|", Delta, "S in interval ["~.(Object$a), ",", ~.(Object$b), "])"))}
    par(mfrow=c(1,3), oma=c(0,0,2,0))
    # plot Delta T = -1
    T_hier <- Object$P_Delta_T_min1
    if (Xlab==TRUE){
    hist(T_hier, xlab=expression(paste("P(", Delta, "T = -1|", Delta, "S)")), main="", col=Col, xlim=c(0,1))
    } else {hist(T_hier, xlab="", main="", col=Col, xlim=c(0,1))}
  
    # plot Delta T = 0
    T_hier <- Object$P_Delta_T_0
    if (Xlab==TRUE){
    hist(T_hier, xlab=expression(paste("P(", Delta, "T = 0|", Delta, "S)")), main="", col=Col, xlim=c(0,1))
      } else {hist(T_hier, xlab="", main="", col=Col, xlim=c(0,1))}
    
    title(Main, outer = T)
    
    # plot Delta T = 1
    T_hier <- Object$P_Delta_T_1
    if (Xlab==TRUE){
      hist(T_hier, xlab=expression(paste("P(", Delta, "T = 1|", Delta, "S)")), main="", col=Col, xlim=c(0,1))
    } else {hist(T_hier, xlab="", main="", col=Col, xlim=c(0,1))
        
      }
    
    par(mfrow=c(1,1))
  }
  
  if (Type=="Percentage"){
    if (missing(Main)){Main <- bquote(paste("P(", Delta, "T|", Delta, "S in interval ["~.(Object$a), ",", ~.(Object$b), "])"))}
    par(mfrow=c(1,3), oma=c(0,0,2,0))
    
  # Max voor y-range
  T_hier <- Object$P_Delta_T_min1
  h <- hist(T_hier, plot = FALSE)
  h$density <- h$counts/sum(h$counts)
  max_1 <- max(h$density)
  T_hier <- Object$P_Delta_T_0
  h <- hist(T_hier, plot = FALSE)
  h$density <- h$counts/sum(h$counts)
  max_2 <- max(h$density)
  T_hier <- Object$P_Delta_T_1
  h <- hist(T_hier, plot = FALSE)
  h$density <- h$counts/sum(h$counts)
  max_3 <- max(h$density)
  maxy <- max(max_1, max_2, max_3)
  ylim_val <- maxy * 1.15
    
  # plot Delta T = -1
  T_hier <- Object$P_Delta_T_min1
  h <- hist(T_hier, plot = FALSE)
  h$density <- h$counts/sum(h$counts)
  cumulMidPoint <- ecdf(x=T_hier)(h$mids)
  if (Xlab==TRUE){
  plot(h, freq=F, xlab=expression(paste("P(", Delta, "T = -1|", Delta, "S)")), ylab="Percentage", main="", col=Col, xlim=c(0,1), ylim=c(0, ylim_val))
  } else {plot(h, freq=F, xlab="", ylab="Percentage", main="", col=Col, xlim=c(0,1), ylim=c(0, ylim_val))
}
  
  # plot Delta T = 0
  T_hier <- Object$P_Delta_T_0
  h <- hist(T_hier, plot = FALSE)
  h$density <- h$counts/sum(h$counts)
  cumulMidPoint <- ecdf(x=T_hier)(h$mids)
  if (Xlab==TRUE){
    plot(h, freq=F, xlab=expression(paste("P(", Delta, "T = 0|", Delta, "S)")), ylab="Percentage", main="", col=Col, xlim=c(0,1), ylim=c(0, ylim_val))
  } else {
    plot(h, freq=F, xlab="", ylab="Percentage", main="", col=Col, xlim=c(0,1), ylim=c(0, ylim_val))  
  }
    title(Main, outer = T)
  
  
    # plot Delta T = 1
  T_hier <- Object$P_Delta_T_1
  h <- hist(T_hier, plot = FALSE)
  h$density <- h$counts/sum(h$counts)
  cumulMidPoint <- ecdf(x=T_hier)(h$mids)
  if (Xlab==TRUE){
  plot(h, freq=F, xlab=expression(paste("P(", Delta, "T = 1|", Delta, "S)")), ylab="Percentage", main="", col=Col, xlim=c(0,1), ylim=c(0, ylim_val))
  } else {
    plot(h, freq=F, xlab=" ", ylab="Percentage", main="", col=Col, xlim=c(0,1), ylim=c(0, ylim_val))
  } 
    
  par(mfrow=c(1,1))
  }
  
  
  }


    