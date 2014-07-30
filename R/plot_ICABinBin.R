plot.ICA.BinBin <- function(x, C3=FALSE, R2_H=TRUE, R_H=FALSE, Theta_T=FALSE, Theta_S=FALSE, 
                            Type="Percent", Labels=FALSE, Xlab.C3, Main.C3, Xlab.R2_H, Main.R2_H, Xlab.R_H, Main.R_H,
                            Xlab.Theta_S, Main.Theta_S, Xlab.Theta_T, Main.Theta_T, Cex.Legend=1, All.Densities=FALSE,
                            col, Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), ...){
  Object <- x 
  
  if (C3==TRUE){
    dev.new()
    par=Par 
    if (missing(Xlab.C3)) {Xlab.C3 <- expression(C[3])}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.C3)) {Main.C3=expression(C[3])}  
    
    
    if (Type=="Density"){
      plot(density(x$C3), xlab=Xlab.C3, ylab="Density", main=Main.C3, lwd=2)
    }
    
    if (Type=="All.Densities"){
      resul <- cbind(x$Pi.Vectors, x$C3, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
      colnames(resul) <-  
        c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
          "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "C3", "R2_H", 
          "Theta_T", "Theta_S", "H_Delta_T")
      C3_General <- resul$C3
      C3_No <- resul$C3[resul$Monotonicity=="No"]
      C3_Surr <- resul$C3[resul$Monotonicity=="Surr"]
      C3_True <- resul$C3[resul$Monotonicity=="True"]
      C3_SurrTrue <- resul$C3[resul$Monotonicity=="SurrTrue"]
      max_val <- max(max(density(C3_General)$y), max(density(C3_No)$y),  max(density(C3_Surr)$y),  max(density(C3_True)$y),  max(density(C3_SurrTrue)$y))   
      plot(density(x$C3), xlab=Xlab.C3, ylab="Density", main=Main.C3, lwd=2, ylim = c(0, max_val), col=0)
      try(lines(density((C3_No)), lty=1, col=1, lwd=3), silent=TRUE) 
      try(lines(density((C3_Surr)), lty=2, col=2, lwd=3), silent=TRUE)
      try(lines(density((C3_True)), lty=3, col=3, lwd=3), silent=TRUE)
      try(lines(density((C3_SurrTrue)), lty=4, col=4, lwd=2), silent=TRUE)
      try(lines(density((C3_General)), lty=5, col=6, lwd=2), silent=TRUE)
      
      legend("topleft", lwd=c(3, 3, 3, 3, 3), col=c(1, 2, 3, 4, 6), lty=c(1, 2, 3, 4, 5), cex = Cex.Legend,
             legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T", "General analysis"))
    }
    
    
    if (Type=="Freq"){
      h <- hist(Object$C3, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$C3)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.C3, ylab="Frequency", col=col, main=Main.C3)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.C3, ylab="Frequency", col=col, main=Main.C3, labels=labs)
      }
    }
    
    if (Type=="Percent"){
      h <- hist(Object$C3, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$C3)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=F, xlab=Xlab.C3, ylab="Percentage", col=col, main=Main.C3)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=Xlab.C3, ylab="Percentage", col=col, main=Main.C3, labels=labs)
      }
    }
    
    if (Type=="CumPerc"){
      h <- hist(Object$C3, breaks=length(Object$C3), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab.C3, ylab="Cumulative percentage", col=0, main=Main.C3)
      lines(x=h$mids, y=cumulative)
    }    
  }
  
  if (R2_H==TRUE){
    dev.new()
    par=Par 
    if (missing(Xlab.R2_H)) {Xlab.R2_H <- expression(R[H]^2)}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.R2_H)) {Main.R2_H <- expression(R[H]^2)}  
    
    
    if (Type=="Density"){
      plot(density(x$R2_H), xlab=Xlab.R2_H, ylab="Density", main=Main.R2_H, lwd=2)
    }
    
    if (Type=="All.Densities"){
      resul <- cbind(x$Pi.Vectors, x$C3, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
      colnames(resul) <-  
        c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
          "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "C3", "R2_H", 
          "Theta_T", "Theta_S", "H_Delta_T")
      R2_H_General <- resul$R2_H
      R2_H_No <- resul$R2_H[resul$Monotonicity=="No"]
      R2_H_Surr <- resul$R2_H[resul$Monotonicity=="Surr"]
      R2_H_True <- resul$R2_H[resul$Monotonicity=="True"]
      R2_H_SurrTrue <- resul$R2_H[resul$Monotonicity=="SurrTrue"]
      max_val <- max(max(density(R2_H_General)$y), max(density(R2_H_No)$y),  max(density(R2_H_Surr)$y),  max(density(R2_H_True)$y),  max(density(R2_H_SurrTrue)$y))   
      plot(density(x$R2_H), xlab=Xlab.R2_H, ylab="Density", main=Main.R2_H, lwd=2, ylim = c(0, max_val), col=0)
      try(lines(density((R2_H_No)), lty=1, col=1, lwd=3), silent=TRUE) 
      try(lines(density((R2_H_Surr)), lty=2, col=2, lwd=3), silent=TRUE)
      try(lines(density((R2_H_True)), lty=3, col=3, lwd=3), silent=TRUE)
      try(lines(density((R2_H_SurrTrue)), lty=4, col=4, lwd=2), silent=TRUE)
      try(lines(density((R2_H_General)), lty=5, col=6, lwd=2), silent=TRUE)
      
      legend("topleft", lwd=c(3, 3, 3, 3, 3), col=c(1, 2, 3, 4, 6), lty=c(1, 2, 3, 4, 5), cex = Cex.Legend,
             legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T", "General analysis"))
    }
    
    
    
    if (Type=="Freq"){
      h <- hist(Object$R2_H, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main=Main.R2_H)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main=Main.R2_H, labels=labs)
      }
    }
    
    if (Type=="Percent"){
      h <- hist(Object$R2_H, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=F, xlab=Xlab.R2_H, ylab="Percentage", col=col, main=Main.R2_H)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=Xlab.R2_H, ylab="Percentage", col=col, main=Main.R2_H, labels=labs)
      }
    }
    
    if (Type=="CumPerc"){
      h <- hist(Object$R2_H, breaks=length(Object$R2_H), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab.R2_H, ylab="Cumulative percentage", col=0, main=Main.R2_H)
      lines(x=h$mids, y=cumulative)
    }    
  }
  
  
  if (R_H==TRUE){
    dev.new()
    par=Par 
    if (missing(Xlab.R_H)) {Xlab.R_H <- expression(R[H])}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.R_H)) {Main.R_H <- expression(R[H])}  
    
    
    if (Type=="Density"){
      plot(density(sqrt(x$R2_H)), xlab=Xlab.R_H, ylab="Density", main=Main.R_H, lwd=2)
    }
    
    if (Type=="All.Densities"){
      resul <- cbind(x$Pi.Vectors, x$C3, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
      colnames(resul) <-  
        c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
          "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "C3", "R2_H", 
          "Theta_T", "Theta_S", "H_Delta_T")
      R2_H_General <- resul$R2_H
      R2_H_No <- resul$R2_H[resul$Monotonicity=="No"]
      R2_H_Surr <- resul$R2_H[resul$Monotonicity=="Surr"]
      R2_H_True <- resul$R2_H[resul$Monotonicity=="True"]
      R2_H_SurrTrue <- resul$R2_H[resul$Monotonicity=="SurrTrue"]
      max_val <- max(max(density(sqrt(R2_H_General))$y), max(density(sqrt(R2_H_No))$y),  max(density(sqrt(R2_H_Surr))$y),
                     max(density(sqrt(R2_H_True))$y),  max(density(sqrt(R2_H_SurrTrue))$y))   
      
      plot(density(sqrt(x$R2_H)), xlab=Xlab.R_H, ylab="Density", main=Main.R_H, lwd=2, ylim = c(0, max_val), col=0)
      try(lines(density(sqrt(R2_H_No)), lty=1, col=1, lwd=3), silent=TRUE) 
      try(lines(density(sqrt(R2_H_Surr)), lty=2, col=2, lwd=3), silent=TRUE)
      try(lines(density(sqrt(R2_H_True)), lty=3, col=3, lwd=3), silent=TRUE)
      try(lines(density(sqrt(R2_H_SurrTrue)), lty=4, col=4, lwd=2), silent=TRUE)
      try(lines(density(sqrt(R2_H_General)), lty=5, col=6, lwd=2), silent=TRUE)
      
      legend("topleft", lwd=c(3, 3, 3, 3, 3), col=c(1, 2, 3, 4, 6), lty=c(1, 2, 3, 4, 5), cex = Cex.Legend,
             legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T", "General analysis"))
    }
    
    
    if (Type=="Freq"){
      h <- hist(sqrt(Object$R2_H), ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=sqrt(Object$R2_H))(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main=Main.R_H)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main=Main.R_H, labels=labs)
      }
    }
    
    if (Type=="Percent"){
      h <- hist(sqrt(Object$R2_H), ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=sqrt(Object$R2_H))(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=F, xlab=Xlab.R_H, ylab="Percentage", col=col, main=Main.R_H)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=Xlab.R_H, ylab="Percentage", col=col, main=Main.R_H, labels=labs)
      }
    }
    
    if (Type=="CumPerc"){
      h <- hist(sqrt(Object$R2_H), breaks=length(sqrt(Object$R2_H)), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab.R_H, ylab="Cumulative percentage", col=0, main=Main.R_H)
      lines(x=h$mids, y=cumulative)
    }    
  }
  
  
  if (Theta_S==TRUE){
    dev.new()
    par=Par 
    if (missing(Xlab.Theta_S)) {Xlab.Theta_S <- expression(theta[S])}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.Theta_S)) {Main.Theta_S  <- expression(theta[S])}  
    
    if (Type=="Freq"){
      h <- hist(Object$Theta_S, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$Theta_S)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main=Main.Theta_S)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main=Main.Theta_S, labels=labs)
      }
    }
    
    if (Type=="Percent"){
      h <- hist(Object$Theta_S, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$Theta_S)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Percentage", col=col, main=Main.Theta_S)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Percentage", col=col, main=Main.Theta_S, labels=labs)
      }
    }
    
    if (Type=="CumPerc"){
      h <- hist(Object$Theta_S, breaks=length(Object$Theta_S), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_S, ylab="Cumulative percentage", col=0, main=Main.Theta_S)
      lines(x=h$mids, y=cumulative)
    }    
  }
  
  if (Theta_T==TRUE){
    dev.new()
    par=Par 
    if (missing(Xlab.Theta_T)) {Xlab.Theta_T <- expression(theta[T])}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.Theta_T)) {Main.Theta_T  <- expression(theta[T])}  
    
    if (Type=="Freq"){
      h <- hist(Object$Theta_T, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$Theta_T)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main=Main.Theta_T)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main=Main.Theta_T, labels=labs)
      }
    }
    
    if (Type=="Percent"){
      h <- hist(Object$Theta_T, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$Theta_T)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Percentage", col=col, main=Main.Theta_T)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Percentage", col=col, main=Main.Theta_T, labels=labs)
      }
    }
    
    if (Type=="CumPerc"){
      h <- hist(Object$Theta_T, breaks=length(Object$Theta_T), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_T, ylab="Cumulative percentage", col=0, main=Main.Theta_T)
      lines(x=h$mids, y=cumulative)
    }    
  }
  
}
