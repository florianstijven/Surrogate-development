#' @export
plot.comb27.BinBin <- function(x,lab, ...){
y=as.data.frame(x[1:5])
a1=with(y,aggregate(y$Pe,list(index=y$index,Monotonicity=y$Monotonicity),min))
a2=merge(a1, y, by = c("index","Monotonicity"))
a3=as.data.frame(subset(a2,a2$Pe==a2$x))
a2=a2[,c(-3)]

#save counts
counts<-as.data.frame(table(a3$combo))
nonzerofreq=subset(counts,counts$Freq>0)
order.freq=order(-nonzerofreq$Freq)
z=nonzerofreq[order.freq,]

Mono=as.character(unique(x$Monotonicity))



  pos=barplot(z$Freq, main="Prediction Functions",
              xlab="Frequency", names=z$Var1,ylim=c(0,z$Freq[1]*1.4))
  barplot(z$Freq, main=paste("Prediction Functions",lab),
              xlab="", las=2,names=z$Var1,ylim=c(0,z$Freq[1]*1.4))

  text(x=pos,y=z$Freq+round(z$Freq[1]/10),labels=as.character(z$Freq))
}
