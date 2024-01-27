plot.Fano.BinBin<-function(x,Type="Density",Xlab.R2_HL,main.R2_HL,ylab="density",Par=par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5.1,4.1,4.1,2.1)),
                           Cex.Legend=1,Cex.Position="top", lwd=3,linety=c(5,6,7),color=c(8,9,3), ...){

Delta=as.numeric(unique(x$delta))
uncert=as.numeric(unique(x$uncertainty))


if (missing(Xlab.R2_HL)) {Xlab.R2_HL<-expression(R[HL]^2)}
if (missing(main.R2_HL)) {main.R2_HL<-" "}
vector_legend<-NULL


if (Type=="Scatter"){

#sort to make line plots work
y=as.data.frame(x[1:12])
order.pi10=order(y$delta,y$uncertainty,y$pi_10)
x=y[order.pi10,]

pi10=as.numeric(unique(x$samplepi10))
meanvector <- matrix(, ncol = length(Delta), nrow = length(pi10))


  for (i in 1:length(pi10)){
   for (j in 1:length(Delta)){
     meanvector[i,j]=mean(x$R2_HL[x$delta==Delta[j] & x$samplepi10==pi10[i]])
      }
   }



if (missing(color)) {col<-seq(8:8+length(Delta))}
if (length(color)!=length(Delta)) {
  cat("\nNote. color vector and Delta vector have unequal length \n")
}
for (j in 1:length(Delta)){
  x$kleur[x$delta==Delta[j]]=color[j]
 }

for (j in 1:length(Delta)){
      helplegend<-paste("delta=",Delta[j])
      vector_legend<-cbind(vector_legend,helplegend)
  }

  Xlab<-expression(pi[10])
  Ylab<-expression(R[HL]^2)


  plot(x$samplepi10,x$R2_HL,xlab=Xlab,main="", ylab=Ylab,type="n")
if (length(uncert)==1) {
  for (j in 1:length(Delta)){
   for (i in 1:length(uncert)){
    points(x$samplepi10[x$delta==Delta[j] & x$uncertainty==uncert[i]] , x$R2_HL[x$delta==Delta[j] & x$uncertainty==uncert[i]],pch=1,col=color[j])
   }
  }
 }


if (length(uncert)>1) {
  for (j in 1:length(Delta)){
   for (i in 1:length(uncert)){
  lines(x$samplepi10[x$delta==Delta[j] & x$uncertainty==uncert[i]] , x$R2_HL[x$delta==Delta[j] & x$uncertainty==uncert[i]], lwd=lwd,col=color[j])
   }
  lines(pi10 , meanvector[,j], lwd=lwd+5,col=color[j])
  }
 }

  legend(bty="n",Cex.Position, col=color, pch=20, cex = Cex.Legend,legend=vector_legend)
}


if (Type=="Density"){
max=Delta*0
par=Par
if (missing(color)) {col<-seq(8:8+length(Delta))}
if (missing(lwd)) {lwd<-2}
if (missing(linety)) {lty<-seq(1:1+length(Delta))}

if (length(color)!=length(Delta)) {
  cat("\nNote. color vector and Delta vector have unequal length \n")
}
if (length(linety)!=length(Delta)) {
  cat("\nNote. linetype vector and Delta vector have unequal length \n")
}

for (j in 1:length(Delta)){
      try(max[j] <- max(density(x$R2_HL[x$delta==Delta[j]], na.rm = T)$y), silent=TRUE)
      helplegend<-paste("delta=",Delta[j])
      vector_legend<-cbind(vector_legend,helplegend)
  }
max_val <- max(max)

plot(density(x$R2_HL[x$delta==Delta[1]]),lwd=lwd,col=color[1],xlab=Xlab.R2_HL,main=main.R2_HL, ylab=ylab,xlim=c(0,1),ylim=c(0,max_val),lty=linety[1])

if (length(Delta)>1){
for (i in 2:length(Delta)){
     lines(density(x$R2_HL[x$delta==Delta[i]]),lwd=lwd,col=color[i],lty=linety[i])
  }
 }

legend(bty="n",Cex.Position, lwd=lwd, col=color, lty=linety, cex = Cex.Legend,legend=vector_legend)
#legend(bty="n",x=0.1,y=10, lwd=lwd, col=color, lty=linety, cex = Cex.Legend,legend=vector_legend)

}

if (Type=="Freq"){
par=Par
graphs=c(par("mfrow")[1],par("mfrow")[2])
graphs_per_page=graphs[1]*graphs[2]
goma=c(par("oma")[1],par("oma")[2],par("oma")[3],par("oma")[4])
gmar=c(par("mar")[1],par("mar")[2],par("mar")[3],par("mar")[4])
par(mfrow=graphs, oma=goma, mar=gmar)

graphnum=0


for (j in 1:length(Delta)){
graphnum=graphnum+1;
if (graphnum>graphs_per_page){graphnum=1;par(mfrow=graphs, oma=goma,mar=gmar);}

hist(x$R2_H[x$delta==Delta[j]],xlab=Xlab.R2_HL,ylab="Frequency",main=main.R2_HL,
     xlim=c(0,1), freq=TRUE,breaks=seq(0,1,0.01))
legend(Cex.Position,  cex = Cex.Legend, legend=paste("Delta=",Delta[j]))
#legend(bty="n",x=0.1,y=500,  cex = Cex.Legend, legend=paste("Delta=",Delta[j]))

}
}
}







