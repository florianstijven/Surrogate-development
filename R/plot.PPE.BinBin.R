

plot.PPE.BinBin<-function(x,Type="Density",Param="PPE",Xlab.PE,main.PE,ylab="density",Par=par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5.1,4.1,4.1,2.1)),
                           Cex.Legend=1,Cex.Position="bottomright",lwd=3,linety=1,color=1,Breaks=0.05,xlimits=c(0,1), ...){


if (missing(Xlab.PE)) {Xlab.PE<-ifelse(Param=="PPE",expression("Probability of a Prediction Error"),ifelse(Param=="RPE",
                                                    expression("Reduction in Prediction Error"),expression("Individual Causal Association")))}
if (missing(main.PE)) {main.PE<-" "}
if (missing(color)) {col<-c(1)}
if (missing(lwd)) {lwd<-3}
if (missing(linety)) {lty<-c(1)}

par=Par
graphs=c(par("mfrow")[1],par("mfrow")[2])
graphs_per_page=graphs[1]*graphs[2]
goma=c(par("oma")[1],par("oma")[2],par("oma")[3],par("oma")[4])
gmar=c(par("mar")[1],par("mar")[2],par("mar")[3],par("mar")[4])
par(mfrow=graphs, oma=goma, mar=gmar)
graphnum=0


if (Type=="Density"){


graphnum=graphnum+1; 
if (graphnum>graphs_per_page){dev.new();graphnum=1;par(mfrow=graphs, oma=goma,mar=gmar);}

if (Param=="PPE"){
plot(density(x$PPE),lwd=lwd,col=color[1],
           xlab=Xlab.PE,main=main.PE, ylab=ylab,xlim=xlimits,lty=linety[1])}

if (Param=="RPE"){
plot(density(x$RPE),lwd=lwd,col=color[1],
           xlab=Xlab.PE,main=main.PE, ylab=ylab,xlim=xlimits,lty=linety[1])}

if (Param=="ICA"){
plot(density(x$R2_H),lwd=lwd,col=color[1],
           xlab=Xlab.PE,main=main.PE, ylab=ylab,xlim=xlimits,lty=linety[1])}


}

if (Type=="Freq"){


graphnum=graphnum+1; 
if (graphnum>graphs_per_page){dev.new();graphnum=1;par(mfrow=graphs, oma=goma,mar=gmar);}

if (Param=="PPE"){
hist(x$PPE,xlab=Xlab.PE,ylab="Frequency",main=main.PE,xlim=xlimits, freq=TRUE,breaks=seq(0,1,Breaks),col=color)}

if (Param=="RPE"){
hist(x$RPE,xlab=Xlab.PE,ylab="Frequency",main=main.PE,xlim=xlimits, freq=TRUE,breaks=seq(0,1,Breaks),col=color)}

if (Param=="ICA"){
hist(x$R2_H,xlab=Xlab.PE,ylab="Frequency",main=main.PE,xlim=xlimits, freq=TRUE,breaks=seq(0,1,Breaks),col=color)}


}
}





