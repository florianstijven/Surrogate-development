plot.MaxEntContCont <- function(x, Type="Freq", Xlab, col, Main, Entropy.By.ICA=FALSE, ...){
  Object <- x

  if (missing(Xlab)) {Xlab <- expression(rho[Delta])}
  if (missing(col)) {col <- c(8)}
  if (missing(Main)) {Main <- " "}

  if (Type=="Density"){
    plot(density(Object$ICA.Fit$ICA, na.rm = T), xlab=Xlab, ylab="Density", main=Main, lwd=2, col=col, ...)
    abline(v=Object$ICA.Max.Ent, lwd=2)
  }

  if (Type=="Freq"){

    h <- hist(Object$ICA.Fit$ICA, plot = FALSE)
    h$density <- h$counts/sum(h$counts)
    cumulMidPoint <- ecdf(x=Object$ICA.Fit$ICA)(h$mids)
    labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
    plot(h,freq=T, xlab=Xlab, ylab="Frequency", col=col, main=Main, ...)
    abline(v=Object$ICA.Max.Ent, lwd=2)
  }

  if (Entropy.By.ICA == TRUE){
    plot(x = Object$Table.ICA.Entropy$Entropy, y=Object$Table.ICA.Entropy$ICA, xlab="Entropy", ylab = expression(rho[Delta]), ...)
  }
}
