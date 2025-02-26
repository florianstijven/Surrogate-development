#' @export
plot.SPF.BinCont <- function(x, Histogram.SPF=TRUE, Causal.necessity=TRUE, Best.pred=TRUE, Max.psi=TRUE, ...)
  {

  Object <- x

  if (Histogram.SPF==TRUE){
    dev.new()

    par(mfrow=c(3,3))
    hist(Object$r_min1_min1,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = -1|", Delta, "S ", epsilon, "(-", infinity,",", a, ")]")),
                              list(a = Object$a))),...)
    hist(Object$r_min1_0,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = -1|", Delta, "S ", epsilon, "(", a,",", b, ")]")),
                              list(a = Object$a, b = Object$b))),...)
    hist(Object$r_min1_1,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = -1|", Delta, "S ", epsilon, "(", b,",", infinity, ")]")),
                              list(b = Object$b))),...)
    hist(Object$r_0_min1,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = 0|", Delta, "S ", epsilon, "(-", infinity,",", a, ")]")),
                              list(a = Object$a))),...)
    hist(Object$r_0_0,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = 0|", Delta, "S ", epsilon, "(", a,",", b, ")]")),
                              list(a = Object$a, b = Object$b))),...)
    hist(Object$r_0_1,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = 0|", Delta, "S ", epsilon, "(", b,",", infinity, ")]")),
                              list(b = Object$b))),...)
    hist(Object$r_1_min1,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = 1|", Delta, "S ", epsilon, "(-", infinity,",", a, ")]")),
                              list(a = Object$a))),...)
    hist(Object$r_1_0,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = 1|", Delta, "S ", epsilon, "(", a,",", b, ")]")),
                              list(a = Object$a, b = Object$b))),...)
    hist(Object$r_1_1,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=eval(substitute(expression(paste("P[", Delta, "T = 1|", Delta, "S ", epsilon, "(", b,",", infinity, ")]")),
                              list(b = Object$b))),...)
    par(mfrow=c(1,1))

  } # end of histogram

 if (Causal.necessity==TRUE){
   dev.new()

   hist(Object$P_DT_0_DS_0,
        main = "",
        xlim = c(0,1), ylim = NULL,
        xlab=expression(paste("P(", Delta, "T = 0|", Delta, "S = 0)")),...)

  }

  if (Best.pred==TRUE){
    dev.new()

    par(mfrow=c(1,3))
    barplot(table(Object$best.pred.min1),
            main = "",
            xlab = eval(substitute(expression(paste(psi["ab"],"(-", infinity,",", a, ")]")),
                                   list(a = Object$a))),
            ylab = "Frequency",...)
    barplot(table(Object$best.pred.0),
            main = "",
            xlab = eval(substitute(expression(paste(psi["ab"],"(", a,",", b, ")]")),
                                   list(a = Object$a, b = Object$b))),
            ylab = "Frequency",...)
    barplot(table(Object$best.pred.1),
            main = "",
            xlab = eval(substitute(expression(paste(psi["ab"],"(", b,",", infinity, ")]")),
                                   list(b = Object$b))),
            ylab = "Frequency",...)
    par(mfrow=c(1,1))

  }

  if (Max.psi==TRUE){
    dev.new()

    hist(Object$P_DT_psi_DS_max,
         main = "",
         xlim = c(0,1), ylim = NULL,
         xlab=expression(paste("P(", Delta, "T = ", psi, "(", Delta, "S))")),...)

  }

}

