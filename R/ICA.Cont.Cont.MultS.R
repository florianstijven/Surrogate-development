ICA.ContCont.MultS <- function(M = 500, N, Sigma,
  G = seq(from=-1, to=1, by = .00001),
  Seed=c(123),  Show.Progress=FALSE){

SDs <- sqrt(diag(Sigma)); mu <- rep(0, times=length(SDs))
results_here <- results <- all_delta_S_T <- all_delta_S_T_here <- Lower.Dig.Corrs.All <- NULL

found <- 0
set.seed(Seed)
for (i in 1: M) {

  Sigma_c <- Sigma_c <- cov2cor(Sigma)
  num_elem <- dim(Sigma_c)[1] ** 2
  Sigma_c_orig <- Sigma_c
  size <- row_num <- col_num_ind <- 3; total_size <- dim(Sigma_c)[1]
  here <- Sigma_c

  while (size <= total_size) {
    here <- Sigma_c[(1:size), (1:size)] #LS!
    here[is.na(here)] <- sample(x = G, size = length(here[is.na(here)]), replace = TRUE) #LS!
    here[upper.tri(here)] = t(here)[upper.tri(here)]  #LS!

      while (det(here) < 0){
      here <- Sigma_c[(1:size), (1:size)]
      here[is.na(here)] <- sample(x = G, size = length(here[is.na(here)]), replace = TRUE)
      here[upper.tri(here)] = t(here)[upper.tri(here)]
    }

    Sigma_c[1:row_num, 1:col_num_ind] <- here
    row_num <- row_num + 1;
    col_num_ind <- col_num_ind + 1
    size <- size + 1
    }

  Sigma_c <- here

  Min.Eigen.Sigma <- try(min(eigen(Sigma_c)$values), TRUE)    # lowest eigenvalue

  if ((inherits(x = Min.Eigen.Sigma, what = "try-error")==FALSE) & (Min.Eigen.Sigma >= 0.00000000001)) {
    found <- found + 1
    if (Show.Progress==TRUE){cat((i/M)*100, "% done... ", sep="")}
    corMat <- Sigma_c
    varVec <- SDs**2
    n = nrow(corMat)
    sdMat = diag(sqrt(varVec))
    rtn = sdMat %*% corMat %*% t(sdMat)

    # some functions to compute variances and covariances in Sigma_{Delta}
    var_diff <- function(cov_mat){
      cov_val <- cov_mat[1,1] + cov_mat[2,2] - (2 * cov_mat[1,2])
      fit <- c(cov_val); fit
    }

    cov_2_diffs <- function(cov_mat){
      cov_val <- (cov_mat[2,2] - cov_mat[1,2]) - (cov_mat[2,1] - cov_mat[1,1])
      fit <- c(cov_val)
      fit
    }



    # Delta T
    A <- matrix(var_diff(cov_mat = rtn[1:2, 1:2]), nrow = 1)

    # Delta S Delta T
    B <- NULL
    Aantal <- (dim(Sigma_c)[1]-2) / 2
    rtn_part <- rtn[c(3:dim(rtn)[1]), c(1,2)]
    for (z in 1:Aantal){
      cov_mat_here <-
        rtn_part[c((z*2)-1, z*2), c(1:2)]
      B <- rbind(B, cov_2_diffs(cov_mat_here))
    }

    # Delta S Delta S
    Dim <- dim(Sigma_c)[1]
    Sub_mat_var <- rtn[c(3:Dim), c(3:Dim)]

    C <- matrix(NA, nrow = Aantal, ncol = Aantal)
    for (l in 1:Aantal){    # l refers to columns
      for (k in 1:Aantal){  # k refers to rows
        Sub_mat_var_hier <-
          Sub_mat_var[c((k*2)-1, k*2), c((l*2)-1, l*2)]
        C[k,l] <-
          cov_2_diffs(cov_mat = Sub_mat_var_hier)
      }
    }

    Delta <- cbind(rbind(A, B), rbind(t(B), C))

    #ICA <- 1 - (det(Delta) / (det(A) * det(C)))

    ICA <- (t(B) %*% solve(C) %*% B) / A
    Adj.ICA <- 1 - (1 - ICA) * ((N - 1) / (N - Aantal - 1))

    # Save lower diagonal Sigma, contains correlations
    Lower.Dig.Corrs.Here <- Sigma_c[lower.tri(Sigma_c)]

    # Fit models
    results_here <-
      (c(ICA, Adj.ICA))

    results <- rbind(results, results_here)
    Lower.Dig.Corrs.All <- data.frame(rbind(Lower.Dig.Corrs.All, Lower.Dig.Corrs.Here))

  }
   Sigma_c <- Sigma_c_orig  #LS!!
  }

row.names(results) <- NULL
row.names(Lower.Dig.Corrs.All) <- NULL
results <- data.frame(results, row.names = NULL)
names(results) <- c("ICA", "Adj.ICA")

fit <-
  list(R2_H=(as.numeric(results$ICA)), Corr.R2_H=(as.numeric(results$Adj.ICA)),  #All_Delta_S_T = all_delta_S_T,
       Lower.Dig.Corrs.Sigma=Lower.Dig.Corrs.All, Call=match.call())

class(fit) <- "ICA.ContCont.MultS"
fit
}




# Summary

summary.ICA.ContCont.MultS <- function(object, ..., Object){
  if (missing(Object)){Object <- object}

  mode <- function(data) {
    x <- data
    z <- density(x)
    mode_val <- z$x[which.max(z$y)]
    fit <- list(mode_val= mode_val)
  }

  cat("\nFunction call:\n\n")
  print(Object$Call)

  cat("\n\n\n# Uncorrected R2_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat("Mean R2_H: ", format(round(mean(Object$R2_H), 4), nsmall = 4), " (", format(round(sd(Object$R2_H), 4), nsmall = 4), ")",
      "  [min: ", format(round(min(Object$R2_H), 4), nsmall = 4), "; max: ",  format(round(max(Object$R2_H), 4), nsmall = 4), "]", sep="")
  cat("\nMode R2_H: ", format(round(mode(Object$R2_H)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(Object$R2_H, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)

  cat("\n\n\n# Bias-corrected R2_H results summary")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
  cat("Mean adjusted R2_H: ", format(round(mean(Object$Corr.R2_H), 4), nsmall = 4), " (", format(round(sd(Object$Corr.R2_H), 4), nsmall = 4), ")",
      "  [min: ", format(round(max(min(Object$Corr.R2_H), 0), 4), nsmall = 4), "; max: ",  format(round(max(Object$Corr.R2_H), 4), nsmall = 4), "]", sep="")
  cat("\nMode adjusted R2_H: ", format(round(mode(Object$Corr.R2_H)$mode_val, 4), nsmall = 4))
  cat("\n\nQuantiles of the distribution: \n\n")
  quant <- quantile(Object$Corr.R2_H, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)

  }




# Plot

plot.ICA.ContCont.MultS <- function(x, R2_H=FALSE, Corr.R2_H=TRUE, Type="Percent", Labels=FALSE,
                                    Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), col,
                                    Prediction.Error.Reduction=FALSE, ...){

  Object <- x
  if (missing(col)) {col <- c(8)}
  options(warn=-100)

  if (R2_H==TRUE){
    par=Par
    if (Type=="Freq"){
      h <- hist(Object$R2_H, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs_vals <- (h$counts)/sum((h$counts))
      labs <- paste(round((labs_vals), digits=4)*100, "%", sep="")

      if (Labels==FALSE){
        plot(h,freq=T, xlab=expression(R[H]^2), ylab="Frequency", col=col, main=" ", ...)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=expression(R[H]^2), ylab="Frequency", col=col, main=" ", labels=labs, ...)
      }
    }

    if (Type=="Percent"){
      h <- hist(Object$R2_H, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs_vals <- (h$counts)/sum((h$counts))
      labs <- paste(round((labs_vals), digits=4)*100, "%", sep="")

      if (Labels==FALSE){
        plot(h,freq=F, xlab=expression(R[H]^2), ylab="Percentage", col=col, main=" ", ...)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=expression(R[H]^2), ylab="Percentage", col=col, main=" ", labels=labs, ...)
      }
    }

    if (Type=="CumPerc"){
      h <- hist(Object$R2_H, breaks=length(Object$R2_H), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=expression(R[H]^2), ylab="Cumulative percentage", col=0, main=" ", ...)
      lines(x=h$mids, y=cumulative)
    }
  }

  if (Corr.R2_H==TRUE){
    par=Par
    if (Type=="Freq"){
      h <- hist(Object$Corr.R2_H, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$Corr.R2_H)(h$mids)
      labs_vals <- (h$counts)/sum((h$counts))
      labs <- paste(round((labs_vals), digits=4)*100, "%", sep="")


      if (Labels==FALSE){
        plot(h,freq=T, xlab=expression(paste("Adj. ", R[H]^2)), ylab="Frequency", col=col, main=" ", ...)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=expression(paste("Adj. ", R[H]^2)), ylab="Frequency", col=col, main=" ", labels=labs, ...)
      }
    }

    if (Type=="Percent"){
      h <- hist(Object$Corr.R2_H, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$Corr.R2_H)(h$mids)
      labs_vals <- (h$counts)/sum((h$counts))
      labs <- paste(round((labs_vals), digits=4)*100, "%", sep="")

      if (Labels==FALSE){
        plot(h,freq=F, xlab=expression(paste("Corr. ", R)[H]^2), ylab="Percentage", col=col, main=" ", ...)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=expression(paste("Corr. ", R)[H]^2), ylab="Percentage", col=col, main=" ", labels=labs, ...)
      }
    }

    if (Type=="CumPerc"){
      h <- hist(Object$Corr.R2_H, breaks=length(Object$Corr.R2_H), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=expression(paste("Corr. ", R)[H]^2), ylab="Cumulative percentage", col=0, main=" ", ...)
      lines(x=h$mids, y=cumulative)
    }
  }

  if (Prediction.Error.Reduction==TRUE){

    if (exists("Object$Res_Err_Delta_T")==FALSE){
      stop("To obtain a plot which shows the prediction error reduction, use the function ICA.ContCont.MultS_alt \ninstead of ICA.ContCont.MultS\n ")
    }

    max_x <- max(max(density(Object$Res_Err_Delta_T)$x), max(density(Object$Res_Err_Delta_T_Given_S)$x))
    max_y <- max(max(density(Object$Res_Err_Delta_T)$y), max(density(Object$Res_Err_Delta_T_Given_S)$y))

    plot(density(Object$Res_Err_Delta_T), xlim=c(0, max_x*1.1), ylim=c(0, max_y*1.1),
         lwd=2, main="", xlab=expression(paste("Prediction error")))
    lines(density(Object$Res_Err_Delta_T_Given_S), col="black", lwd=2, lty=2)
    legend("topright", legend=expression(paste(Delta, "T"), paste(Delta, "T|", Delta, "S")),
           lwd=c(2, 2), lty=c(1, 2), col=c("black", "black"), ...)
  }

  options(warn=0)
  }

