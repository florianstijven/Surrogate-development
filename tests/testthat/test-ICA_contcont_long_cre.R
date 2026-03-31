library(testthat)
library(pracma)
library(fastmatrix)

ICA_contcont_long_cre <- function(p,
                                  time_points,
                                  R_structure = c("AR1", "CS"),
                                  R_rho = numeric(),
                                  VT0S0,
                                  VT1S1,
                                  VT0T0,
                                  VT1T1,
                                  VS0S0,
                                  VS1S1,
                                  VT0T1 = seq(-1, 1, by = 0.1),
                                  VT0S1 = seq(-1, 1, by = 0.1),
                                  VT1S0 = seq(-1, 1, by = 0.1),
                                  VS0S1 = seq(-1, 1, by = 0.1),
                                  D,
                                  N,
                                  M = 10000,
                                  prob = NULL,
                                  Seed = 123,
                                  Show.Progress = FALSE) {
  build_covariance_matrix <- function(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1,
                                      T0T0, T1T1, S0S0, S1S1) {
    Sigma_c <- diag(c(T0T0, T1T1, S0S0, S1S1))

    Sigma_c[2, 1] <- Sigma_c[1, 2] <- T0T1 * sqrt(T0T0 * T1T1)
    Sigma_c[3, 1] <- Sigma_c[1, 3] <- T0S0
    Sigma_c[4, 1] <- Sigma_c[1, 4] <- T0S1 * sqrt(T0T0 * S1S1)
    Sigma_c[3, 2] <- Sigma_c[2, 3] <- T1S0 * sqrt(T1T1 * S0S0)
    Sigma_c[4, 2] <- Sigma_c[2, 4] <- T1S1
    Sigma_c[4, 3] <- Sigma_c[3, 4] <- S0S1 * sqrt(S0S0 * S1S1)

    Sigma_c
  }

  is_positive_definite <- function(mat) {
    chol_result <- try(chol(mat), silent = TRUE)
    !inherits(chol_result, "try-error")
  }

  build_residual_correlation <- function(p, structure, rho) {
    if (structure == "AR1") {
      return(fastmatrix::corAR1(rho, p))
    }

    corCS(rho, p = p)
  }

  p <- as.numeric(p)
  R_structure <- match.arg(R_structure)
  R_rho <- as.numeric(R_rho)
  time_points <- as.numeric(time_points)
  p <- as.integer(p)

  if (length(time_points) != p) {
    stop("`time_points` must have length `p`.")
  }

  cor_matrix <- build_residual_correlation(p = p, structure = R_structure, rho = R_rho)
  v_grid <- expand.grid(
    T0T1 = VT0T1,
    T0S0 = VT0S0,
    T0S1 = VT0S1,
    T1S0 = VT1S0,
    T1S1 = VT1S1,
    S0S1 = VS0S1,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  complete_V <- list()
  if (nrow(v_grid) > 0) {
    for (i in seq_len(nrow(v_grid))) {
      V_i <- build_covariance_matrix(
        T0T1 = v_grid$T0T1[i],
        T0S0 = v_grid$T0S0[i],
        T0S1 = v_grid$T0S1[i],
        T1S0 = v_grid$T1S0[i],
        T1S1 = v_grid$T1S1[i],
        S0S1 = v_grid$S0S1[i],
        T0T0 = VT0T0,
        T1T1 = VT1T1,
        S0S0 = VS0S0,
        S1S1 = VS1S1
      )

      if (is_positive_definite(V_i)) {
        complete_V[[length(complete_V) + 1]] <- V_i
      }
    }
  }

  D_corr <- cov2cor(D)
  lower_template <- D_corr[lower.tri(D_corr)]
  missing_lower <- which(is.na(lower_template))
  D_sd <- sqrt(diag(D))
  complete_D <- list()

  if (length(missing_lower) == 0) {
    D_matrix <- D
    D_matrix[upper.tri(D_matrix)] <- t(D_matrix)[upper.tri(D_matrix)]
    if (is_positive_definite(D_matrix)) {
      complete_D[[1]] <- D_matrix
    }
  } else {
    ica_mpc <- suppressWarnings(
      ICA.ContCont.MultS.MPC(
        M = M,
        N = N,
        Sigma = D,
        prob = prob,
        Seed = Seed,
        Save.Corr = TRUE,
        Show.Progress = Show.Progress
      )
    )

    valid_rows <- which(complete.cases(ica_mpc$Lower.Dig.Corrs.All[, missing_lower, drop = FALSE]))

    for (i in valid_rows) {
      D_corr_i <- D_corr
      lower_values <- D_corr_i[lower.tri(D_corr_i)]
      lower_values[missing_lower] <- as.numeric(ica_mpc$Lower.Dig.Corrs.All[i, missing_lower])
      D_corr_i[lower.tri(D_corr_i)] <- lower_values
      D_corr_i[upper.tri(D_corr_i)] <- t(D_corr_i)[upper.tri(D_corr_i)]
      diag(D_corr_i) <- 1

      D_matrix <- D_corr_i * outer(D_sd, D_sd)
      D_matrix[upper.tri(D_matrix)] <- t(D_matrix)[upper.tri(D_matrix)]

      if (is_positive_definite(D_matrix)) {
        complete_D[[length(complete_D) + 1]] <- D_matrix
      }
    }
  }

  num_pos_def_v <- length(complete_V)
  num_pos_def_d <- length(complete_D)
  r2_lambda <- numeric(0)
  Z0 <- cbind(1, time_points)
  Q <- matrix(c(-1, 0, 1, 0,
                0, -1, 0, 1), nrow = 2, ncol = 4)
  Z_full <- kronecker(diag(4), Z0)
  A <- kronecker(Q, diag(p)) %*% Z_full

  if (num_pos_def_v > 0 && num_pos_def_d > 0) {
    for (i in seq_len(num_pos_def_v)) {
      sigma_u <- Q %*% complete_V[[i]] %*% t(Q)

      for (j in seq_len(num_pos_def_d)) {
        sigma_delta <- A %*% complete_D[[j]] %*% t(A) + kronecker(sigma_u, cor_matrix)
        sigma_delta <- round(sigma_delta, digits = 6)
        sigma_tt <- sigma_delta[1:p, 1:p, drop = FALSE]
        sigma_ss <- sigma_delta[(p + 1):(2 * p), (p + 1):(2 * p), drop = FALSE]

        if (is_positive_definite(sigma_delta) && is_positive_definite(sigma_u)) {
          pair_value <- 1 - det(sigma_delta) / (det(sigma_tt) * det(sigma_ss))
          r2_lambda <- c(r2_lambda, pair_value)
        }
      }
    }
  }

  fit <- list(
    Num.Pos.Def.V = num_pos_def_v,
    Num.Pos.Def.D = num_pos_def_d,
    Num.Pos.Def.Pairs = length(r2_lambda),
    R2_Lambda = r2_lambda,
    Call = match.call()
  )

  class(fit) <- "ICA_contcont_long_cre"
  fit
}


D <- matrix(0, nrow = 8, ncol = 8)
un11 <- D[1, 1] <- 0.04162
un22 <- D[3, 3] <- 0.04140
un31 <- D[5, 1] <- D[1, 5] <- 0.02306
un33 <- D[5, 5] <- 0.03991
un42 <- D[7, 3] <- D[3, 7] <- 0.02896
un44 <- D[7, 7] <- 0.04811
un51 <- D[1, 2] <- D[2, 1] <- 0.001000
un53 <- D[2, 5] <- D[5, 2] <- 0.001224
un55 <- D[2, 2] <- 0.001086
un62 <- D[3, 4] <- D[4, 3] <- 0.000439
un64 <- D[4, 7] <- D[7, 4] <- 0.000025
un66 <- D[4, 4] <- 0.001452
un71 <- D[1, 6] <- D[6, 1] <- 0.001114
un73 <- D[5, 6] <- D[6, 5] <- 0.000374
un75 <- D[2, 6] <- D[6, 2] <- 0.000684
un77 <- D[6, 6] <- 0.000919
un82 <- D[3, 8] <- D[8, 3] <- 0.000918
un84 <- D[7, 8] <- D[8, 7] <- -0.00053
un86 <- D[4, 8] <- D[8, 4] <- 0.001472
un88 <- D[8, 8] <- 0.002100


test_that("ICA_contcont_long_cre",
          {p=6
          time_points = c(0, 1, 2, 4, 6, 8)
          R_rho = 0.000236
          R_structure = "CS"
          VT0S0 = 0.009203
          VT1S1 = 0.008759
          VT0T0 = 0.01365
          VT1T1 = 0.01379
          VS0S0 = 0.01819
          VS1S1 = 0.01854
          VT0T1 = 0
          VT0S1 = 0
          VT1S0 = 0
          VS0S1 = 0
          D = D
          N = 568
          M = 1000
          Seed = 123
          loglik= ICA_contcont_long_cre(p=6,
                                       time_points = c(0, 1, 2, 4, 6, 8),
                                       R_rho = 0.000236,
                                       R_structure = "CS",
                                       VT0S0 = 0.009203, VT1S1 = 0.008759,
                                       VT0T0 = 0.01365,  VT1T1 = 0.01379,
                                       VS0S0 = 0.01819,  VS1S1 = 0.01854,
                                       VT0T1 = 0, VT0S1 = 0,
                                       VT1S0 = 0, VS0S1 = 0,
                                       D = D,
                                       N = 568,
                                       M = 10,
                                       Seed = 123)
          expect_equal(loglik$R2_Lambda, 0.9440744, tolerance = 1e-6)})
