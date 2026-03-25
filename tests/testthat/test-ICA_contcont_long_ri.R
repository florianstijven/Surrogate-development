library(testthat)
library(pracma)

ICA_contcont_long_ri <- function(p = numeric(),
                                 AR1 = TRUE,
                                 AR1_rho = numeric(),
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
                                 DT0S0,
                                 DT1S1,
                                 DT0T0,
                                 DT1T1,
                                 DS0S0,
                                 DS1S1,
                                 DT0T1 = seq(-1, 1, by = 0.1),
                                 DT0S1 = seq(-1, 1, by = 0.1),
                                 DT1S0 = seq(-1, 1, by = 0.1),
                                 DS0S1 = seq(-1, 1, by = 0.1)) {

  validate_scalar <- function(x, name, positive = FALSE) {
    if (length(x) != 1 || !is.finite(x)) {
      stop("`", name, "` must be a single finite numeric value.")
    }

    if (positive && x <= 0) {
      stop("`", name, "` must be strictly positive.")
    }
  }

  validate_covariance_vector <- function(x, name) {
    if (length(x) < 1 || any(!is.finite(x))) {
      stop("`", name, "` must contain at least one finite value.")
    }
  }

  build_covariance_matrix <- function(T0T1, T0S0, T0S1, T1S0, T1S1, S0S1,
                                      T0T0, T1T1, S0S0, S1S1) {
    Sigma_c <- diag(c(T0T0, T1T1, S0S0, S1S1))

    Sigma_c[2, 1] <- Sigma_c[1, 2] <- T0T1
    Sigma_c[3, 1] <- Sigma_c[1, 3] <- T0S0
    Sigma_c[4, 1] <- Sigma_c[1, 4] <- T0S1
    Sigma_c[3, 2] <- Sigma_c[2, 3] <- T1S0
    Sigma_c[4, 2] <- Sigma_c[2, 4] <- T1S1
    Sigma_c[4, 3] <- Sigma_c[3, 4] <- S0S1

    Sigma_c
  }

  min_eigenvalue <- function(mat) {
    min_eval <- try(
      min(eigen(mat, symmetric = TRUE, only.values = TRUE)$values),
      silent = TRUE
    )

    if (inherits(min_eval, "try-error")) {
      return(NA_real_)
    }

    as.numeric(min_eval)
  }

  enumerate_valid_matrices <- function(T0S0_cov, T1S1_cov, T0T0, T1T1, S0S0, S1S1,
                                       T0T1, T0S1, T1S0, S0S1) {
    combins <- expand.grid(
      T0T1 = T0T1,
      T0S0 = T0S0_cov,
      T0S1 = T0S1,
      T1S0 = T1S0,
      T1S1 = T1S1_cov,
      S0S1 = S0S1,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    if (nrow(combins) == 0) {
      return(list(total = 0, valid = combins[0, , drop = FALSE]))
    }

    keep <- logical(nrow(combins))

    for (i in seq_len(nrow(combins))) {
      Sigma_c <- build_covariance_matrix(
        T0T1 = combins$T0T1[i],
        T0S0 = combins$T0S0[i],
        T0S1 = combins$T0S1[i],
        T1S0 = combins$T1S0[i],
        T1S1 = combins$T1S1[i],
        S0S1 = combins$S0S1[i],
        T0T0 = T0T0,
        T1T1 = T1T1,
        S0S0 = S0S0,
        S1S1 = S1S1
      )

      min_eval <- min_eigenvalue(Sigma_c)

      if (is.finite(min_eval) && min_eval > 0) {
        keep[i] <- TRUE
      }
    }

    list(total = nrow(combins), valid = combins[keep, , drop = FALSE])
  }

  p <- as.numeric(p)
  AR1_rho <- as.numeric(AR1_rho)

  if (!isTRUE(AR1)) {
    stop("Currently only the AR(1) structure is implemented, so `AR1` must be TRUE.")
  }

  if (length(p) != 1 || !is.finite(p) || p < 1 || p != round(p)) {
    stop("`p` must be a single positive integer.")
  }

  if (length(AR1_rho) != 1 || !is.finite(AR1_rho) || AR1_rho <= -1 || AR1_rho >= 1) {
    stop("`AR1_rho` must be a single finite value in (-1, 1).")
  }

  scalar_values <- list(
    VT0S0 = VT0S0,
    VT1S1 = VT1S1,
    VT0T0 = VT0T0,
    VT1T1 = VT1T1,
    VS0S0 = VS0S0,
    VS1S1 = VS1S1,
    DT0S0 = DT0S0,
    DT1S1 = DT1S1,
    DT0T0 = DT0T0,
    DT1T1 = DT1T1,
    DS0S0 = DS0S0,
    DS1S1 = DS1S1
  )

  variance_names <- c("VT0T0", "VT1T1", "VS0S0", "VS1S1", "DT0T0", "DT1T1", "DS0S0", "DS1S1")

  for (name in names(scalar_values)) {
    validate_scalar(scalar_values[[name]], name, positive = name %in% variance_names)
  }

  covariance_vectors <- list(
    VT0T1 = VT0T1,
    VT0S1 = VT0S1,
    VT1S0 = VT1S0,
    VS0S1 = VS0S1,
    DT0T1 = DT0T1,
    DT0S1 = DT0S1,
    DT1S0 = DT1S0,
    DS0S1 = DS0S1
  )

  for (name in names(covariance_vectors)) {
    validate_covariance_vector(covariance_vectors[[name]], name)
  }

  alpha <- (p - (p * AR1_rho) + (2 * AR1_rho)) / (1 + AR1_rho)
  Q <- rbind(c(-1, 1, 0, 0), c(0, 0, -1, 1))

  v_candidates <- enumerate_valid_matrices(
    T0S0_cov = VT0S0,
    T1S1_cov = VT1S1,
    T0T0 = VT0T0,
    T1T1 = VT1T1,
    S0S0 = VS0S0,
    S1S1 = VS1S1,
    T0T1 = VT0T1,
    T0S1 = VT0S1,
    T1S0 = VT1S0,
    S0S1 = VS0S1
  )

  d_candidates <- enumerate_valid_matrices(
    T0S0_cov = DT0S0,
    T1S1_cov = DT1S1,
    T0T0 = DT0T0,
    T1T1 = DT1T1,
    S0S0 = DS0S0,
    S1S1 = DS1S1,
    T0T1 = DT0T1,
    T0S1 = DT0S1,
    T1S0 = DT1S0,
    S0S1 = DS0S1
  )

  r2_lambda <- numeric(0)
  num_pos_def_v <- nrow(v_candidates$valid)
  num_pos_def_d <- nrow(d_candidates$valid)

  if (nrow(v_candidates$valid) > 0 && nrow(d_candidates$valid) > 0) {
    for (i in seq_len(nrow(v_candidates$valid))) {
      V_c <- build_covariance_matrix(
        T0T1 = v_candidates$valid$T0T1[i],
        T0S0 = v_candidates$valid$T0S0[i],
        T0S1 = v_candidates$valid$T0S1[i],
        T1S0 = v_candidates$valid$T1S0[i],
        T1S1 = v_candidates$valid$T1S1[i],
        S0S1 = v_candidates$valid$S0S1[i],
        T0T0 = VT0T0,
        T1T1 = VT1T1,
        S0S0 = VS0S0,
        S1S1 = VS1S1
      )

      sigma_u <- Q %*% V_c %*% t(Q)
      min_eval_u <- min_eigenvalue(sigma_u)

      if (!is.finite(min_eval_u) || min_eval_u <= 0) {
        next
      }

      for (j in seq_len(nrow(d_candidates$valid))) {
        D_c <- build_covariance_matrix(
          T0T1 = d_candidates$valid$T0T1[j],
          T0S0 = d_candidates$valid$T0S0[j],
          T0S1 = d_candidates$valid$T0S1[j],
          T1S0 = d_candidates$valid$T1S0[j],
          T1S1 = d_candidates$valid$T1S1[j],
          S0S1 = d_candidates$valid$S0S1[j],
          T0T0 = DT0T0,
          T1T1 = DT1T1,
          S0S0 = DS0S0,
          S1S1 = DS1S1
        )

        sigma_d <- Q %*% D_c %*% t(Q)
        min_eval_d <- min_eigenvalue(sigma_d)

        if (!is.finite(min_eval_d) || min_eval_d <= 0) {
          next
        }

        u11 <- sigma_u[1, 1]
        u12 <- sigma_u[1, 2]
        u22 <- sigma_u[2, 2]
        d11 <- sigma_d[1, 1]
        d12 <- sigma_d[1, 2]
        d22 <- sigma_d[2, 2]

        rho_u <- u12 / sqrt(u11 * u22)
        rho_u2 <- rho_u^2

        denom_ud <- (u11 + alpha * d11) * (u22 + alpha * d22)
        if (!is.finite(denom_ud) || denom_ud <= 0) {
          next
        }

        rho_ud <- (u12 + alpha * d12) / sqrt(denom_ud)
        rho_ud2 <- rho_ud^2
        R2_Lambda <- 1 - (1 - rho_u2)^(p - 1) * (1 - rho_ud2)

        if (all(is.finite(c(rho_u, rho_u2, rho_ud, rho_ud2, R2_Lambda)))) {
          r2_lambda <- c(r2_lambda, R2_Lambda)
        }
      }
    }
  }

  num_pos_def_pairs <- length(r2_lambda)

  fit <- list(
    Num.Pos.Def.V = num_pos_def_v,
    Num.Pos.Def.D = num_pos_def_d,
    Num.Pos.Def.Pairs = num_pos_def_pairs,
    R2_Lambda = r2_lambda,
    Call = match.call()
  )

  class(fit) <- "ICA_contcont_long_ri"
  fit
}


test_that("ICA_contcont_long_ri",
          {p=6
          AR1 = TRUE
          AR1_rho = 0.6223
          VT0S0 = 0.03401
          VT1S1 = 0.03640
          VT0T0 = 0.03804
          VT1T1 = 0.04099
          VS0S0 = 0.03485
          VS1S1 = 0.03653
          VT0T1 = 0
          VT0S1 = 0
          VT1S0 = 0
          VS0S1 = 0
          DT0S0 = 0.03497
          DT1S1 = 0.03412
          DT0T0 = 0.04034
          DT1T1 = 0.03470
          DS0S0 = 0.03492
          DS1S1 = 0.03800
          DT0T1 = 0
          DT0S1 = 0
          DT1S0 = 0
          DS0S1 = 0
          loglik= ICA_contcont_long_ri(p=6,
                                       AR1 = TRUE, AR1_rho = 0.6223,
                                       VT0S0 = 0.03401, VT1S1 = 0.03640,
                                       VT0T0 = 0.03804, VT1T1 = 0.04099,
                                       VS0S0 = 0.03485, VS1S1 = 0.03653,
                                       VT0T1 = 0, VT0S1 = 0,
                                       VT1S0 = 0, VS0S1 = 0,
                                       DT0S0 = 0.03497, DT1S1 = 0.03412,
                                       DT0T0 = 0.04034, DT1T1 = 0.03470,
                                       DS0S0 = 0.03492, DS1S1 = 0.03800,
                                       DT0T1 = 0, DT0S1 = 0,
                                       DT1S0 = 0, DS0S1 = 0 )
          expect_equal(loglik$R2_Lambda, 0.9999967, tolerance = 1e-6)})
