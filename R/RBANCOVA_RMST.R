
#--------------------------------------------------------------------------------#
# Helper Functions
#--------------------------------------------------------------------------------#
#' Build K x L interval-length matrix for RMST
#'
#' Given timepoints τ (length L), total horizon tau_max, and K intervals of
#' equal length, returns a K x L matrix whose column l has `interval_length`
#' in the first floor/round((τ_l / tau_max) * K) rows and 0 below.
#' @keywords internal
#' @noRd
.build_t_tilde <- function(timepoints, tau_max, K, interval_length) {
  L <- length(timepoints)
  cnt <- round((timepoints / tau_max) * K)
  cnt <- pmin(pmax(cnt, 0), K)
  Tmat <- sapply(cnt, function(cn) c(rep(interval_length, cn),
                                     rep(0, K - cn)))
  as.matrix(Tmat)  # K x L
}


#--------------------------------------------------------------------------------#
# Function to calculate single-arm RMST and RMST variance + covariates
#--------------------------------------------------------------------------------#

calculate_rmst_with_covariates  <- function(dataset, time_var, event_var, covariate_vars, timepoints, num_intervals) {

  #---------------------------------------------------------------------#
  # STEP 1
  #---------------------------------------------------------------------#  
  tau_max <- max(timepoints)
  nints <- num_intervals
  endpoints <- c(0, seq(tau_max / num_intervals, tau_max, length.out = num_intervals)) # Specify the exact endpoints you want
  endpts <- endpoints
  
  
  # Convert results to data frame
  #treatment_df <- data.frame(time = time_var, event = event_var)
  treatment_df <- data.frame(time = dataset[[time_var]], event = dataset[[event_var]])
  covariate_df <- dataset[, covariate_vars, drop = FALSE]
  
  # Calculate x_bar: mean of each covariate
  x_bar <- colMeans(covariate_df)  
  
  # Loop over each row of the dataframe
  for (row in 1:nrow(treatment_df)) {
    if (treatment_df$event[row] == 0) {
      # Loop over the intervals to add indicators - now instead of just 1 for at-risk/surviving, we only include the proportion of time followed for those censored
      for (i in seq_along(endpoints)) {
        endpoint <- endpoints[i]
        previous_endpoint <- if (i == 1) 0 else endpoints[i - 1]
        
        # Create y_i column
        treatment_df[row, paste0("y_", i)] <- as.numeric(treatment_df$time[row] >= endpoint)
        
        # Create z_i column
        if (i == 1) {
          treatment_df[row, paste0("z_", i)] <- 1
        } else {
          y_i <- treatment_df[row, paste0("y_", i)]
          y_prev <- treatment_df[row, paste0("y_", i - 1)]
          
          treatment_df[row, paste0("z_", i)] <- ifelse(
            y_i == 0 & y_prev > 0,
            (treatment_df$time[row] - previous_endpoint) / (endpoint - previous_endpoint),
            ifelse(y_prev == 0, 0, 1)
          )
        }
        
      }
    } else {
      # Loop over the intervals to add indicators
      for (i in seq_along(endpoints)) {
        endpoint <- endpoints[i]
        previous_endpoint <- if (i == 1) 0 else endpoints[i - 1]
        
        # Create y_i and z_i columns based on the new logic
        treatment_df[row, paste0("y_", i)] <- ifelse(treatment_df$time[row] >= endpoint, 1,
                                                     ifelse(treatment_df[row, paste0("y_", i - 1)] < 1, 0,
                                                            ifelse(treatment_df$event[row] == 1, 0, 0.5)))
        
        treatment_df[row, paste0("z_", i)] <- ifelse(treatment_df$time[row] >= endpoint, 1,
                                                     ifelse(treatment_df[row, paste0("y_", i - 1)] < 1, 0,
                                                            ifelse(treatment_df$event[row] == 1, 1, 0.5)))
        
      }
    }
  }
  
  # After the above processing, update y_i to match z_i if 0 < z_i < 1
  for (row in 1:nrow(treatment_df)) {
    for (i in seq_along(endpoints)) {
      z_value <- treatment_df[row, paste0("z_", i)]
      if (z_value > 0 & z_value < 1) {
        treatment_df[row, paste0("y_", i)] <- z_value
      }
    }
  }
  
  #---------------------------------------------------------------------#
  # STEP 2
  #---------------------------------------------------------------------#     
  
  ###############################
  # Calculate RMST and RMST_var
  ###############################
  # Combine survival indicators and covariates into a single matrix f_j+
  f_i1s <- treatment_df[, -c(1:2)]  # Exclude time and event columns
  f_bar_1 <- colMeans(f_i1s)
  f_bar_j_plus <- c(f_bar_1, x_bar)
  

  # Check if the last element of f_bar_1 is 0 (nobody at risk in last interval)
  intervals_dropped_z <- 0
  while (length(f_bar_1) > 1 && f_bar_1[length(f_bar_1)] == 0) {
    # Drop the last two columns of f_i1s and the last element of f_bar_1
    f_i1s <- f_i1s[, -((ncol(f_i1s) - 1):ncol(f_i1s))]
    f_bar_1 <- f_bar_1[-((length(f_bar_1) - 1):length(f_bar_1))]
    # Update the number of intervals dropped
    intervals_dropped_z <- intervals_dropped_z + 1
    # Adjust the number of intervals
    num_intervals <- num_intervals - 1
    endpoints <- endpoints[-length(endpoints)]
  }
  
  # Check if the second to last element of f_bar_1 is greater than 0 (meaning the last ybar value)
  intervals_dropped_y <- 0
  second_last_valid <- TRUE
  while (length(f_bar_1) > 1 && f_bar_1[length(f_bar_1) - 1] == 0) {
    # Drop the last two columns of f_i1s and the last element of f_bar_1
    f_i1s <- f_i1s[, -((ncol(f_i1s) - 1):ncol(f_i1s))]
    f_bar_1 <- f_bar_1[-((length(f_bar_1) - 1):length(f_bar_1))]
    # Update the number of intervals dropped
    intervals_dropped_y <- intervals_dropped_y + 1
    # Adjust the number of intervals
    num_intervals <- num_intervals - 1
    endpoints <- endpoints[-length(endpoints)]
    second_last_valid <- FALSE
  }
  
  # Calculate V_f_j+ matrix 
  f_ij <- as.matrix(cbind(f_i1s, covariate_df))  # Matrix of both survival and covariate values
  V_f_j_plus <- matrix(0, ncol = length(f_bar_j_plus), nrow = length(f_bar_j_plus))
  for (i in 1:nrow(f_ij)) {
    diff_vector <- f_ij[i, ] - f_bar_j_plus
    outer_product <- diff_vector %*% t(diff_vector)
    V_f_j_plus <- V_f_j_plus + outer_product
  }
  V_f_j_plus <- V_f_j_plus / (nrow(f_ij) * (nrow(f_ij) - 1))
  
  
  # Get D_1j's
  #D_11 <- diag(f_bar_1)
  D1j <- diag(c(f_bar_1, rep(1, length(x_bar))))
  
  #---------------------------------------------------------------------#
  # STEP 3
  #---------------------------------------------------------------------#     
  # Get A_tilde (A1 here)
  T <- matrix(1, nrow = length(endpoints), ncol = length(endpoints))
  T[upper.tri(T)] <- 0
  vector <- c(1,-1)
  A_tilde <- matrix(0, nrow = length(endpoints), ncol = length(endpoints) * 2)
  for (i in 1:length(endpoints)) {
    for (j in 1:length(endpoints)) {
      A_tilde[i, (2 * j - 1):(2 * j)] <- vector * T[i, j]
    }
  }
  
  # Function to create block diagonal matrix A2
  # Get dimensions of A_tilde
  dim_A_tilde <- dim(A_tilde)
  K <- dim_A_tilde[1]  # Number of rows in A_tilde
  num_cols <- dim_A_tilde[2]  # Number of columns in A_tilde
  
  # Create an identity matrix of size m x m
  m <- length(x_bar)
  I_m <- diag(1, nrow = m, ncol = m)
  
  # Construct A2 matrix by combining A_tilde and I_m
  # A2 will have dimensions (K + m) x (num_cols + m)
  A2 <- rbind(
    cbind(A_tilde, matrix(0, nrow = K, ncol = m)),  # Upper left: A_tilde, Upper right: zero matrix
    cbind(matrix(0, nrow = m, ncol = num_cols), I_m) # Lower left: zero matrix, Lower right: identity matrix
  )
  
  
  # Get S & D_2j's
  S_1 <- as.vector(exp(A_tilde %*% log(f_bar_1)))
  D2j <- diag(c(S_1, rep(1, length(x_bar)))) 
  
  # Get V_G_j's
  J <- D2j %*% A2 %*% solve(D1j)
  V_G_j <- J %*% V_f_j_plus %*% t(J)
  #V_G_1 <- D_2j %*% A_tilde %*% solve(D_11) %*% V_f_1 %*% solve(D_11) %*% t(A_tilde) %*% D_21
  
  #---------------------------------------------------------------------#
  # STEP 4
  #---------------------------------------------------------------------#   
  
  tau_max <- max(timepoints)
  
  L <- length(timepoints)
  
  tau_min_prop = min(timepoints)/tau_max
  
  # Get RMST & V_rmst
  #if (second_last_valid) {
    
    # Create vector of interval lengths
    interval_length <- diff(endpoints)[1]  
    K_eff <- num_intervals                
    
    t_tilde <- .build_t_tilde(
      timepoints = timepoints,
      tau_max    = tau_max,
      K          = K_eff,
      interval_length = interval_length
    )
    #interval_length = endpoints[2]
    #t_tilde <- matrix(c(rep(interval_length, num_intervals)), nrow = num_intervals, ncol = 1)
    #t_tilde <- matrix(c(rep(interval_length, num_intervals*tau_min_prop), rep(0, num_intervals-(num_intervals*tau_min_prop)), rep(interval_length, num_intervals)), nrow = num_intervals, ncol = L)
    
    #A31:
    # Create an identity matrix of size m x m
    m <- length(x_bar)
    I_m <- diag(1, nrow = m, ncol = m)
    A31 <- rbind(
      cbind(t(t_tilde)/2, matrix(0, nrow = L, ncol = m)),  
      cbind(matrix(0, nrow = m, ncol = L), I_m)  # Lower-left: Zero matrix, Lower-right: Identity matrix
    )
    
    # Original A2 matrix
    A32 <- matrix(0, nrow = K_eff, ncol = length(endpoints))
    for (i in 1:K_eff) {
      for (j in 1:length(endpoints)) {
        if ((j == i) || (j == i + 1)) {
          A32[i, j] <- 1
        }
      }
    }
    
    A_32 <- rbind(
      cbind(A32, matrix(0, nrow = nrow(A32), ncol = m)),  
      cbind(matrix(0, nrow = m, ncol = ncol(A32)), I_m)  # Lower-left: Zero matrix, Lower-right: Identity matrix
    )
    
    A_3 <- A31%*%A_32
    
    d_j <- A_3 %*% matrix(c(S_1, x_bar), ncol=1)
    

    RMST <- as.numeric(d_j[seq_len(L)])   #First L entries
    #RMST <- (A_3 %*% matrix(c(S_1, x_bar), ncol=1))[1:2]
    
    V_M_j <- A_3 %*% V_G_j %*% t(A_3)
    RMST_var <- diag(V_M_j)[seq_len(L)]
    #RMST_var <- c(V_M_j[1, 1], V_M_j[2, 2])
    
    
    
    
    
  # } else {
  #   
  #   # Create vector of interval lengths
  #   interval_length = endpoints[2]
  #   t_tilde <- matrix(c(rep(interval_length, num_intervals*tau_min_prop), rep(0, (num_intervals-num_intervals*tau_min_prop)+1), rep(interval_length, (num_intervals + 1))), nrow = (num_intervals + 1), ncol = L)
  #   
  #   #A31:
  #   # Create an identity matrix of size m x m
  #   m <- length(x_bar)
  #   I_m <- diag(1, nrow = m, ncol = m)
  #   A31 <- rbind(
  #     cbind(t(t_tilde)/2, matrix(0, nrow = 1, ncol = m)),  
  #     cbind(matrix(0, nrow = m, ncol = ncol(t(t_tilde))), I_m)  # Lower-left: Zero matrix, Lower-right: Identity matrix
  #   )
  #   
  #   # Original A2 matrix
  #   A32 <- matrix(0, nrow = num_intervals+1, ncol = length(endpoints))
  #   for (i in 1:num_intervals) {
  #     for (j in 1:length(endpoints)) {
  #       if ((j == i) || (j == i + 1)) {
  #         A32[i, j] <- 1
  #       }
  #     }
  #   }
  #   for (i in (num_intervals + 1)) {
  #     for (j in (length(endpoints))) {
  #       A32[i, j] <- 1
  #     }
  #   }
  #   
  #   A_32 <- rbind(
  #     cbind(A32, matrix(0, nrow = nrow(A32), ncol = m)),  
  #     cbind(matrix(0, nrow = m, ncol = ncol(A32)), I_m)  # Lower-left: Zero matrix, Lower-right: Identity matrix
  #   )
  #   
  #   A_3 <- A31%*%A_32
  #   
  #   d_j <- A_3 %*% matrix(c(S_1, x_bar), ncol=1)
  #   RMST <- (A_3 %*% matrix(c(S_1, x_bar), ncol=1))[1]
  #   #RMST <- (t(t_tilde) %*% (head(S_1, -1) + tail(S_1, -1))) / 2 
  #   
  #   V_M_j <- A_3 %*% V_G_j %*% t(A_3)
  #   RMST_var <- c(V_M_j[1, 1], V_M_j[2, 2])
  #   
  #   
  # }
  
  # Now calculate RMST based on continuous time
  # Note that this keeps all intervals (even the ones dropped for our calcs)
  RMST_cont <- as.data.frame(simtrial:::rmst_single_arm(
    time_var = treatment_df$time,
    event_var = treatment_df$event,
    tau = tau_max
  ))
  
  # Store results in the results_df
  results_df_covs <- data.frame(
    sample_size = nrow(treatment_df),
    num_intervals = num_intervals,
    nints = nints,
    RMST = RMST,
    RMST_var = RMST_var,
    RMST_cont = RMST_cont$rmst,
    RMST_cont_var = RMST_cont$variance,
    intervals_dropped_z = intervals_dropped_z,
    intervals_dropped_y = intervals_dropped_y)
  
  
  
  #return(results_df_covs)
  return(list(results_df_covs = results_df_covs, V_M_j = V_M_j, d_j = d_j))
  
  
}


# ------------------------
# Utilities
# ------------------------
`%||%` <- function(x,y) if(!is.null(x)) x else y

.adapt_result <- function(obj) {
  est <- obj$d_j %||% obj$RMSTs
  V   <- obj$V_M_j %||% obj$RMST_vars
  if (is.null(est) || is.null(V))
    stop("Need (d_j, V_M_j) or (RMSTs, RMST_vars) in trt/ctl objects.")
  list(d = as.numeric(est), V = as.matrix(V))
}

.safe_solve <- function(M) {
  out <- try(solve(M), silent = TRUE)
  if (inherits(out, "try-error")) {
    if (!requireNamespace("MASS", quietly = TRUE))
      stop("Matrix is singular. Install 'MASS' for ginv() fallback.")
    return(MASS::ginv(M))
  }
  out
}

._one_ci <- function(est, var, alpha = 0.05) {
  se <- sqrt(max(var, 0))
  z  <- if (se > 0) est / se else NA_real_
  p  <- if (is.na(z)) NA_real_ else 2 * (1 - pnorm(abs(z)))
  crit <- qnorm(1 - alpha/2)
  c(lcl = est - crit*se, ucl = est + crit*se, se = se, z = z, p = p)
}

._fmt <- function(e, l, u, digits = 2) {
  sprintf(paste0("%.",digits,"f (%.",digits,"f, %.",digits,"f)"), e, l, u)
}

#------------------------------------------------------------------------------------#
# Function to estimate difference in RMSTs + CIs, unadjusted and covariate-adjusted
#------------------------------------------------------------------------------------#
calc_diff_unadj_and_wls <- function(trt, ctl,
                                    time_labels,      # length L; defines L
                                    alpha = 0.05,
                                    digits = 2,
                                    contrast_idx = c(1, length(time_labels)),
                                    print_table = TRUE) {
  
  if (length(time_labels) < 1) stop("Provide time_labels (length = L).")
  L <- length(time_labels)
  
  # Pull stacked vector and covariance from each arm
  TRT <- .adapt_result(trt)
  CTL <- .adapt_result(ctl)
  
  # d is length L+K; V is (L+K)x(L+K)
  d_full  <- TRT$d - CTL$d
  V_full  <- TRT$V + CTL$V
  if (length(d_full) != nrow(V_full) || nrow(V_full) != ncol(V_full))
    stop("Dimensions of d and V are inconsistent.")
  
  K <- length(d_full) - L
  if (K < 0) stop("time_labels longer than d vector.")
  
  # --- (1) Unadjusted per-timepoint (use the first L entries)
  d_L  <- d_full[1:L]
  V_LL <- V_full[1:L, 1:L, drop = FALSE]
  
  unadj_rows <- do.call(rbind, lapply(seq_len(L), function(j) {
    s <- ._one_ci(d_L[j], V_LL[j,j], alpha)
    data.frame(
      method   = "unadjusted",
      contrast = "Treated - Control",
      time     = as.character(time_labels[j]),
      estimate = d_L[j],
      est_CI   = ._fmt(d_L[j], s["lcl"], s["ucl"], digits),
      p_value  = unname(s["p"]),
      stringsAsFactors = FALSE
    )
  }))
  
  # Optional unadjusted Δ (between two timepoints)
  delta_unadj <- NULL
  if (length(contrast_idx) == 2) {
    i <- contrast_idx[1]; k <- contrast_idx[2]
    cvec <- numeric(L); cvec[i] <- -1; cvec[k] <-  1
    est  <- as.numeric(c(cvec %*% d_L))
    var  <- as.numeric(t(cvec) %*% V_LL %*% cvec)
    s    <- ._one_ci(est, var, alpha)
    delta_unadj <- data.frame(
      method   = "unadjusted",
      contrast = "Δ(t2 - t1) of (Treated - Control)",
      time     = paste0(time_labels[k], " - ", time_labels[i]),
      estimate = est,
      est_CI   = ._fmt(est, s["lcl"], s["ucl"], digits),
      p_value  = unname(s["p"]),
      stringsAsFactors = FALSE
    )
  }
  
  # --- (2) WLS/GLS-adjusted for covariates using full (L+K)
  # C selects the first L components (timepoint effects) out of L+K
  C <- rbind(diag(L), matrix(0, nrow = K, ncol = L))  # (L+K) x L
  Vinv <- .safe_solve(V_full)
  M    <- t(C) %*% Vinv %*% C
  Minv <- .safe_solve(M)
  
  b_w   <- Minv %*% t(C) %*% Vinv %*% d_full    # Lx1
  V_bw  <- Minv                                  # LxL
  
  b_w   <- as.numeric(b_w)
  
  adj_rows <- do.call(rbind, lapply(seq_len(L), function(j) {
    s <- ._one_ci(b_w[j], V_bw[j,j], alpha)
    data.frame(
      method   = "covariate-adjusted",
      contrast = "Treated - Control",
      time     = as.character(time_labels[j]),
      estimate = b_w[j],
      est_CI   = ._fmt(b_w[j], s["lcl"], s["ucl"], digits),
      p_value  = unname(s["p"]),
      stringsAsFactors = FALSE
    )
  }))
  
  # Optional adjusted Δ
  delta_adj <- NULL
  if (length(contrast_idx) == 2) {
    i <- contrast_idx[1]; k <- contrast_idx[2]
    g <- numeric(L); g[i] <- -1; g[k] <-  1
    est <- as.numeric(c(g %*% b_w))
    var <- as.numeric(t(g) %*% V_bw %*% g)
    s   <- ._one_ci(est, var, alpha)
    delta_adj <- data.frame(
      method   = "covariate-adjusted",
      contrast = "Δ(t2 - t1) of (Treated - Control)",
      time     = paste0(time_labels[k], " - ", time_labels[i]),
      estimate = est,
      est_CI   = ._fmt(est, s["lcl"], s["ucl"], digits),
      p_value  = unname(s["p"]),
      stringsAsFactors = FALSE
    )
  }
  
  # --- Combine pretty table
  pretty_tbl <- rbind(unadj_rows, delta_unadj, adj_rows, delta_adj)
  rownames(pretty_tbl) <- NULL
  
  if (isTRUE(print_table)) {
    print(pretty_tbl[, c("method","contrast","time","est_CI","p_value")], row.names = FALSE)
  }
  
  # Return everything you might want downstream
  list(
    table_pretty = pretty_tbl[, c("method","contrast","time","est_CI","p_value")],
    unadjusted   = list(d = d_L, V = V_LL,
                        per_time = unadj_rows, delta = delta_unadj),
    adjusted     = list(b = b_w, V = V_bw,
                        per_time = adj_rows,  delta = delta_adj),
    meta         = list(L = L, K = K, time_labels = time_labels,
                        contrast_idx = contrast_idx)
  )
}





#------------------------------------------------------------------------------------#
# Pipeline function that estimates RMSTs, Differences, and outputs table of results
#------------------------------------------------------------------------------------#

#' RB-ANCOVA for RMST
#' @param data_treated A data frame subset to the treated units
#' @param data_control A data frame subset to the control units
#' @param event_var Name of the event indicator (0/1).
#' @param covariate_vars list of covariates to adjust for
#' @param timepoints Numeric vector of taus for RMST.
#' @param num_intervals Integer: number of intervals for discretization.
#'
#' @return A data frame with one row per difference estimate
#' @export
#' 
rbancova_rmst <- function(data_treated, data_control, time_var, event_var, covariate_vars, timepoints, num_intervals){
  trt <- calculate_rmst_with_covariates(data_treated, time_var, event_var, covariate_vars, timepoints, num_intervals)
  ctl <- calculate_rmst_with_covariates(data_control, time_var, event_var, covariate_vars, timepoints, num_intervals)
  
  L <- length(timepoints)
  # keep the single-τ contrasts (Treated − Control), skip Δ if L == 1
  ci_idx <- if (L >= 2) c(1, L) else NULL
  
  calc_diff_unadj_and_wls(
    trt, ctl,
    time_labels  = timepoints,  # length L (works with L == 1)
    alpha        = 0.05,
    contrast_idx = ci_idx       # NULL => no Δ row
  )
}




# #diff_tc <- calc_diff_unadj_and_wls( trt, ctl, time_labels = timepoints, 
#   alpha = 0.05,
#   contrast_idx = ifelse(length(timepoints)>1), c(1, 2), NULL )   # Δ(τ2 - τ1); change as needed)
