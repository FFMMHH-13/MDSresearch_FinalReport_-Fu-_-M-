library(tibble)
library(dplyr)

# Baseline FOVS simualtion

simulate_fovs_once <- function(
    n_target = 10000,
    N1 = 1,
    tablet_mean = 5000,
    tablet_cv = 0,
    region_side = 1,
    fov_size = 0.1,
    n_cal = 20,
    n_ext = 20,
    omega = 2
) {
  
  # Sample area
  V <- region_side^2
  
  # True concentration
  true_conc <- n_target / V
  
  # Generate target particles
  points <- tibble(
    x = runif(n_target, min = 0, max = region_side),
    y = runif(n_target, min = 0, max = region_side)
  )
  
  # Known average marker dose used in the estimator
  known_marker_dose <- N1 * tablet_mean
  
  # Generate per-tablet marker doses, then sum them
  tablet_sd <- tablet_mean * tablet_cv
  tablet_doses <- rnorm(N1, mean = tablet_mean, sd = tablet_sd)
  tablet_doses <- pmax(round(tablet_doses), 1)
  
  true_marker_total <- sum(tablet_doses)
  
  # Generate marker particles
  markers <- tibble(
    x = runif(true_marker_total, min = 0, max = region_side),
    y = runif(true_marker_total, min = 0, max = region_side)
  )
  
  # Build grid of candidate FOV locations
  n_side <- region_side / fov_size
  if (abs(n_side - round(n_side)) > 1e-10) {
    stop("fov_size must divide the region side exactly.")
  }
  
  all_grid <- expand.grid(
    x = seq(0, region_side - fov_size, by = fov_size),
    y = seq(0, region_side - fov_size, by = fov_size)
  )
  
  # Remove outermost edge cells
  inner_grid <- subset(
    all_grid,
    x > 0 & x < (region_side - fov_size) &
      y > 0 & y < (region_side - fov_size)
  )
  
  if (nrow(inner_grid) < n_cal + n_ext) {
    stop("Not enough inner grid cells for the requested number of FOVs.")
  }
  
  # Define a shared local window so calibration and extrapolation
  # FOVs come from nearby locations
  window_cols <- 7
  window_rows <- 6
  
  x_vals <- sort(unique(inner_grid$x))
  y_vals <- sort(unique(inner_grid$y))
  
  inner_nx <- length(x_vals)
  inner_ny <- length(y_vals)
  
  if (window_cols > inner_nx || window_rows > inner_ny) {
    stop("Shared FOV window is too large for the inner grid.")
  }
  
  start_col <- sample(1:(inner_nx - window_cols + 1), size = 1)
  start_row <- sample(1:(inner_ny - window_rows + 1), size = 1)
  
  window_x <- x_vals[start_col:(start_col + window_cols - 1)]
  window_y <- y_vals[start_row:(start_row + window_rows - 1)]
  
  shared_grid <- subset(
    inner_grid,
    x %in% window_x & y %in% window_y
  )
  
  if (nrow(shared_grid) < n_cal + n_ext) {
    stop("Not enough cells in the shared local FOV region.")
  }
  
  # Select calibration and extrapolation FOVs from the same local region
  cal_idx <- sample(1:nrow(shared_grid), size = n_cal, replace = FALSE)
  cal_fovs <- shared_grid[cal_idx, ]
  
  remaining_grid <- shared_grid[-cal_idx, ]
  ext_idx <- sample(1:nrow(remaining_grid), size = n_ext, replace = FALSE)
  ext_fovs <- remaining_grid[ext_idx, ]
  
  # Count fossils in calibration FOVs
  cal_fossil <- numeric(n_cal)
  
  for (i in 1:n_cal) {
    cal_fossil[i] <- sum(
      points$x >= cal_fovs$x[i] &
        points$x <  cal_fovs$x[i] + fov_size &
        points$y >= cal_fovs$y[i] &
        points$y <  cal_fovs$y[i] + fov_size
    )
  }
  
  # Count markers in extrapolation FOVs
  ext_marker <- numeric(n_ext)
  
  for (i in 1:n_ext) {
    ext_marker[i] <- sum(
      markers$x >= ext_fovs$x[i] &
        markers$x <  ext_fovs$x[i] + fov_size &
        markers$y >= ext_fovs$y[i] &
        markers$y <  ext_fovs$y[i] + fov_size
    )
  }
  
  # Calculate effort (time proxy)
  x_total <- sum(cal_fossil)
  n_count <- sum(ext_marker)
  effort <- (omega * n_cal) + x_total + (omega * n_ext) + n_count
  
  # Calculate concentration and error
  ybar <- mean(cal_fossil)
  xhat <- ybar * n_ext
  
  failed <- (n_count == 0 || ybar == 0)
  
  if (failed) {
    conc <- NA
    tablet_term <- NA
    calibration_term <- NA
    marker_term <- NA
    total_error <- NA
    bias <- NA
    percent_bias <- NA
    s3p <- NA
  } else {
    conc <- (xhat * known_marker_dose) / (n_count * V)
    
    s3p <- sd(cal_fossil) / mean(cal_fossil)
    
    # Error decomposition
    tablet_term <- 100 * (tablet_cv / sqrt(N1))
    calibration_term <- 100 * (s3p / sqrt(n_cal))
    marker_term <- 100 * (sqrt(n_count) / n_count)
    
    total_error <- sqrt(
      tablet_term^2 +
        calibration_term^2 +
        marker_term^2
    )
    
    bias <- conc - true_conc
    percent_bias <- 100 * bias / true_conc
  }
  
  tibble(
    true_conc = true_conc,
    known_marker_dose = known_marker_dose,
    true_marker_total = true_marker_total,
    counted_markers = n_count,
    ybar = ybar,
    s3p = s3p,
    conc = conc,
    tablet_term = tablet_term,
    calibration_term = calibration_term,
    marker_term = marker_term,
    total_error = total_error,
    bias = bias,
    percent_bias = percent_bias,
    effort = effort,
    failed = failed
  )
  
}



# Repeated simulations
simulate_fovs_repeated <- function(
    n_rep = 5000,
    n_target = 10000,
    N1 = 1,
    tablet_mean = 5000,
    tablet_cv = 0,
    region_side = 1,
    fov_size = 0.1,
    n_cal = 20,
    n_ext = 20,
    omega = 2
) {
  
  results <- vector("list", n_rep)
  
  for (k in 1:n_rep) {
    results[[k]] <- simulate_fovs_once(
      n_target = n_target,
      N1 = N1,
      tablet_mean = tablet_mean,
      tablet_cv = tablet_cv,
      region_side = region_side,
      fov_size = fov_size,
      n_cal = n_cal,
      n_ext = n_ext,
      omega = omega
    )
  }
  
  results_df <- bind_rows(results)
  
  mean_conc <- mean(results_df$conc, na.rm = TRUE)
  sd_conc <- sd(results_df$conc, na.rm = TRUE)
  mean_error <- mean(results_df$total_error, na.rm = TRUE)
  sim_cv <- 100 * sd_conc / mean_conc
  rmse <- sqrt(mean((results_df$conc - results_df$true_conc)^2, na.rm = TRUE))
  failure_rate <- mean(results_df$failed)
  
  summary_df <- tibble(
    N1 = N1,
    tablet_mean = tablet_mean,
    tablet_cv = tablet_cv,
    mean_true_conc = mean(results_df$true_conc, na.rm = TRUE),
    mean_estimated_conc = mean_conc,
    empirical_sd = sd_conc,
    empirical_cv = sim_cv,
    mean_tablet_term = mean(results_df$tablet_term, na.rm = TRUE),
    mean_calibration_term = mean(results_df$calibration_term, na.rm = TRUE),
    mean_marker_term = mean(results_df$marker_term, na.rm = TRUE),
    mean_theoretical_error = mean_error,
    mean_bias = mean(results_df$bias, na.rm = TRUE),
    mean_percent_bias = mean(results_df$percent_bias, na.rm = TRUE),
    rmse = rmse,
    mean_effort = mean(results_df$effort, na.rm = TRUE),
    failure_rate = failure_rate
  )
  
  
  list(
    raw = results_df,
    summary = summary_df
  )
}


# Baseline 

set.seed(456)

baseline_result <- simulate_fovs_repeated(
  n_rep = 5000,
  n_target = 10000,
  N1 = 1,
  tablet_mean = 5000,
  tablet_cv = 0,
  region_side = 1,
  fov_size = 0.1,
  n_cal = 20,
  n_ext = 20,
  omega = 2
)

baseline_summary <- baseline_result$summary


# Tablet CV scenarios
set.seed(456)

tablet_cv_levels <- c(0, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11)

results_cv <- vector("list", length(tablet_cv_levels))

for (i in seq_along(tablet_cv_levels)) {
  results_cv[[i]] <- simulate_fovs_repeated(
    n_rep = 5000,
    n_target = 10000,
    N1 = 1,
    tablet_mean = 5000,
    tablet_cv = tablet_cv_levels[i],
    region_side = 1,
    fov_size = 0.1,
    n_cal = 20,
    n_ext = 20,
    omega = 2
  )$summary
}

summary_cv <- bind_rows(results_cv)


# N1 scenarios
set.seed(456)

N1_levels <- c(1, 2, 4, 6)

total_marker_dose <- 5000
tablet_mean_levels <- total_marker_dose / N1_levels

results_N1 <- vector("list", length(N1_levels))

for (i in seq_along(N1_levels)) {
  results_N1[[i]] <- simulate_fovs_repeated(
    n_rep = 5000,
    n_target = 10000,
    N1 = N1_levels[i],
    tablet_mean = tablet_mean_levels[i],
    tablet_cv = 0.04,
    region_side = 1,
    fov_size = 0.1,
    n_cal = 20,
    n_ext = 20,
    omega = 2
  )$summary
}

summary_N1 <- bind_rows(results_N1)


# Target-to-marker ratio scenarios
set.seed(456)

ratio_levels <- c(1, 2, 4, 6, 8, 10, 15)
n_target <- 10000

results_ratio <- vector("list", length(ratio_levels))

for (i in seq_along(ratio_levels)) {
  current_ratio <- ratio_levels[i]
  current_tablet_mean <- n_target / current_ratio
  
  results_ratio[[i]] <- simulate_fovs_repeated(
    n_rep = 5000,
    n_target = n_target,
    N1 = 1,
    tablet_mean = current_tablet_mean,
    tablet_cv = 0.04,
    region_side = 1,
    fov_size = 0.1,
    n_cal = 20,
    n_ext = 20,
    omega = 2
  )$summary
}

summary_ratio <- bind_rows(results_ratio)
summary_ratio$ratio <- ratio_levels


if (!dir.exists("simulation_outputs")) {
  dir.create("simulation_outputs")
}

write.csv(baseline_summary, "simulation_outputs/baseline_summary.csv", row.names = FALSE)
write.csv(summary_cv, "simulation_outputs/summary_cv.csv", row.names = FALSE)
write.csv(summary_N1, "simulation_outputs/summary_N1.csv", row.names = FALSE)
write.csv(summary_ratio, "simulation_outputs/summary_ratio.csv", row.names = FALSE)

