library(tibble)
library(dplyr)
library(ggplot2)


baseline_summary <- read.csv("simulation_outputs/baseline_summary.csv")
summary_cv <- read.csv("simulation_outputs/summary_cv.csv")
summary_N1 <- read.csv("simulation_outputs/summary_N1.csv")
summary_ratio <- read.csv("simulation_outputs/summary_ratio.csv")


# Baseline analysis
baseline_summary %>%
  select(
    tablet_cv,
    mean_true_conc,
    mean_estimated_conc,
    mean_theoretical_error,
    empirical_cv,
    mean_bias,
    rmse,
    mean_effort,
    failure_rate
  ) %>%
  knitr::kable(digits = 2, caption = "Baseline simulation results")

# Tablet CV analysis
## Simulation results across tablet CV levels
summary_cv %>%
  select(
    tablet_cv,
    mean_true_conc,
    mean_estimated_conc,
    mean_theoretical_error,
    empirical_cv,
    mean_bias,
    rmse,
    mean_effort,
    failure_rate
  ) %>%
  knitr::kable(digits = 2, caption = "Simulation results across tablet CV levels")

## Plot: tablet CV and total error
ggplot(summary_cv, aes(x = tablet_cv * 100, y = mean_theoretical_error)) +
  geom_line(linewidth = 1, colour = "blue") +
  geom_point(size = 3, colour = "blue") +
  labs(
    title = "Effect of tablet CV on total error",
    x = "Tablet CV (%)",
    y = "Mean total error (%)"
  )

## Error decomposition across tablet CV levels
summary_cv %>%
  select(
    tablet_cv,
    mean_tablet_term,
    mean_calibration_term,
    mean_marker_term,
    mean_theoretical_error
  ) %>%
  knitr::kable(digits = 2, caption = "Error decomposition across tablet CV levels")


summary_cv <- summary_cv %>%
  mutate(
    tablet_share = mean_tablet_term^2 / mean_theoretical_error^2,
    calibration_share = mean_calibration_term^2 / mean_theoretical_error^2,
    marker_share = mean_marker_term^2 / mean_theoretical_error^2
  )

summary_cv %>%
  select(
    tablet_cv,
    tablet_share,
    calibration_share,
    marker_share
  ) %>%
  knitr::kable(digits = 3, caption = "Relative contributions of error components across tablet CV levels")


summary_cv_share_long <- summary_cv %>%
  select(
    tablet_cv,
    tablet_share,
    calibration_share,
    marker_share
  ) %>%
  tidyr::pivot_longer(
    cols = c(tablet_share, calibration_share, marker_share),
    names_to = "component",
    values_to = "share"
  )

ggplot(summary_cv_share_long, aes(x = factor(tablet_cv), y = share, fill = component)) +
  geom_col() +
  labs(
    title = "Relative contributions of error components across tablet CV levels",
    x = "Tablet CV",
    y = "Proportion of total error",
    fill = "Component"
  )


# N1 analysis
## Simulation results across numbers of tablets
summary_N1 %>%
  select(
    N1,
    tablet_cv,
    tablet_mean,
    mean_true_conc,
    mean_estimated_conc,
    mean_theoretical_error,
    empirical_cv,
    mean_bias,
    rmse,
    mean_effort,
    failure_rate
  ) %>%
  knitr::kable(digits = 2, caption = "Simulation results across numbers of tablets")

## Plot: number of tablets and total error
ggplot(summary_N1, aes(x = N1, y = mean_theoretical_error)) +
  geom_line(linewidth = 1, colour = "darkgreen") +
  geom_point(size = 3, colour = "darkgreen") +
  labs(
    title = "Effect of number of tablets on total error",
    x = "Number of tablets (N1)",
    y = "Mean total error (%)"
  )

## Error decomposition across numbers of tablets
summary_N1 %>%
  select(
    N1,
    mean_tablet_term,
    mean_calibration_term,
    mean_marker_term,
    mean_theoretical_error
  ) %>%
  knitr::kable(digits = 2, caption = "Error decomposition across numbers of tablets")

summary_N1 <- summary_N1 %>%
  mutate(
    tablet_share = mean_tablet_term^2 / mean_theoretical_error^2,
    calibration_share = mean_calibration_term^2 / mean_theoretical_error^2,
    marker_share = mean_marker_term^2 / mean_theoretical_error^2
  )

summary_N1 %>%
  select(
    N1,
    tablet_share,
    calibration_share,
    marker_share
  ) %>%
  knitr::kable(digits = 3, caption = "Relative contributions of error components across numbers of tablets")


summary_N1_share_long <- summary_N1 %>%
  select(
    N1,
    tablet_share,
    calibration_share,
    marker_share
  ) %>%
  tidyr::pivot_longer(
    cols = c(tablet_share, calibration_share, marker_share),
    names_to = "component",
    values_to = "share"
  )

ggplot(summary_N1_share_long, aes(x = factor(N1), y = share, fill = component)) +
  geom_col() +
  labs(
    title = "Relative contributions of error components across numbers of tablets",
    x = "Number of tablets (N1)",
    y = "Proportion of total error",
    fill = "Component"
  )


# Target-to-marker ratio analysis

## Simulation results across target-to-marker ratio levels
summary_ratio %>%
  select(
    ratio,
    N1,
    tablet_cv,
    tablet_mean,
    mean_true_conc,
    mean_estimated_conc,
    mean_theoretical_error,
    empirical_cv,
    mean_bias,
    rmse,
    mean_effort,
    failure_rate
  ) %>%
  knitr::kable(digits = 2, caption = "Simulation results across target-to-marker ratio levels")

## Plot: target-to-marker ratio vs total error
ggplot(summary_ratio, aes(x = ratio, y = mean_theoretical_error)) +
  geom_line(linewidth = 1, colour = "purple") +
  geom_point(size = 3, colour = "purple") +
  labs(
    title = "Effect of target-to-marker ratio on total error",
    x = "Target-to-marker ratio",
    y = "Mean total error (%)"
  )

## Error decomposition across target-to-marker ratio levels
summary_ratio %>%
  select(
    ratio,
    mean_tablet_term,
    mean_calibration_term,
    mean_marker_term,
    mean_theoretical_error
  ) %>%
  knitr::kable(digits = 2, caption = "Error decomposition across target-to-marker ratio levels")


summary_ratio <- summary_ratio %>%
  mutate(
    tablet_share = mean_tablet_term^2 / mean_theoretical_error^2,
    calibration_share = mean_calibration_term^2 / mean_theoretical_error^2,
    marker_share = mean_marker_term^2 / mean_theoretical_error^2
  )

summary_ratio %>%
  select(
    ratio,
    tablet_share,
    calibration_share,
    marker_share
  ) %>%
  knitr::kable(digits = 3, caption = "Relative contributions of error components across target-to-marker ratio levels")


summary_ratio_share_long <- summary_ratio %>%
  select(
    ratio,
    tablet_share,
    calibration_share,
    marker_share
  ) %>%
  tidyr::pivot_longer(
    cols = c(tablet_share, calibration_share, marker_share),
    names_to = "component",
    values_to = "share"
  )

ggplot(summary_ratio_share_long, aes(x = factor(ratio), y = share, fill = component)) +
  geom_col() +
  labs(
    title = "Relative contributions of error components across target-to-marker ratio levels",
    x = "Target-to-marker ratio",
    y = "Proportion of total error",
    fill = "Component"
  )






