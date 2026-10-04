# =============================================================================
# 01_raw_series_diagnostics.R
#
# Replication script for raw-series diagnostics and figures.
#
# Input:
#   data/Total_Tonkm.xlsx
# =============================================================================

# ---- Packages ---------------------------------------------------------------

required_packages <- c(
  "tseries",
  "moments",
  "forecast",
  "ggplot2",
  "FinTS",
  "strucchange",
  "gridExtra"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

library(tseries)
library(moments)
library(forecast)
library(ggplot2)
library(FinTS)
library(strucchange)
library(gridExtra)

# ---- Data -------------------------------------------------------------------

df_full <- Total_Tonkm
df_train <- Total_Tonkm[Total_Tonkm$Year <= 2018, ]

ts_total_values <- df_full$Total
ts_train <- df_train$Total

totaltkm_diff1 <- diff(ts_total_values)
totaltkm_diff2 <- diff(totaltkm_diff1)

diff1_train <- diff(ts_train)
diff2_train <- diff(diff1_train)

# =============================================================================
# FIGURE 1
# Total freight volumes (1990–2023), training/test split
# =============================================================================

ts_total <- data.frame(
  Year = seq(1990, by = 1, length.out = length(df_full$Year)),
  Total_tkm = df_full$Total
)

figure_1 <- ggplot(ts_total, aes(x = Year, y = Total_tkm)) +
  geom_line(color = "grey20", linewidth = 0.8) +
  geom_point(color = "grey20", size = 1.8) +
  geom_vline(
    xintercept = 2018,
    linetype = "dashed",
    color = "grey20",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 2014,
    y = max(ts_total$Total_tkm, na.rm = TRUE) * 0.95,
    label = "Training set",
    color = "grey20",
    size = 3.8
  ) +
  annotate(
    "text",
    x = 2021,
    y = max(ts_total$Total_tkm, na.rm = TRUE) * 0.95,
    label = "Test set",
    color = "grey20",
    size = 3.8
  ) +
  labs(
    x = "Year",
    y = "Total freight activity (million tonne-km)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()
  )

print(figure_1)

# =============================================================================
# TABLE 1
# Descriptive statistics: full and training samples
# =============================================================================

summary_stats_full <- data.frame(
  Mean = mean(ts_total_values),
  Median = median(ts_total_values),
  Min = min(ts_total_values),
  Max = max(ts_total_values),
  Variance = var(ts_total_values),
  Standard_Deviation = sd(ts_total_values),
  Skewness = moments::skewness(ts_total_values),
  Kurtosis = moments::kurtosis(ts_total_values)
)

summary_stats_train <- data.frame(
  Mean = mean(ts_train),
  Median = median(ts_train),
  Min = min(ts_train),
  Max = max(ts_train),
  Variance = var(ts_train),
  Standard_Deviation = sd(ts_train),
  Skewness = moments::skewness(ts_train),
  Kurtosis = moments::kurtosis(ts_train)
)

print(summary_stats_full)
print(summary_stats_train)

# =============================================================================
# APPENDIX B, FIGURE 6
# (a) raw series; (b) first-differenced series
# =============================================================================

acf_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Figure 6a
acf_plot_level <- ggAcf(ts_total_values) +
  labs(title = "ACF", x = "Lag", y = "ACF") +
  acf_theme

pacf_plot_level <- ggPacf(ts_total_values) +
  labs(title = "PACF", x = "Lag", y = "PACF") +
  acf_theme

grid.arrange(acf_plot_level, pacf_plot_level, ncol = 2)

# Figure 6b
acf_plot_diff1 <- ggAcf(totaltkm_diff1) +
  labs(title = "ACF", x = "Lag", y = "ACF") +
  acf_theme

pacf_plot_diff1 <- ggPacf(totaltkm_diff1) +
  labs(title = "PACF", x = "Lag", y = "PACF") +
  acf_theme

grid.arrange(acf_plot_diff1, pacf_plot_diff1, ncol = 2)

# =============================================================================
# APPENDIX B, FIGURE 7
# Training-sample ACF/PACF after d = 1 and d = 2
# =============================================================================

# Figure 7a: d = 1
acf_plot_d1 <- ggAcf(diff1_train) +
  labs(title = "ACF", x = "Lag", y = "ACF") +
  acf_theme

pacf_plot_d1 <- ggPacf(diff1_train) +
  labs(title = "PACF", x = "Lag", y = "PACF") +
  acf_theme

grid.arrange(acf_plot_d1, pacf_plot_d1, ncol = 2)

# Figure 7b: d = 2
acf_plot_d2 <- ggAcf(diff2_train) +
  labs(title = "ACF", x = "Lag", y = "ACF") +
  acf_theme

pacf_plot_d2 <- ggPacf(diff2_train) +
  labs(title = "PACF", x = "Lag", y = "PACF") +
  acf_theme

grid.arrange(acf_plot_d2, pacf_plot_d2, ncol = 2)

# =============================================================================
# TABLE 2
# Structural diagnostics: full and training samples
# =============================================================================

# ---- Shapiro-Wilk ------------------------------------------------------------

shapiro_full <- shapiro.test(ts_total_values)
shapiro_train <- shapiro.test(ts_train)

# ---- ADF ---------------------------------------------------------------------

# Augmented Dickey-Fuller test with the package default lag selection.
adf_full <- tseries::adf.test(ts_total_values)
adf_train <- tseries::adf.test(ts_train)

# ---- KPSS --------------------------------------------------------------------

# KPSS level-stationarity test with package defaults.
kpss_full <- tseries::kpss.test(ts_total_values)
kpss_train <- tseries::kpss.test(ts_train)

# ---- Box-Pierce --------------------------------------------------------------

# Box-Pierce tests; lag = floor(sqrt(n)).

bp_level_full <- Box.test(
  ts_total_values,
  lag = floor(sqrt(length(ts_total_values)))
)

bp_level_train <- Box.test(
  ts_train,
  lag = floor(sqrt(length(ts_train)))
)

bp_diff1_full <- Box.test(
  totaltkm_diff1,
  lag = floor(sqrt(length(totaltkm_diff1)))
)

bp_diff1_train <- Box.test(
  diff1_train,
  lag = floor(sqrt(length(diff1_train)))
)

# ---- Ljung-Box ---------------------------------------------------------------

# Ljung-Box tests; lag = 10.

lb_level_full <- Box.test(
  ts_total_values,
  lag = 10,
  type = "Ljung-Box"
)

lb_level_train <- Box.test(
  ts_train,
  lag = 10,
  type = "Ljung-Box"
)

lb_diff1_full <- Box.test(
  totaltkm_diff1,
  lag = 10,
  type = "Ljung-Box"
)

lb_diff1_train <- Box.test(
  diff1_train,
  lag = 10,
  type = "Ljung-Box"
)

# ---- ARCH LM -----------------------------------------------------------------

arch_full <- FinTS::ArchTest(ts_total_values, lags = 5)
arch_train <- FinTS::ArchTest(ts_train, lags = 5)

# ---- Bai-Perron breakpoints --------------------------------------------------

ts_full_bp <- ts(
  ts_total_values,
  start = 1990,
  frequency = 1
)

ts_train_bp <- ts(
  ts_train,
  start = 1990,
  frequency = 1
)

breaks_full_ts <- strucchange::breakpoints(ts_full_bp ~ 1)
breaks_ts_train <- strucchange::breakpoints(ts_train_bp ~ 1)

print(summary(breaks_full_ts))
print(summary(breaks_ts_train))

# 95% confidence intervals for the selected breakpoints.
ci_breaks_full <- confint(
  breaks_full_ts,
  breaks = 3,
  level = 0.95
)

ci_breaks_train <- confint(
  breaks_ts_train,
  breaks = 2,
  level = 0.95
)

print(ci_breaks_full)
print(breakdates(ci_breaks_full))

print(ci_breaks_train)
print(breakdates(ci_breaks_train))

# ---- CUSUM -------------------------------------------------------------------

# Recursive CUSUM with intercept-only specification.

# Full sample
years_full <- df_full$Year
cusum_full <- strucchange::efp(
  ts_total_values ~ 1,
  type = "Rec-CUSUM"
)

plot(cusum_full)
print(sctest(cusum_full))

cusum_values_full <- as.vector(cusum_full$process)
cusum_bounds_full <- as.numeric(boundary(cusum_full))

cusum_df_full <- data.frame(
  year = years_full,
  cusum = cusum_values_full,
  bound = cusum_bounds_full,
  out_of_bounds = abs(cusum_values_full) > cusum_bounds_full
)

print(cusum_df_full)
print(cusum_df_full[cusum_df_full$out_of_bounds, ])

# Training sample
years_train <- df_train$Year
cusum_train <- strucchange::efp(
  ts_train ~ 1,
  type = "Rec-CUSUM"
)

plot(cusum_train)
print(sctest(cusum_train))

cusum_values_train <- as.vector(cusum_train$process)
cusum_bounds_train <- as.numeric(boundary(cusum_train))

cusum_df_train <- data.frame(
  year = years_train,
  cusum = cusum_values_train,
  bound = cusum_bounds_train,
  out_of_bounds = abs(cusum_values_train) > cusum_bounds_train
)

print(cusum_df_train)
print(cusum_df_train[cusum_df_train$out_of_bounds, ])

table2_check <- data.frame(
  Diagnostic = c(
    "Shapiro-Wilk p",
    "ADF p",
    "KPSS statistic",
    "Box-Pierce level p",
    "Ljung-Box level p",
    "Box-Pierce diff1 p",
    "Ljung-Box diff1 p",
    "ARCH LM p"
  ),
  Full = c(
    shapiro_full$p.value,
    adf_full$p.value,
    as.numeric(kpss_full$statistic),
    bp_level_full$p.value,
    lb_level_full$p.value,
    bp_diff1_full$p.value,
    lb_diff1_full$p.value,
    arch_full$p.value
  ),
  Training = c(
    shapiro_train$p.value,
    adf_train$p.value,
    as.numeric(kpss_train$statistic),
    bp_level_train$p.value,
    lb_level_train$p.value,
    bp_diff1_train$p.value,
    lb_diff1_train$p.value,
    arch_train$p.value
  )
)

print(table2_check)

# End of script
