# =============================================================================
# 04_loess_screening.R
#
# LOESS screening on the training sample, 1990-2018.
#
# Input:
#   data/Total_Tonkm.xlsx
# =============================================================================

# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "forecast",
  "tseries",
  "FinTS",
  "ggplot2",
  "strucchange"
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

library(forecast)
library(tseries)
library(FinTS)
library(ggplot2)
library(strucchange)

# ---- Data --------------------------------------------------------------------



Total_Tonkm_train <- Total_Tonkm[Total_Tonkm$Year <= 2018, ]
years_train <- seq(1990, by = 1, length.out = nrow(Total_Tonkm_train))

ts_train <- ts(
  Total_Tonkm_train$Total,
  start = 1990,
  frequency = 1
)

# =============================================================================
# LOESS CONFIGURATIONS
# =============================================================================

# span = 0.15, degree = 1
loess_trailing_train151 <- loess(
  Total ~ Year,
  data = Total_Tonkm_train,
  span = 0.15,
  degree = 1,
  control = loess.control(surface = "direct")
)

trend_trailing_train151 <- predict(loess_trailing_train151)

ts_trend_loess_trailing_train151 <- ts(
  trend_trailing_train151,
  start = 1990,
  frequency = 1
)

# span = 0.25, degree = 1
loess_trailing_train251 <- loess(
  Total_Tonkm_train$Total ~ years_train,
  span = 0.25,
  degree = 1,
  control = loess.control(surface = "direct")
)

trend_trailing_train251 <- predict(loess_trailing_train251)
residuals_trailing_train251 <- Total_Tonkm_train$Total - trend_trailing_train251

ts_trend_loess_trailing_train251 <- ts(
  trend_trailing_train251,
  start = 1990,
  frequency = 1
)

ts_resid_loess_trailing_train251 <- ts(
  residuals_trailing_train251,
  start = 1990,
  frequency = 1
)

# span = 0.30, degree = 2
loess_trailing_train032 <- loess(
  Total_Tonkm_train$Total ~ years_train,
  span = 0.3,
  degree = 2,
  control = loess.control(surface = "direct")
)

trend_trailing_train032 <- predict(loess_trailing_train032)

ts_trend_loess_trailing_train032 <- ts(
  trend_trailing_train032,
  start = 1990,
  frequency = 1
)

# =============================================================================
# LOESS TREND PLOT
# =============================================================================

figure_loess <- ggplot(
  data.frame(
    Year = as.numeric(time(ts_train)),
    Raw_series = as.numeric(ts_train),
    LOESS_trend = as.numeric(ts_trend_loess_trailing_train251)
  ),
  aes(x = Year)
) +
  geom_line(
    aes(y = Raw_series, color = "Raw series"),
    linewidth = 0.8
  ) +
  geom_line(
    aes(y = LOESS_trend, color = "LOESS trend"),
    linewidth = 0.9
  ) +
  scale_color_manual(
    name = "Series",
    values = c(
      "Raw series" = "grey25",
      "LOESS trend" = "#1F78B4"
    )
  ) +
  labs(
    x = "Year",
    y = "Total freight activity (million tonne-km)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(figure_loess)

# =============================================================================
# DIAGNOSTICS ON THE SELECTED LOESS TREND
# =============================================================================

adf_loess251 <- adf.test(ts_trend_loess_trailing_train251)
kpss_loess251 <- kpss.test(ts_trend_loess_trailing_train251)

print(adf_loess251)
print(kpss_loess251)

# Bai-Perron breakpoints
breaks_train251 <- breakpoints(
  ts_trend_loess_trailing_train251 ~ 1
)

print(summary(breaks_train251))
plot(breaks_train251)

# Recursive CUSUM
years_loess251 <- 1990:2018

cusum_train251 <- efp(
  ts_trend_loess_trailing_train251 ~ 1,
  type = "Rec-CUSUM"
)

plot(cusum_train251)
print(sctest(cusum_train251))

cusum_values_251 <- as.vector(cusum_train251$process)
cusum_bounds_251 <- as.numeric(boundary(cusum_train251))

cusum_df_251 <- data.frame(
  year = years_loess251,
  cusum = cusum_values_251,
  bound = cusum_bounds_251,
  out_of_bounds = abs(cusum_values_251) > cusum_bounds_251
)

print(cusum_df_251)
print(cusum_df_251[cusum_df_251$out_of_bounds, ])

# =============================================================================
# RESIDUAL DIAGNOSTICS
# =============================================================================

checkresiduals(ts_resid_loess_trailing_train251)

lb_resid_lag6 <- Box.test(
  ts_resid_loess_trailing_train251,
  lag = 6,
  type = "Ljung-Box"
)

lb_resid_lag10 <- Box.test(
  ts_resid_loess_trailing_train251,
  lag = 10,
  type = "Ljung-Box"
)

adf_resid_loess251 <- adf.test(ts_resid_loess_trailing_train251)
kpss_resid_loess251 <- kpss.test(ts_resid_loess_trailing_train251)
arch_resid_loess251 <- ArchTest(ts_resid_loess_trailing_train251)
shapiro_resid_loess251 <- shapiro.test(ts_resid_loess_trailing_train251)


print(lb_resid_lag6)
print(lb_resid_lag10)
print(adf_resid_loess251)
print(kpss_resid_loess251)
print(arch_resid_loess251)
print(shapiro_resid_loess251)

summary(ts_resid_loess_trailing_train251)
mean(ts_resid_loess_trailing_train251)

# End of script
