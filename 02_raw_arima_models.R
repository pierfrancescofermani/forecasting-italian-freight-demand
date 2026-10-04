# =============================================================================
# 02_raw_arima_models.R
#
# ARIMA estimation, diagnostics, forecasting, and differencing sensitivity
# for annual freight activity, 1990-2023.
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
  "dplyr"
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
library(dplyr)

# ---- Data --------------------------------------------------------------------


ts_tonkm <- ts(
  Total_Tonkm$Total,
  start = 1990,
  end = 2023,
  frequency = 1
)

ts_train <- window(ts_tonkm, start = 1990, end = 2018)
ts_test  <- window(ts_tonkm, start = 2019, end = 2023)

# =============================================================================
# BASELINE MODEL SELECTION
# =============================================================================

auto_arima_aicc <- auto.arima(
  ts_train,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "aicc"
)

auto_arima_aic <- auto.arima(
  ts_train,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "aic"
)

auto_arima_bic <- auto.arima(
  ts_train,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "bic"
)

print(summary(auto_arima_aicc))
print(summary(auto_arima_aic))
print(summary(auto_arima_bic))

selection_summary <- data.frame(
  Criterion = c("AICc", "AIC", "BIC"),
  Model = c(
    paste0("ARIMA(", paste(arimaorder(auto_arima_aicc)[1:3], collapse = ","), ")"),
    paste0("ARIMA(", paste(arimaorder(auto_arima_aic)[1:3], collapse = ","), ")"),
    paste0("ARIMA(", paste(arimaorder(auto_arima_bic)[1:3], collapse = ","), ")")
  ),
  AIC = c(
    AIC(auto_arima_aicc),
    AIC(auto_arima_aic),
    AIC(auto_arima_bic)
  ),
  AICc = c(
    auto_arima_aicc$aicc,
    auto_arima_aic$aicc,
    auto_arima_bic$aicc
  ),
  BIC = c(
    BIC(auto_arima_aicc),
    BIC(auto_arima_aic),
    BIC(auto_arima_bic)
  ),
  Residual_Variance = c(
    auto_arima_aicc$sigma2,
    auto_arima_aic$sigma2,
    auto_arima_bic$sigma2
  )
)

print(selection_summary)

# Distinct candidate models selected by the information criteria
arima_100 <- auto_arima_aicc
arima_103 <- auto_arima_aic

# =============================================================================
# RESIDUAL DIAGNOSTICS
# =============================================================================

# ARIMA(1,0,0)
checkresiduals(arima_100)

res_100 <- residuals(arima_100)

shapiro_100 <- shapiro.test(res_100)
ks_100 <- ks.test(
  res_100,
  "pnorm",
  mean(res_100),
  sd(res_100)
)
adf_100 <- adf.test(res_100)
kpss_100 <- kpss.test(res_100, null = "Level")
arch_100 <- ArchTest(res_100)
var_100 <- var(res_100)

print(shapiro_100)
print(ks_100)
print(adf_100)
print(kpss_100)
print(arch_100)
print(var_100)

# ARIMA(1,0,3)
checkresiduals(arima_103)

res_103 <- residuals(arima_103)

shapiro_103 <- shapiro.test(res_103)
ks_103 <- ks.test(
  res_103,
  "pnorm",
  mean(res_103),
  sd(res_103)
)
adf_103 <- adf.test(res_103)
kpss_103 <- kpss.test(res_103, null = "Level")
arch_103 <- ArchTest(res_103)
var_103 <- var(res_103)

print(shapiro_103)
print(ks_103)
print(adf_103)
print(kpss_103)
print(arch_103)
print(var_103)

# =============================================================================
# BASELINE OUT-OF-SAMPLE FORECASTS, 2019-2023
# =============================================================================

fc_100 <- forecast(
  arima_100,
  h = length(ts_test),
  level = 95
)

fc_103 <- forecast(
  arima_103,
  h = length(ts_test),
  level = 95
)

baseline_forecasts <- data.frame(
  Year = as.numeric(time(ts_test)),
  Observed = as.numeric(ts_test),
  Forecast_100 = as.numeric(fc_100$mean),
  Lower_100 = as.numeric(fc_100$lower[, "95%"]),
  Upper_100 = as.numeric(fc_100$upper[, "95%"]),
  Forecast_103 = as.numeric(fc_103$mean),
  Lower_103 = as.numeric(fc_103$lower[, "95%"]),
  Upper_103 = as.numeric(fc_103$upper[, "95%"])
)

print(baseline_forecasts)

mae <- function(e) mean(abs(e))
rmse <- function(e) sqrt(mean(e^2))
mape <- function(e, y) mean(abs(e / y)) * 100
forecast_error_variance <- function(e) var(e)

error_100 <- baseline_forecasts$Observed - baseline_forecasts$Forecast_100
error_103 <- baseline_forecasts$Observed - baseline_forecasts$Forecast_103

baseline_error_metrics <- data.frame(
  Model = c("ARIMA(1,0,0)", "ARIMA(1,0,3)"),
  MAE = c(
    mae(error_100),
    mae(error_103)
  ),
  RMSE = c(
    rmse(error_100),
    rmse(error_103)
  ),
  MAPE = c(
    mape(error_100, baseline_forecasts$Observed),
    mape(error_103, baseline_forecasts$Observed)
  ),
  Forecast_Error_Variance = c(
    forecast_error_variance(error_100),
    forecast_error_variance(error_103)
  )
)

print(baseline_error_metrics)

# Baseline fitted values
df_fitted_baseline <- bind_rows(
  data.frame(
    Year = as.numeric(time(ts_train)),
    Observed = as.numeric(ts_train),
    Fitted = as.numeric(fitted(arima_100)),
    Model = "ARIMA(1,0,0)"
  ),
  data.frame(
    Year = as.numeric(time(ts_train)),
    Observed = as.numeric(ts_train),
    Fitted = as.numeric(fitted(arima_103)),
    Model = "ARIMA(1,0,3)"
  )
)

p_fitted_baseline <- ggplot(df_fitted_baseline, aes(x = Year)) +
  geom_line(aes(y = Observed, color = "Observed"), linewidth = 0.8) +
  geom_line(
    aes(y = Fitted, color = "Fitted"),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  facet_wrap(~Model, scales = "fixed") +
  scale_color_manual(
    name = "Series",
    values = c(
      "Observed" = "grey25",
      "Fitted" = "#1F78B4"
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
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(p_fitted_baseline)

# Baseline forecast plot
years_all <- as.numeric(time(ts_tonkm))
years_test <- as.numeric(time(ts_test))

make_forecast_df <- function(fc, model_name) {
  df <- data.frame(
    Year = years_all,
    Observed = as.numeric(ts_tonkm),
    Forecast = NA_real_,
    Lower = NA_real_,
    Upper = NA_real_,
    Model = model_name
  )

  df[df$Year %in% years_test, "Forecast"] <- as.numeric(fc$mean)
  df[df$Year %in% years_test, "Lower"] <- as.numeric(fc$lower[, "95%"])
  df[df$Year %in% years_test, "Upper"] <- as.numeric(fc$upper[, "95%"])

  df
}

df_forecast_baseline <- bind_rows(
  make_forecast_df(fc_100, "ARIMA(1,0,0)"),
  make_forecast_df(fc_103, "ARIMA(1,0,3)")
)

p_forecast_baseline <- ggplot(df_forecast_baseline, aes(x = Year)) +
  geom_vline(
    xintercept = 2018.5,
    linetype = "dashed",
    color = "grey25",
    linewidth = 0.7
  ) +
  geom_ribbon(
    data = df_forecast_baseline %>% filter(!is.na(Forecast)),
    aes(ymin = Lower, ymax = Upper),
    fill = "#D95F02",
    alpha = 0.15
  ) +
  geom_line(
    aes(y = Observed, color = "Observed"),
    linewidth = 0.8
  ) +
  geom_line(
    data = df_forecast_baseline %>% filter(!is.na(Forecast)),
    aes(y = Forecast, color = "Forecasted"),
    linewidth = 0.8
  ) +
  facet_wrap(~Model, scales = "fixed") +
  scale_color_manual(
    name = "Series",
    values = c(
      "Observed" = "grey25",
      "Forecasted" = "#D95F02"
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
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(p_forecast_baseline)

# =============================================================================
# FIXED DIFFERENCING: d = 1
# =============================================================================

auto_arima_d1_aicc <- auto.arima(
  ts_train,
  d = 1,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "aicc"
)

auto_arima_d1_aic <- auto.arima(
  ts_train,
  d = 1,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "aic"
)

auto_arima_d1_bic <- auto.arima(
  ts_train,
  d = 1,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "bic"
)

print(summary(auto_arima_d1_aicc))
print(summary(auto_arima_d1_aic))
print(summary(auto_arima_d1_bic))

# ARIMA(0,1,0), without drift
arima_010 <- Arima(
  ts_train,
  order = c(0, 1, 0),
  include.drift = FALSE
)

fc_010 <- forecast(
  arima_010,
  h = length(ts_test),
  level = 95
)

error_010 <- as.numeric(ts_test) - as.numeric(fc_010$mean)

d1_error_metrics <- data.frame(
  Model = "ARIMA(0,1,0)",
  MAE = mae(error_010),
  RMSE = rmse(error_010),
  MAPE = mape(error_010, as.numeric(ts_test)),
  Forecast_Error_Variance = forecast_error_variance(error_010)
)

print(summary(arima_010))
print(
  data.frame(
    Year = as.numeric(time(ts_test)),
    Observed = as.numeric(ts_test),
    Forecast = as.numeric(fc_010$mean),
    Lower = as.numeric(fc_010$lower[, "95%"]),
    Upper = as.numeric(fc_010$upper[, "95%"])
  )
)
print(d1_error_metrics)

# =============================================================================
# FIXED DIFFERENCING: d = 2
# =============================================================================

auto_arima_d2_aicc <- auto.arima(
  ts_train,
  d = 2,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "aicc"
)

auto_arima_d2_aic <- auto.arima(
  ts_train,
  d = 2,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "aic"
)

auto_arima_d2_bic <- auto.arima(
  ts_train,
  d = 2,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "bic"
)

print(summary(auto_arima_d2_aicc))
print(summary(auto_arima_d2_aic))
print(summary(auto_arima_d2_bic))

arima_220 <- Arima(
  ts_train,
  order = c(2, 2, 0)
)

arima_223 <- Arima(
  ts_train,
  order = c(2, 2, 3)
)

print(summary(arima_220))
print(summary(arima_223))

fc_220 <- forecast(
  arima_220,
  h = length(ts_test),
  level = 95
)

fc_223 <- forecast(
  arima_223,
  h = length(ts_test),
  level = 95
)

error_220 <- as.numeric(ts_test) - as.numeric(fc_220$mean)
error_223 <- as.numeric(ts_test) - as.numeric(fc_223$mean)

d2_error_metrics <- data.frame(
  Model = c("ARIMA(2,2,0)", "ARIMA(2,2,3)"),
  MAE = c(
    mae(error_220),
    mae(error_223)
  ),
  RMSE = c(
    rmse(error_220),
    rmse(error_223)
  ),
  MAPE = c(
    mape(error_220, as.numeric(ts_test)),
    mape(error_223, as.numeric(ts_test))
  ),
  Forecast_Error_Variance = c(
    forecast_error_variance(error_220),
    forecast_error_variance(error_223)
  )
)

print(d2_error_metrics)

d2_forecasts <- data.frame(
  Year = as.numeric(time(ts_test)),
  Observed = as.numeric(ts_test),
  Forecast_220 = as.numeric(fc_220$mean),
  Lower_220 = as.numeric(fc_220$lower[, "95%"]),
  Upper_220 = as.numeric(fc_220$upper[, "95%"]),
  Forecast_223 = as.numeric(fc_223$mean),
  Lower_223 = as.numeric(fc_223$lower[, "95%"]),
  Upper_223 = as.numeric(fc_223$upper[, "95%"])
)

print(d2_forecasts)

# Fitted values for d = 2 specifications
df_fitted_d2 <- bind_rows(
  data.frame(
    Year = as.numeric(time(ts_train)),
    Observed = as.numeric(ts_train),
    Fitted = as.numeric(fitted(arima_220)),
    Model = "ARIMA(2,2,0)"
  ),
  data.frame(
    Year = as.numeric(time(ts_train)),
    Observed = as.numeric(ts_train),
    Fitted = as.numeric(fitted(arima_223)),
    Model = "ARIMA(2,2,3)"
  )
)

p_fitted_d2 <- ggplot(df_fitted_d2, aes(x = Year)) +
  geom_line(aes(y = Observed, color = "Observed"), linewidth = 0.8) +
  geom_line(
    aes(y = Fitted, color = "Fitted"),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  facet_wrap(~Model, scales = "fixed") +
  scale_color_manual(
    name = "Series",
    values = c(
      "Observed" = "grey25",
      "Fitted" = "#1F78B4"
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
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(p_fitted_d2)

df_forecast_d2 <- bind_rows(
  make_forecast_df(fc_220, "ARIMA(2,2,0)"),
  make_forecast_df(fc_223, "ARIMA(2,2,3)")
)

p_forecast_d2 <- ggplot(df_forecast_d2, aes(x = Year)) +
  geom_vline(
    xintercept = 2018.5,
    linetype = "dashed",
    color = "grey25",
    linewidth = 0.7
  ) +
  geom_ribbon(
    data = df_forecast_d2 %>% filter(!is.na(Forecast)),
    aes(ymin = Lower, ymax = Upper),
    fill = "#D95F02",
    alpha = 0.15
  ) +
  geom_line(
    aes(y = Observed, color = "Observed"),
    linewidth = 0.8
  ) +
  geom_line(
    data = df_forecast_d2 %>% filter(!is.na(Forecast)),
    aes(y = Forecast, color = "Forecasted"),
    linewidth = 0.8,
    linetype = "solid"
  ) +
  facet_wrap(~Model, scales = "fixed") +
  scale_color_manual(
    name = "Series",
    values = c(
      "Observed" = "grey25",
      "Forecasted" = "#D95F02"
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
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(p_forecast_d2)

# =============================================================================
# FORECAST-ACCURACY SUMMARY
# =============================================================================

forecast_accuracy <- rbind(
  baseline_error_metrics,
  d1_error_metrics,
  d2_error_metrics
)

rownames(forecast_accuracy) <- NULL
print(forecast_accuracy)

# End of script
