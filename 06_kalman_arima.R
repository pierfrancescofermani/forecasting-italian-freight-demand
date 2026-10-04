# =============================================================================
# 06_kalman_arima.R
#
# ARIMA models on the baseline Kalman-smoothed trend.
#
# Input:
#   data/Total_Tonkm.xlsx
# =============================================================================

# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "forecast",
  "KFAS",
  "tseries",
  "FinTS",
  "Metrics",
  "ggplot2"
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
library(KFAS)
library(tseries)
library(FinTS)
library(Metrics)
library(ggplot2)

# ---- Data --------------------------------------------------------------------

ts_total <- ts(
  Total_Tonkm$Total,
  start = Total_Tonkm$Year[1],
  frequency = 1
)

ts_kf_train <- window(
  ts_total,
  start = 1990,
  end = 2018
)

ts_kf_full <- window(
  ts_total,
  start = 1990,
  end = 2023
)

# =============================================================================
# BASELINE KALMAN-SMOOTHED TRENDS
# =============================================================================

# ---- Training sample ----------------------------------------------------------

Q_level <- 9999999
Q_slope <- 900000
fixed_H <- 9999999

initial_level <- 191099.7
initial_slope <- 450

P1_level <- 200
P1_slope <- 10

model_kf_train_manual2 <- SSModel(
  ts_kf_train ~ SSMtrend(
    degree = 2,
    Q = list(
      matrix(Q_level, 1, 1),
      matrix(Q_slope, 1, 1)
    )
  ),
  H = array(
    as.double(fixed_H),
    dim = c(1, 1, 1)
  )
)

model_kf_train_manual2$a1 <- matrix(
  as.double(c(initial_level, initial_slope)),
  ncol = 1
)

model_kf_train_manual2$P1 <- diag(
  as.double(c(P1_level, P1_slope))
)

model_kf_train_manual2$P1inf <- matrix(
  0,
  2,
  2
)

smoothed_kf_train_manual2 <- KFS(
  model_kf_train_manual2,
  smoothing = "state"
)

train_kf_manual2 <- smoothed_kf_train_manual2$alphahat[, 1]

ts_train_kf_manual2 <- ts(
  train_kf_manual2,
  start = 1990,
  frequency = 1
)

# ---- Full sample --------------------------------------------------------------

Q_level <- 9999999
Q_slope <- 500000
fixed_H <- 9999999

initial_level <- 189781.9
initial_slope <- 350

P1_level <- 500
P1_slope <- 10

model_kf_full_manual2 <- SSModel(
  ts_kf_full ~ SSMtrend(
    degree = 2,
    Q = list(
      matrix(as.double(Q_level)),
      matrix(as.double(Q_slope))
    )
  ),
  H = array(
    as.double(fixed_H),
    dim = c(1, 1, 1)
  )
)

model_kf_full_manual2$a1 <- matrix(
  as.double(c(initial_level, initial_slope)),
  ncol = 1
)

model_kf_full_manual2$P1 <- diag(
  as.double(c(P1_level, P1_slope))
)

model_kf_full_manual2$P1inf <- matrix(
  0,
  2,
  2
)

smoothed_kf_full_manual2 <- KFS(
  model_kf_full_manual2,
  smoothing = "state"
)

full_kf_manual2 <- smoothed_kf_full_manual2$alphahat[, 1]

ts_full_kf_manual2 <- ts(
  full_kf_manual2,
  start = 1990,
  frequency = 1
)

# =============================================================================
# PRE-FIT DIAGNOSTICS
# =============================================================================

train_kf_manual2_diff1 <- diff(
  train_kf_manual2
)

train_kf_manual2_diff2 <- diff(
  diff(train_kf_manual2)
)

adf.test(train_kf_manual2)
kpss.test(train_kf_manual2)
acf(train_kf_manual2)
pacf(train_kf_manual2)
Box.test(train_kf_manual2)

adf.test(train_kf_manual2_diff1)
kpss.test(train_kf_manual2_diff1)
acf(train_kf_manual2_diff1)
pacf(train_kf_manual2_diff1)
Box.test(train_kf_manual2_diff1)

adf.test(train_kf_manual2_diff2)
kpss.test(train_kf_manual2_diff2)
acf(train_kf_manual2_diff2)
pacf(train_kf_manual2_diff2)
Box.test(train_kf_manual2_diff2)

# =============================================================================
# AUTOMATIC MODEL SELECTION
# =============================================================================

# ---- d selected automatically -------------------------------------------------

model_aic <- auto.arima(
  train_kf_manual2,
  ic = "aic",
  stepwise = FALSE,
  approximation = FALSE
)

cat("AIC-selected model:\n")
print(model_aic)

model_aicc <- auto.arima(
  train_kf_manual2,
  ic = "aicc",
  stepwise = FALSE,
  approximation = FALSE
)

cat("\nAICc-selected model:\n")
print(model_aicc)

model_bic <- auto.arima(
  train_kf_manual2,
  ic = "bic",
  stepwise = FALSE,
  approximation = FALSE
)

cat("\nBIC-selected model:\n")
print(model_bic)

# ---- d = 1 -------------------------------------------------------------------

model_d1_aic <- auto.arima(
  train_kf_manual2,
  d = 1,
  ic = "aic",
  stepwise = FALSE,
  approximation = FALSE
)

cat("\nAIC-selected model with d = 1:\n")
print(model_d1_aic)

model_d1_aicc <- auto.arima(
  train_kf_manual2,
  d = 1,
  ic = "aicc",
  stepwise = FALSE,
  approximation = FALSE
)

cat("\nAICc-selected model with d = 1:\n")
print(model_d1_aicc)

model_d1_bic <- auto.arima(
  train_kf_manual2,
  d = 1,
  ic = "bic",
  stepwise = FALSE,
  approximation = FALSE
)

cat("\nBIC-selected model with d = 1:\n")
print(model_d1_bic)

# =============================================================================
# RETAINED ARIMA SPECIFICATIONS
# =============================================================================

fit_auto200 <- Arima(
  train_kf_manual2,
  order = c(2, 0, 0)
)

fit_autoaic013 <- Arima(
  train_kf_manual2,
  order = c(0, 1, 3)
)

fit_auto110 <- Arima(
  train_kf_manual2,
  order = c(1, 1, 0)
)

fit_211 <- Arima(
  train_kf_manual2,
  order = c(2, 1, 1)
)

cat("\n============================================================\n")
cat("ARIMA MODEL SUMMARIES\n")
cat("============================================================\n")

print(fit_auto200)
print(fit_autoaic013)
print(fit_auto110)
print(fit_211)

# =============================================================================
# RESIDUAL DIAGNOSTICS
# =============================================================================

# ---- ARIMA(0,1,3) -------------------------------------------------------------

resid_013 <- residuals(
  fit_autoaic013
)

cat("\n============================================================\n")
cat("ARIMA(0,1,3) RESIDUAL DIAGNOSTICS\n")
cat("============================================================\n")

checkresiduals(resid_013)
summary(resid_013)
acf(
  resid_013,
  main = "ACF of ARIMA(0,1,3) Residuals"
)
pacf(
  resid_013,
  main = "PACF of ARIMA(0,1,3) Residuals"
)
adf.test(resid_013)
kpss.test(resid_013)
ArchTest(resid_013)
shapiro.test(resid_013)
var(resid_013)

# ---- ARIMA(2,0,0) -------------------------------------------------------------

resid_200 <- residuals(
  fit_auto200
)

cat("\n============================================================\n")
cat("ARIMA(2,0,0) RESIDUAL DIAGNOSTICS\n")
cat("============================================================\n")

checkresiduals(resid_200)
summary(resid_200)
acf(
  resid_200,
  main = "ACF of ARIMA(2,0,0) Residuals"
)
pacf(
  resid_200,
  main = "PACF of ARIMA(2,0,0) Residuals"
)
adf.test(resid_200)
kpss.test(resid_200)
ArchTest(resid_200)
shapiro.test(resid_200)
var(resid_200)

# ---- ARIMA(2,1,1) -------------------------------------------------------------

resid_211 <- residuals(
  fit_211
)

cat("\n============================================================\n")
cat("ARIMA(2,1,1) RESIDUAL DIAGNOSTICS\n")
cat("============================================================\n")

checkresiduals(resid_211)
summary(resid_211)
acf(
  resid_211,
  main = "ACF of ARIMA(2,1,1) Residuals"
)
pacf(
  resid_211,
  main = "PACF of ARIMA(2,1,1) Residuals"
)
adf.test(resid_211)
kpss.test(resid_211)
ArchTest(resid_211)
shapiro.test(resid_211)
var(resid_211)

# =============================================================================
# FORECASTS
# =============================================================================

forecast_013 <- forecast(
  fit_autoaic013,
  h = 5
)

forecast_200 <- forecast(
  fit_auto200,
  h = 5
)

forecast_211 <- forecast(
  fit_211,
  h = 5
)

cat("\n============================================================\n")
cat("ARIMA FORECASTS 2019-2023\n")
cat("============================================================\n")

print(forecast_013)
print(forecast_200)
print(forecast_211)

cat("\nFull-sample Kalman-smoothed trend, 2019-2023:\n")
print(
  tail(
    ts_full_kf_manual2,
    5
  )
)

# =============================================================================
# FORECAST-ERROR MEASURES
# =============================================================================

true_values <- c(
  200504.8,
  200829.0,
  203835.7,
  202202.7,
  200731.6
)

forecast_values_013 <- c(
  196200.3,
  197461.9,
  197397.9,
  197397.9,
  197397.9
)

forecast_values_200 <- c(
  197832.8,
  201604.5,
  204749.2,
  207073.1,
  208527.1
)

forecast_values_211 <- c(
  195624.6,
  196906.4,
  197568.5,
  198008.0,
  198245.1
)

calculate_errors <- function(true, predicted) {
  
  errors <- true - predicted
  
  rmse_value <- Metrics::rmse(true, predicted)
  mae_value  <- Metrics::mae(true, predicted)
  mape_value <- Metrics::mape(true, predicted) * 100
  fev_value  <- var(errors)
  
  return(c(
    MAE = mae_value,
    RMSE = rmse_value,
    MAPE = mape_value,
    FEV = fev_value
  ))
}

errors_013 <- calculate_errors(
  true_values,
  forecast_values_013
)

errors_200 <- calculate_errors(
  true_values,
  forecast_values_200
)

errors_211 <- calculate_errors(
  true_values,
  forecast_values_211
)

results_table <- rbind(
  "ARIMA(0,1,3)" = errors_013,
  "ARIMA(2,0,0)" = errors_200,
  "ARIMA(2,1,1)" = errors_211
)

print(round(results_table, 2))

# =============================================================================
# ARIMA(2,1,1) FITTED VALUES AND FORECAST PLOT
# =============================================================================

arima_fitted211 <- fitted(
  fit_211
)

data_plot <- data.frame(
  Year = c(
    time(ts_full_kf_manual2),
    time(arima_fitted211),
    time(forecast_211$mean)
  ),
  Value = c(
    as.numeric(ts_full_kf_manual2),
    as.numeric(arima_fitted211),
    as.numeric(forecast_211$mean)
  ),
  Series = c(
    rep(
      "Kalman-Smoothed Trend",
      length(ts_full_kf_manual2)
    ),
    rep(
      "ARIMA(2,1,1) Fitted",
      length(arima_fitted211)
    ),
    rep(
      "ARIMA(2,1,1) Forecast",
      length(forecast_211$mean)
    )
  )
)

p_arima_211 <- ggplot(
  data_plot,
  aes(
    x = Year,
    y = Value,
    colour = Series,
    linetype = Series
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  scale_colour_manual(
    name = "Series",
    values = c(
      "Kalman-Smoothed Trend" = "#1F78B4",
      "ARIMA(2,1,1) Fitted" = "grey25",
      "ARIMA(2,1,1) Forecast" = "#D95F02"
    ),
    breaks = c(
      "Kalman-Smoothed Trend",
      "ARIMA(2,1,1) Fitted",
      "ARIMA(2,1,1) Forecast"
    )
  ) +
  scale_linetype_manual(
    name = "Series",
    values = c(
      "Kalman-Smoothed Trend" = "solid",
      "ARIMA(2,1,1) Fitted" = "dotted",
      "ARIMA(2,1,1) Forecast" = "dashed"
    ),
    breaks = c(
      "Kalman-Smoothed Trend",
      "ARIMA(2,1,1) Fitted",
      "ARIMA(2,1,1) Forecast"
    )
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        linewidth = 1.1
      )
    ),
    linetype = guide_legend(
      override.aes = list(
        linewidth = 1.1
      )
    )
  ) +
  labs(
    x = "Year",
    y = "Total freight activity (million tonne-km)"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    axis.title = element_text(
      size = 12
    ),
    axis.text = element_text(
      size = 10
    ),
    legend.position = "top",
    legend.title = element_text(
      size = 10
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.key.width = unit(
      0.8,
      "cm"
    ),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(
  p_arima_211
)

# End of script
