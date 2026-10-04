# =============================================================================
# 08_gdp_arimax_benchmark.R
#
# GDP-augmented ARIMAX robustness analysis for annual freight activity.
#
# Input:
#   data/Tkm_Gdp.xlsx
# =============================================================================


# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "tseries",
  "forecast",
  "urca",
  "FinTS",
  "strucchange",
  "ggplot2",
  "dplyr",
  "gridExtra",
  "moments",
  "zoo",
  "vars",
  "aTSA"
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
library(forecast)
library(urca)
library(FinTS)
library(strucchange)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(moments)
library(zoo)
library(vars)


# ==============================================================================
# 1. TKM AND GDP SERIES
# ==============================================================================

# Raw series

ggplot(Tkm_Gdp, aes(x = Year)) +
  geom_line(aes(y = Tkm, colour = "Tkm")) +
  geom_line(aes(y = GdP, colour = "GDP")) +
  scale_colour_manual(values = c("Tkm" = "steelblue", "GDP" = "darkred")) +
  labs(
    title = "Freight Activity and GDP - Italy, 1990-2023",
    x = "Year",
    y = "Value",
    colour = "Series"
  ) +
  theme_minimal()

# Indexed series: 1990 = 100

Tkm_Gdp_norm <- Tkm_Gdp %>%
  mutate(
    Tkm_index = Tkm / Tkm[Year == 1990] * 100,
    GdP_index = GdP / GdP[Year == 1990] * 100
  )

ggplot(Tkm_Gdp_norm, aes(x = Year)) +
  geom_line(aes(y = Tkm_index, colour = "Tkm")) +
  geom_line(aes(y = GdP_index, colour = "GDP")) +
  scale_colour_manual(values = c("Tkm" = "steelblue", "GDP" = "darkred")) +
  labs(
    title = "Freight Activity and GDP - Index, 1990 = 100",
    x = "Year",
    y = "Index (1990 = 100)",
    colour = "Series"
  ) +
  theme_minimal()


# ==============================================================================
# 2. FREIGHT-ACTIVITY DIAGNOSTICS - TRAINING SAMPLE 1990-2018
# ==============================================================================

# Training series

ts_tkm_train <- ts(
  Tkm_Gdp$Tkm,
  start = 1990,
  end = 2018,
  frequency = 1
)

ts_gdp_train <- ts(
  Tkm_Gdp$GdP,
  start = 1990,
  end = 2018,
  frequency = 1
)

# ------------------------------------------------------------------------------
# 2.1 Descriptive statistics
# ------------------------------------------------------------------------------

summary_stats_tkm_train <- data.frame(
  Mean = mean(ts_tkm_train),
  Median = median(ts_tkm_train),
  Min = min(ts_tkm_train),
  Max = max(ts_tkm_train),
  Variance = var(ts_tkm_train),
  Standard_Deviation = sd(ts_tkm_train),
  Skewness = skewness(ts_tkm_train),
  Kurtosis = kurtosis(ts_tkm_train)
)

print(summary_stats_tkm_train)
summary(ts_tkm_train)

# ------------------------------------------------------------------------------
# 2.2 Distributional diagnostics
# ------------------------------------------------------------------------------

ggplot(Tkm_Gdp, aes(x = Tkm)) +
  geom_histogram(color = "black", fill = "lightblue", bins = 10) +
  labs(
    title = "Histogram - Freight Activity",
    x = "Tonne-km",
    y = "Frequency"
  ) +
  theme_minimal()

qqnorm(ts_tkm_train, main = "Q-Q Plot - Freight Activity")
qqline(ts_tkm_train, col = "red", lwd = 2)

shapiro.test(ts_tkm_train)
ks.test(scale(ts_tkm_train), "pnorm")

# ------------------------------------------------------------------------------
# 2.3 ACF/PACF and differences
# ------------------------------------------------------------------------------

acf_plot <- ggAcf(ts_tkm_train) +
  ggtitle("ACF - Freight Activity") +
  theme_minimal()

pacf_plot <- ggPacf(ts_tkm_train) +
  ggtitle("PACF - Freight Activity") +
  theme_minimal()

grid.arrange(acf_plot, pacf_plot, ncol = 2)

diff1_tkm_train <- diff(ts_tkm_train, differences = 1)
diff2_tkm_train <- diff(ts_tkm_train, differences = 2)

grid.arrange(
  ggAcf(diff1_tkm_train) + ggtitle("ACF - Freight Activity, d = 1"),
  ggPacf(diff1_tkm_train) + ggtitle("PACF - Freight Activity, d = 1"),
  ncol = 2
)

grid.arrange(
  ggAcf(diff2_tkm_train) + ggtitle("ACF - Freight Activity, d = 2"),
  ggPacf(diff2_tkm_train) + ggtitle("PACF - Freight Activity, d = 2"),
  ncol = 2
)

# ------------------------------------------------------------------------------
# 2.4 Stationarity tests
# ------------------------------------------------------------------------------

adf.test(ts_tkm_train)
adf.test(diff1_tkm_train)
adf.test(diff2_tkm_train)

kpss.test(ts_tkm_train)
kpss.test(diff1_tkm_train)
kpss.test(diff2_tkm_train)

# ------------------------------------------------------------------------------
# 2.5 Serial dependence and ARCH diagnostics
# ------------------------------------------------------------------------------

Box.test(ts_tkm_train, lag = 10, type = "Ljung-Box")
Box.test(diff1_tkm_train, lag = 10, type = "Ljung-Box")

ArchTest(ts_tkm_train, lags = 5)
ArchTest(diff1_tkm_train, lags = 5)

# ------------------------------------------------------------------------------
# 2.6 Structural diagnostics
# ------------------------------------------------------------------------------

breaks_tkm_train <- breakpoints(ts_tkm_train ~ 1)
plot(breaks_tkm_train)
summary(breaks_tkm_train)

roll_slope_tkm_train <- rollapply(
  ts_tkm_train,
  width = 5,
  by = 1,
  FUN = function(x) coef(lm(x ~ seq_along(x)))[2],
  align = "right",
  fill = NA
)

plot(
  roll_slope_tkm_train,
  type = "l",
  main = "Rolling Slope - Freight Activity",
  col = "darkgreen"
)
abline(h = 0, lty = 2)

spectrum(ts_tkm_train, main = "Spectral Density - Freight Activity")

cusum_tkm_train <- efp(ts_tkm_train ~ 1, type = "Rec-CUSUM")
plot(cusum_tkm_train)


# ==============================================================================
# 3. GDP DIAGNOSTICS - TRAINING SAMPLE 1990-2018
# ==============================================================================

# Differences

diff1_gdp_train <- diff(ts_gdp_train, differences = 1)
diff2_gdp_train <- diff(ts_gdp_train, differences = 2)

# ------------------------------------------------------------------------------
# 3.1 Descriptive statistics
# ------------------------------------------------------------------------------

summary_stats_gdp_train <- data.frame(
  Mean = mean(ts_gdp_train),
  Median = median(ts_gdp_train),
  Min = min(ts_gdp_train),
  Max = max(ts_gdp_train),
  Variance = var(ts_gdp_train),
  Standard_Deviation = sd(ts_gdp_train),
  Skewness = skewness(ts_gdp_train),
  Kurtosis = kurtosis(ts_gdp_train)
)

print(summary_stats_gdp_train)
summary(ts_gdp_train)

# ------------------------------------------------------------------------------
# 3.2 Distributional diagnostics
# ------------------------------------------------------------------------------

ggplot(data.frame(GdP = as.numeric(ts_gdp_train)), aes(x = GdP)) +
  geom_histogram(color = "black", fill = "lightblue", bins = 10) +
  labs(
    title = "Histogram - GDP, 1990-2018",
    x = "GDP",
    y = "Frequency"
  ) +
  theme_minimal()

qqnorm(ts_gdp_train, main = "Q-Q Plot - GDP")
qqline(ts_gdp_train, col = "red", lwd = 2)

shapiro.test(ts_gdp_train)
ks.test(scale(ts_gdp_train), "pnorm")

# ------------------------------------------------------------------------------
# 3.3 ACF/PACF and differences
# ------------------------------------------------------------------------------

grid.arrange(
  ggAcf(ts_gdp_train) + ggtitle("ACF - GDP, levels") + theme_minimal(),
  ggPacf(ts_gdp_train) + ggtitle("PACF - GDP, levels") + theme_minimal(),
  ncol = 2
)

grid.arrange(
  ggAcf(diff1_gdp_train) + ggtitle("ACF - GDP, d = 1") + theme_minimal(),
  ggPacf(diff1_gdp_train) + ggtitle("PACF - GDP, d = 1") + theme_minimal(),
  ncol = 2
)

grid.arrange(
  ggAcf(diff2_gdp_train) + ggtitle("ACF - GDP, d = 2") + theme_minimal(),
  ggPacf(diff2_gdp_train) + ggtitle("PACF - GDP, d = 2") + theme_minimal(),
  ncol = 2
)

# ------------------------------------------------------------------------------
# 3.4 Stationarity tests
# ------------------------------------------------------------------------------

adf.test(ts_gdp_train)
adf.test(diff1_gdp_train)
adf.test(diff2_gdp_train)

kpss.test(ts_gdp_train)
kpss.test(diff1_gdp_train)
kpss.test(diff2_gdp_train)

# ------------------------------------------------------------------------------
# 3.5 Serial dependence and ARCH diagnostics
# ------------------------------------------------------------------------------

Box.test(ts_gdp_train, lag = 10, type = "Ljung-Box")
Box.test(diff1_gdp_train, lag = 10, type = "Ljung-Box")

ArchTest(ts_gdp_train, lags = 5)
ArchTest(diff1_gdp_train, lags = 5)

# ------------------------------------------------------------------------------
# 3.6 Structural diagnostics
# ------------------------------------------------------------------------------

breaks_gdp_train <- breakpoints(ts_gdp_train ~ 1)
plot(breaks_gdp_train)
summary(breaks_gdp_train)

roll_slope_gdp_train <- rollapply(
  ts_gdp_train,
  width = 5,
  by = 1,
  FUN = function(x) coef(lm(x ~ seq_along(x)))[2],
  align = "right",
  fill = NA
)

plot(
  roll_slope_gdp_train,
  type = "l",
  main = "Rolling Slope - GDP",
  col = "darkgreen"
)
abline(h = 0, lty = 2)

spectrum(ts_gdp_train, main = "Spectral Density - GDP")

cusum_gdp_train <- efp(ts_gdp_train ~ 1, type = "Rec-CUSUM")
plot(cusum_gdp_train)


# ==============================================================================
# 4. TKM-GDP RELATIONSHIP
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1 Correlations in levels and first differences
# ------------------------------------------------------------------------------

cor_levels <- cor(
  ts_tkm_train,
  ts_gdp_train,
  use = "complete.obs"
)

diff1_tkm_train <- diff(ts_tkm_train)
diff1_gdp_train <- diff(ts_gdp_train)

cor_diff1 <- cor(
  diff1_tkm_train,
  diff1_gdp_train,
  use = "complete.obs"
)

print(
  list(
    correlation_levels = cor_levels,
    correlation_diff1 = cor_diff1
  )
)

# ------------------------------------------------------------------------------
# 4.2 Prewhitened cross-correlation function
# ------------------------------------------------------------------------------

# GDP is used as the prewhitened input series.
# Positive lags indicate GDP changes leading freight-activity changes.

dx <- diff1_gdp_train
dy <- diff1_tkm_train

mx <- auto.arima(
  dx,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE
)

rx <- residuals(mx)
ry <- residuals(Arima(dy, model = mx))

ccf(
  rx,
  ry,
  lag.max = 6,
  main = "Prewhitened CCF: GDP Changes -> Freight-Activity Changes"
)

# ------------------------------------------------------------------------------
# 4.3 Cointegration tests
# ------------------------------------------------------------------------------

# Engle-Granger cointegration test

eg_test <- aTSA::coint.test(
  y = as.numeric(ts_tkm_train),
  X = as.numeric(ts_gdp_train),
  d = 0,
  nlag = 1,
  output = TRUE
)

print(eg_test)

# Phillips-Ouliaris

po_test <- ca.po(
  cbind(ts_tkm_train, ts_gdp_train),
  demean = "constant",
  type = "Pz"
)
summary(po_test)

# Johansen trace and maximum-eigenvalue tests

sel <- VARselect(
  cbind(ts_tkm_train, ts_gdp_train),
  lag.max = 5,
  type = "const"
)$selection

k <- as.integer(sel["AIC(n)"])
if (is.na(k)) k <- 2
Kj <- max(2, k)

joh_trace <- ca.jo(
  cbind(ts_tkm_train, ts_gdp_train),
  type = "trace",
  ecdet = "const",
  K = Kj
)

joh_eigen <- ca.jo(
  cbind(ts_tkm_train, ts_gdp_train),
  type = "eigen",
  ecdet = "const",
  K = Kj
)

summary(joh_trace)
summary(joh_eigen)

# ------------------------------------------------------------------------------
# 4.4 Granger causality on first differences
# ------------------------------------------------------------------------------

lag_sel_d <- VARselect(
  cbind(diff1_tkm_train, diff1_gdp_train),
  lag.max = 5,
  type = "const"
)$selection

p_d <- lag_sel_d["AIC(n)"]
if (is.na(p_d)) p_d <- 2

var_d <- VAR(
  cbind(diff1_tkm_train, diff1_gdp_train),
  p = as.integer(p_d),
  type = "const"
)

# GDP -> freight activity
causality(var_d, cause = "diff1_gdp_train")

# Freight activity -> GDP
causality(var_d, cause = "diff1_tkm_train")


# ==============================================================================
# 5. GDP-AUGMENTED ARIMAX SPECIFICATION SEARCH
# ==============================================================================

# First-difference specifications

fit_aicc <- auto.arima(
  diff1_tkm_train,
  xreg = diff1_gdp_train,
  seasonal = FALSE,
  ic = "aicc"
)

fit_bic <- auto.arima(
  diff1_tkm_train,
  xreg = diff1_gdp_train,
  seasonal = FALSE,
  ic = "bic"
)

fit_aic <- auto.arima(
  diff1_tkm_train,
  xreg = diff1_gdp_train,
  seasonal = FALSE,
  ic = "aic"
)

# Level specifications

fit_aicc1 <- auto.arima(
  ts_tkm_train,
  xreg = ts_gdp_train,
  seasonal = FALSE,
  ic = "aicc"
)

fit_bic1 <- auto.arima(
  ts_tkm_train,
  xreg = ts_gdp_train,
  seasonal = FALSE,
  ic = "bic"
)

fit_aic1 <- auto.arima(
  ts_tkm_train,
  xreg = ts_gdp_train,
  seasonal = FALSE,
  ic = "aic"
)


# ==============================================================================
# 6. GDP-ARIMAX BENCHMARK FORECAST - 2019-2023
# ==============================================================================

# Holdout series

holdout_idx <- Tkm_Gdp$Year >= 2019 & Tkm_Gdp$Year <= 2023

ts_tkm_holdout <- ts(
  Tkm_Gdp$Tkm[holdout_idx],
  start = 2019,
  frequency = 1
)

ts_gdp_holdout <- ts(
  Tkm_Gdp$GdP[holdout_idx],
  start = 2019,
  frequency = 1
)

# GDP-augmented ARIMAX(0,1,0)

fit_gdp_arimax_010 <- forecast::Arima(
  ts_tkm_train,
  order = c(0, 1, 0),
  xreg = ts_gdp_train
)

forecast_gdp_arimax_010 <- forecast::forecast(
  fit_gdp_arimax_010,
  xreg = ts_gdp_holdout,
  h = length(ts_gdp_holdout)
)

# Corresponding univariate ARIMA benchmark

fit_arima_100 <- forecast::auto.arima(
  ts_tkm_train,
  max.p = 7,
  max.q = 7,
  seasonal = FALSE,
  stepwise = FALSE,
  approximation = FALSE,
  ic = "aicc"
)

stopifnot(
  all(
    forecast::arimaorder(fit_arima_100)[1:3] == c(1, 0, 0)
  )
)

forecast_arima_100 <- forecast::forecast(
  fit_arima_100,
  h = length(ts_tkm_holdout),
  level = 95
)

# Forecast table

arimax_forecast_table <- data.frame(
  Year = 2019:2023,
  Observed = as.numeric(ts_tkm_holdout),
  GDP = as.numeric(ts_gdp_holdout),
  ARIMA_100 = as.numeric(forecast_arima_100$mean),
  GDP_ARIMAX_010 = as.numeric(forecast_gdp_arimax_010$mean)
)

# Point-error measures

mae <- function(e) mean(abs(e))
rmse <- function(e) sqrt(mean(e^2))
mape <- function(e, y) mean(abs(e / y)) * 100
forecast_error_variance <- function(e) var(e)

error_arima_100 <-
  arimax_forecast_table$Observed - arimax_forecast_table$ARIMA_100

error_gdp_arimax_010 <-
  arimax_forecast_table$Observed - arimax_forecast_table$GDP_ARIMAX_010

arimax_accuracy_table <- data.frame(
  Model = c(
    "ARIMA(1,0,0)",
    "GDP-ARIMAX(0,1,0)"
  ),
  MAE = c(
    mae(error_arima_100),
    mae(error_gdp_arimax_010)
  ),
  RMSE = c(
    rmse(error_arima_100),
    rmse(error_gdp_arimax_010)
  ),
  MAPE = c(
    mape(error_arima_100, arimax_forecast_table$Observed),
    mape(error_gdp_arimax_010, arimax_forecast_table$Observed)
  ),
  Forecast_Error_Variance = c(
    forecast_error_variance(error_arima_100),
    forecast_error_variance(error_gdp_arimax_010)
  )
)


# ==============================================================================
# 7. OUTPUT TO EXTRACT
# ==============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("OUTPUT TO EXTRACT\n")
cat("####################################################################\n")

cat("\n--- ARIMAX SPECIFICATION SEARCH: FIRST DIFFERENCES ---\n")

cat("\nAICc:\n")
print(fit_aicc)

cat("\nBIC:\n")
print(fit_bic)

cat("\nAIC:\n")
print(fit_aic)


cat("\n--- ARIMAX SPECIFICATION SEARCH: LEVELS ---\n")

cat("\nAICc:\n")
print(fit_aicc1)

cat("\nBIC:\n")
print(fit_bic1)

cat("\nAIC:\n")
print(fit_aic1)

cat("\n--- GDP-ARIMAX(0,1,0): MODEL ---\n")
print(fit_gdp_arimax_010)

cat("\n--- GDP-ARIMAX(0,1,0): FORECAST 2019-2023 ---\n")
print(
  arimax_forecast_table,
  row.names = FALSE,
  digits = 12
)

cat("\n--- ARIMA vs GDP-ARIMAX: FORECAST ACCURACY ---\n")
print(
  arimax_accuracy_table,
  row.names = FALSE,
  digits = 12
)

# End of script
