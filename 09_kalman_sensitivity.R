# =============================================================================
# 09_kalman_sensitivity.R
#
# Sensitivity analysis for the Kalman calibration: baseline vs constrained-ML Test C.
#
# Input:
#   data/Total_Tonkm.xlsx
# =============================================================================


# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "KFAS",
  "forecast",
  "MSwM",
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

library(KFAS)
library(forecast)
library(MSwM)
library(strucchange)


# ==============================================================================
# 1. DATA
# ==============================================================================

stopifnot(
  exists("Total_Tonkm"),
  all(c("Year", "Total") %in% names(Total_Tonkm)),
  identical(as.integer(Total_Tonkm$Year), 1990:2023)
)

ts_total <- ts(
  Total_Tonkm$Total,
  start = 1990,
  end = 2023,
  frequency = 1
)

raw_train <- window(
  ts_total,
  start = 1990,
  end = 2018
)

raw_full <- window(
  ts_total,
  start = 1990,
  end = 2023
)

raw_holdout <- window(
  ts_total,
  start = 2019,
  end = 2023
)


# ==============================================================================
# 2. BASELINE KALMAN CALIBRATION
# ==============================================================================

build_baseline_kf <- function(
  y,
  q_slope,
  initial_level,
  initial_slope,
  p1_level,
  p1_slope
) {

  model <- SSModel(
    y ~ SSMtrend(
      degree = 2,
      Q = list(
        matrix(as.double(9999999)),
        matrix(as.double(q_slope))
      )
    ),
    H = array(
      as.double(9999999),
      dim = c(1, 1, 1)
    )
  )

  model$a1 <- matrix(
    as.double(c(initial_level, initial_slope)),
    ncol = 1
  )

  model$P1 <- diag(
    as.double(c(p1_level, p1_slope))
  )

  model$P1inf <- matrix(
    0,
    2,
    2
  )

  kfs <- KFS(
    model,
    filtering = "state",
    smoothing = "state"
  )

  trend <- ts(
    kfs$alphahat[, 1],
    start = start(y)[1],
    frequency = 1
  )

  list(
    model = model,
    kfs = kfs,
    trend = trend
  )
}


baseline_train_fit <- build_baseline_kf(
  y = raw_train,
  q_slope = 900000,
  initial_level = 191099.7,
  initial_slope = 450,
  p1_level = 200,
  p1_slope = 10
)

baseline_full_fit <- build_baseline_kf(
  y = raw_full,
  q_slope = 500000,
  initial_level = 189781.9,
  initial_slope = 350,
  p1_level = 500,
  p1_slope = 10
)

baseline_train <- baseline_train_fit$trend
baseline_full <- baseline_full_fit$trend


# ==============================================================================
# 3. TEST C: FIXED H, Q11 AND Q22 ESTIMATED BY ML
# ==============================================================================

fit_testC_kf <- function(y) {

  model <- SSModel(
    y ~ SSMtrend(
      degree = 2,
      Q = list(NA, NA)
    ),
    H = 9999999
  )

  fit <- fitSSM(
    model,
    inits = c(0.1, 0.1),
    method = "Nelder-Mead"
  )

  kfs <- KFS(
    fit$model,
    filtering = "state",
    smoothing = "state"
  )

  trend <- ts(
    kfs$alphahat[, 1],
    start = start(y)[1],
    frequency = 1
  )

  Q_mat <- kfs$model$Q[, , 1]
  H_mat <- kfs$model$H[, , 1]

  list(
    fit = fit,
    kfs = kfs,
    trend = trend,
    Q11 = Q_mat[1, 1],
    Q22 = Q_mat[2, 2],
    H = H_mat[1],
    convergence = fit$optim.out$convergence
  )
}


testC_train_fit <- fit_testC_kf(
  raw_train
)

testC_full_fit <- fit_testC_kf(
  raw_full
)

testC_train <- testC_train_fit$trend
testC_full <- testC_full_fit$trend


# ==============================================================================
# 4. SENSITIVITY OF THE EXTRACTED TREND
# ==============================================================================

get_bic_break_years <- function(x) {

  bp <- strucchange::breakpoints(
    x ~ 1
  )

  # For a breakpointsfull object, this is the minimum-BIC partition.
  idx <- bp$breakpoints
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0L) {
    return(character(0))
  }

  as.character(
    as.integer(
      time(x)[idx]
    )
  )
}


trend_sensitivity_metrics <- function(
  raw,
  trend,
  calibration,
  sample_label
) {

  d_raw <- diff(raw)
  d_trend <- diff(trend)

  d2_raw <- diff(d_raw)
  d2_trend <- diff(d_trend)

  residuals <- raw - trend

  directional_agreement <- 100 * mean(
    sign(d_raw) == sign(d_trend)
  )

  retained_d2_variance <- 100 * (
    var(d2_trend) / var(d2_raw)
  )

  residual_variance <- var(
    residuals
  )

  # This sensitivity table uses lag = 7 for the residual Ljung-Box check.
  residual_ljung_box_p <- Box.test(
    residuals,
    lag = 7,
    type = "Ljung-Box"
  )$p.value

  break_years <- get_bic_break_years(
    trend
  )

  data.frame(
    Calibration = calibration,
    Sample = sample_label,
    Directional_agreement_pct = directional_agreement,
    Retained_D2_variance_pct = retained_d2_variance,
    Residual_variance = residual_variance,
    Residual_Ljung_Box_p = residual_ljung_box_p,
    BIC_break_years = paste(
      break_years,
      collapse = ", "
    ),
    check.names = FALSE
  )
}


table_trend_sensitivity <- rbind(
  trend_sensitivity_metrics(
    raw_train,
    baseline_train,
    "Baseline",
    "Training"
  ),
  trend_sensitivity_metrics(
    raw_train,
    testC_train,
    "Test C",
    "Training"
  ),
  trend_sensitivity_metrics(
    raw_full,
    baseline_full,
    "Baseline",
    "Full"
  ),
  trend_sensitivity_metrics(
    raw_full,
    testC_full,
    "Test C",
    "Full"
  )
)


table_trend_sensitivity_print <- data.frame(
  Calibration =
    table_trend_sensitivity$Calibration,
  Sample =
    table_trend_sensitivity$Sample,
  Directional_agreement_pct =
    round(
      table_trend_sensitivity$Directional_agreement_pct,
      1
    ),
  Retained_D2_variance_pct =
    round(
      table_trend_sensitivity$Retained_D2_variance_pct,
      1
    ),
  Residual_variance_million =
    round(
      table_trend_sensitivity$Residual_variance / 1e6,
      2
    ),
  Residual_Ljung_Box_p =
    round(
      table_trend_sensitivity$Residual_Ljung_Box_p,
      3
    ),
  BIC_break_years =
    table_trend_sensitivity$BIC_break_years,
  check.names = FALSE
)


raw_breaks <- data.frame(
  Sample = c(
    "Training",
    "Full"
  ),
  BIC_break_years = c(
    paste(
      get_bic_break_years(raw_train),
      collapse = ", "
    ),
    paste(
      get_bic_break_years(raw_full),
      collapse = ", "
    )
  )
)


# ==============================================================================
# 5. MSR(2,0): BASELINE vs TEST C
# ==============================================================================

msr_control <- list(
  parallel = FALSE,
  maxiter = 5000,
  tol = 1e-8,
  trace = TRUE
)


fit_msr_20 <- function(y) {

  base_lm <- lm(
    as.numeric(y) ~ 1
  )

  msmFit(
    base_lm,
    k = 2,
    p = 0,
    sw = c(TRUE, TRUE),
    control = msr_control
  )
}


msm_baseline_20 <- fit_msr_20(
  baseline_train
)

msm_testC_20 <- fit_msr_20(
  testC_train
)


extract_msr_structure <- function(
  model,
  years = 1990:2018
) {

  means <- as.numeric(
    model@Coef[1:2, 1]
  )

  low_index <- which.min(
    means
  )

  high_index <- which.max(
    means
  )

  P <- model@transMat

  probabilities <- model@Fit@smoProb

  if (nrow(probabilities) == length(years) + 1L) {
    probabilities <- probabilities[
      -1,
      ,
      drop = FALSE
    ]
  }

  if (nrow(probabilities) != length(years)) {
    stop(
      "MSR probability matrix is not aligned with 1990-2018."
    )
  }

  model_class <- apply(
    probabilities,
    1,
    which.max
  )

  aligned_class <- ifelse(
    model_class == low_index,
    "Low",
    "High"
  )

  change_flag <- c(
    FALSE,
    aligned_class[-1] !=
      aligned_class[-length(aligned_class)]
  )

  transition_years <- years[
    change_flag
  ]

  list(
    low_mean =
      means[low_index],
    high_mean =
      means[high_index],
    p_LL =
      P[low_index, low_index],
    p_HH =
      P[high_index, high_index],
    classification =
      aligned_class,
    transition_years =
      transition_years
  )
}


msr_baseline_structure <- extract_msr_structure(
  msm_baseline_20
)

msr_testC_structure <- extract_msr_structure(
  msm_testC_20
)


classification_same <-
  msr_baseline_structure$classification ==
  msr_testC_structure$classification

classification_agreement_n <- sum(
  classification_same
)

classification_agreement_total <- length(
  classification_same
)

classification_agreement_pct <-
  100 *
  classification_agreement_n /
  classification_agreement_total


table_msr_sensitivity <- data.frame(
  Calibration = c(
    "Baseline",
    "Test C"
  ),
  Low_regime_mean = c(
    msr_baseline_structure$low_mean,
    msr_testC_structure$low_mean
  ),
  High_regime_mean = c(
    msr_baseline_structure$high_mean,
    msr_testC_structure$high_mean
  ),
  pLL = c(
    msr_baseline_structure$p_LL,
    msr_testC_structure$p_LL
  ),
  pHH = c(
    msr_baseline_structure$p_HH,
    msr_testC_structure$p_HH
  ),
  Transition_years = c(
    paste(
      msr_baseline_structure$transition_years,
      collapse = ", "
    ),
    paste(
      msr_testC_structure$transition_years,
      collapse = ", "
    )
  )
)


table_msr_sensitivity_print <- table_msr_sensitivity

table_msr_sensitivity_print$Low_regime_mean <- round(
  table_msr_sensitivity_print$Low_regime_mean,
  1
)

table_msr_sensitivity_print$High_regime_mean <- round(
  table_msr_sensitivity_print$High_regime_mean,
  1
)

table_msr_sensitivity_print$pLL <- round(
  table_msr_sensitivity_print$pLL,
  3
)

table_msr_sensitivity_print$pHH <- round(
  table_msr_sensitivity_print$pHH,
  3
)


regime_classification_comparison <- data.frame(
  Year = 1990:2018,
  Baseline =
    msr_baseline_structure$classification,
  Test_C =
    msr_testC_structure$classification,
  Same =
    classification_same
)


# ==============================================================================
# 6. RICHER MSR SPECIFICATIONS
# ==============================================================================

fit_msr_safe <- function(
  y,
  k,
  p,
  sw,
  label
) {

  base_lm <- lm(
    as.numeric(y) ~ 1
  )

  fit <- tryCatch(
    msmFit(
      base_lm,
      k = k,
      p = p,
      sw = sw,
      control = msr_control
    ),
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(fit)) {
    return(
      list(
        label = label,
        fit = NULL,
        fit_ok = FALSE,
        summary_ok = FALSE,
        residual_sd_ok = FALSE,
        trans_ok = FALSE,
        valid_fit = FALSE
      )
    )
  }

  summary_ok <- tryCatch(
    {
      capture.output(
        summary(fit)
      )
      TRUE
    },
    error = function(e) {
      FALSE
    }
  )

  residual_sd_ok <- tryCatch(
    {
      vals <- suppressWarnings(
        as.numeric(fit@std)
      )

      length(vals) > 0L &&
        all(is.finite(vals))
    },
    error = function(e) {
      FALSE
    }
  )

  trans_ok <- tryCatch(
    {
      all(
        is.finite(
          as.numeric(fit@transMat)
        )
      )
    },
    error = function(e) {
      FALSE
    }
  )

  valid_fit <-
    isTRUE(summary_ok) &&
    isTRUE(residual_sd_ok) &&
    isTRUE(trans_ok)

  list(
    label = label,
    fit = fit,
    fit_ok = TRUE,
    summary_ok = summary_ok,
    residual_sd_ok = residual_sd_ok,
    trans_ok = trans_ok,
    valid_fit = valid_fit
  )
}


msr_specs <- list(
  list(
    name = "MSR(2,0)",
    k = 2,
    p = 0,
    sw = c(TRUE, TRUE)
  ),
  list(
    name = "MSR(2,1)",
    k = 2,
    p = 1,
    sw = c(TRUE, TRUE, TRUE)
  ),
  list(
    name = "MSR(3,0)",
    k = 3,
    p = 0,
    sw = c(TRUE, TRUE)
  ),
  list(
    name = "MSR(3,1)",
    k = 3,
    p = 1,
    sw = c(TRUE, TRUE, TRUE)
  )
)


fit_spec_set <- function(
  y,
  calibration
) {

  results <- lapply(
    msr_specs,
    function(spec) {

      fit_msr_safe(
        y = y,
        k = spec$k,
        p = spec$p,
        sw = spec$sw,
        label = paste(
          calibration,
          spec$name
        )
      )
    }
  )

  data.frame(
    Calibration = calibration,
    Specification = vapply(
      msr_specs,
      function(x) x$name,
      character(1)
    ),
    Fit_returned = vapply(
      results,
      function(x) isTRUE(x$fit_ok),
      logical(1)
    ),
    Summary_available = vapply(
      results,
      function(x) isTRUE(x$summary_ok),
      logical(1)
    ),
    Finite_residual_SDs = vapply(
      results,
      function(x) isTRUE(x$residual_sd_ok),
      logical(1)
    ),
    Transition_matrix_OK = vapply(
      results,
      function(x) isTRUE(x$trans_ok),
      logical(1)
    ),
    Valid_fit = vapply(
      results,
      function(x) isTRUE(x$valid_fit),
      logical(1)
    )
  )
}


msr_status_baseline <- fit_spec_set(
  baseline_train,
  "Baseline"
)

msr_status_testC <- fit_spec_set(
  testC_train,
  "Test C"
)

msr_specification_status <- rbind(
  msr_status_baseline,
  msr_status_testC
)


# ==============================================================================
# 7. ARIMA SELECTION ON TEST C TRAINING TREND
# ==============================================================================

testC_aic <- auto.arima(
  testC_train,
  ic = "aic",
  stepwise = FALSE,
  approximation = FALSE
)

testC_aicc <- auto.arima(
  testC_train,
  ic = "aicc",
  stepwise = FALSE,
  approximation = FALSE
)

testC_bic <- auto.arima(
  testC_train,
  ic = "bic",
  stepwise = FALSE,
  approximation = FALSE
)


testC_d1_aic <- auto.arima(
  testC_train,
  d = 1,
  ic = "aic",
  stepwise = FALSE,
  approximation = FALSE
)

testC_d1_aicc <- auto.arima(
  testC_train,
  d = 1,
  ic = "aicc",
  stepwise = FALSE,
  approximation = FALSE
)

testC_d1_bic <- auto.arima(
  testC_train,
  d = 1,
  ic = "bic",
  stepwise = FALSE,
  approximation = FALSE
)


model_name <- function(x) {
  paste0(
    "ARIMA(",
    paste(
      arimaorder(x)[1:3],
      collapse = ","
    ),
    ")"
  )
}


testC_selection <- data.frame(
  Search = c(
    "Automatic d - AIC",
    "Automatic d - AICc",
    "Automatic d - BIC",
    "Fixed d=1 - AIC",
    "Fixed d=1 - AICc",
    "Fixed d=1 - BIC"
  ),
  Model = c(
    model_name(testC_aic),
    model_name(testC_aicc),
    model_name(testC_bic),
    model_name(testC_d1_aic),
    model_name(testC_d1_aicc),
    model_name(testC_d1_bic)
  )
)


# ==============================================================================
# 8. RETAINED ARIMA SPECIFICATIONS
# ==============================================================================

# Baseline specifications retained in the main Kalman-ARIMA analysis.

fit_baseline_200 <- Arima(
  baseline_train,
  order = c(2, 0, 0)
)

fit_baseline_013 <- Arima(
  baseline_train,
  order = c(0, 1, 3)
)

fit_baseline_211 <- Arima(
  baseline_train,
  order = c(2, 1, 1)
)


# Test C selected specifications.

fit_testC_201 <- Arima(
  testC_train,
  order = c(2, 0, 1),
  include.mean = TRUE
)

fit_testC_210 <- Arima(
  testC_train,
  order = c(2, 1, 0)
)


# Forecasts are generated only from the 1990-2018 training information.

fc_baseline_200 <- forecast(
  fit_baseline_200,
  h = 5
)

fc_baseline_013 <- forecast(
  fit_baseline_013,
  h = 5
)

fc_baseline_211 <- forecast(
  fit_baseline_211,
  h = 5
)

fc_testC_201 <- forecast(
  fit_testC_201,
  h = 5
)

fc_testC_210 <- forecast(
  fit_testC_210,
  h = 5
)


# ==============================================================================
# 9. ARIMA SENSITIVITY AGAINST THE COMMON RAW 2019-2023 HOLDOUT
# ==============================================================================

forecast_metrics <- function(
  actual,
  forecast_values
) {

  predicted <- as.numeric(
    forecast_values
  )

  actual <- as.numeric(
    actual
  )

  error <- actual - predicted

  c(
    MAE =
      mean(
        abs(error)
      ),
    RMSE =
      sqrt(
        mean(error^2)
      ),
    MAPE =
      mean(
        abs(error / actual)
      ) * 100
  )
}


make_metric_row <- function(
  calibration,
  specification,
  actual,
  forecast_values
) {

  metrics <- forecast_metrics(
    actual,
    forecast_values
  )

  data.frame(
    Calibration = calibration,
    Specification = specification,
    MAE = unname(metrics["MAE"]),
    RMSE = unname(metrics["RMSE"]),
    MAPE = unname(metrics["MAPE"])
  )
}


table_arima_sensitivity <- rbind(
  make_metric_row(
    "Baseline",
    "ARIMA(2,0,0)",
    raw_holdout,
    fc_baseline_200$mean
  ),
  make_metric_row(
    "Baseline",
    "ARIMA(0,1,3)",
    raw_holdout,
    fc_baseline_013$mean
  ),
  make_metric_row(
    "Baseline",
    "ARIMA(2,1,1)",
    raw_holdout,
    fc_baseline_211$mean
  ),
  make_metric_row(
    "Test C",
    "ARIMA(2,0,1)",
    raw_holdout,
    fc_testC_201$mean
  ),
  make_metric_row(
    "Test C",
    "ARIMA(2,1,0)",
    raw_holdout,
    fc_testC_210$mean
  )
)


table_arima_sensitivity_print <- table_arima_sensitivity

table_arima_sensitivity_print$MAE <- round(
  table_arima_sensitivity_print$MAE,
  2
)

table_arima_sensitivity_print$RMSE <- round(
  table_arima_sensitivity_print$RMSE,
  2
)

table_arima_sensitivity_print$MAPE <- round(
  table_arima_sensitivity_print$MAPE,
  2
)


# ==============================================================================
# 10. DIRECTIONAL CHECK: 2019-2023
# ==============================================================================

testC_full_holdout <- window(
  testC_full,
  start = 2019,
  end = 2023
)

directional_check <- data.frame(
  Period = c(
    "2019-2020",
    "2020-2021",
    "2021-2022",
    "2022-2023"
  ),
  Raw_change =
    as.numeric(
      diff(raw_holdout)
    ),
  TestC_full_change =
    as.numeric(
      diff(testC_full_holdout)
    ),
  ARIMA_201_change =
    diff(
      as.numeric(
        fc_testC_201$mean
      )
    ),
  ARIMA_210_change =
    diff(
      as.numeric(
        fc_testC_210$mean
      )
    )
)

directional_check$Raw_direction <- sign(
  directional_check$Raw_change
)

directional_check$TestC_full_direction <- sign(
  directional_check$TestC_full_change
)

directional_check$ARIMA_201_direction <- sign(
  directional_check$ARIMA_201_change
)

directional_check$ARIMA_210_direction <- sign(
  directional_check$ARIMA_210_change
)

directional_check$ARIMA_201_matches_raw <-
  directional_check$ARIMA_201_direction ==
  directional_check$Raw_direction

directional_check$ARIMA_210_matches_raw <-
  directional_check$ARIMA_210_direction ==
  directional_check$Raw_direction

directional_check$ARIMA_201_matches_TestC <-
  directional_check$ARIMA_201_direction ==
  directional_check$TestC_full_direction

directional_check$ARIMA_210_matches_TestC <-
  directional_check$ARIMA_210_direction ==
  directional_check$TestC_full_direction


# ==============================================================================
# 11. OUTPUT TO EXTRACT
# ==============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("OUTPUT TO EXTRACT\n")
cat("####################################################################\n")


cat("\n--- A. TEST C KALMAN ESTIMATES ---\n")

cat("\nTraining 1990-2018:\n")
cat("Q11 =", testC_train_fit$Q11, "\n")
cat("Q22 =", testC_train_fit$Q22, "\n")
cat("H   =", testC_train_fit$H, "\n")
cat("Convergence code =", testC_train_fit$convergence, "\n")

cat("\nFull 1990-2023:\n")
cat("Q11 =", testC_full_fit$Q11, "\n")
cat("Q22 =", testC_full_fit$Q22, "\n")
cat("H   =", testC_full_fit$H, "\n")
cat("Convergence code =", testC_full_fit$convergence, "\n")


cat("\n--- B. SENSITIVITY OF THE EXTRACTED TREND ---\n")

print(
  table_trend_sensitivity_print,
  row.names = FALSE
)

cat("\nRaw-series BIC break years used as reference:\n")

print(
  raw_breaks,
  row.names = FALSE
)


cat("\n--- C. MSR(2,0) SENSITIVITY ---\n")

print(
  table_msr_sensitivity_print,
  row.names = FALSE
)

cat(
  "\nRegime-classification agreement: ",
  classification_agreement_n,
  "/",
  classification_agreement_total,
  sprintf(
    " (%.1f%%)\n",
    classification_agreement_pct
  ),
  sep = ""
)

cat("\nAnnual aligned classification:\n")

print(
  regime_classification_comparison,
  row.names = FALSE
)


cat("\n--- D. RICHER MSR SPECIFICATION STATUS ---\n")

print(
  msr_specification_status,
  row.names = FALSE
)


cat("\n--- E. TEST C ARIMA SELECTION ---\n")

print(
  testC_selection,
  row.names = FALSE
)

cat("\nRetained Test C ARIMA(2,0,1):\n")
print(
  fit_testC_201
)

cat("\nRetained Test C ARIMA(2,1,0):\n")
print(
  fit_testC_210
)


cat("\n--- F. ARIMA SENSITIVITY AGAINST RAW 2019-2023 HOLDOUT ---\n")

print(
  table_arima_sensitivity_print,
  row.names = FALSE
)


cat("\n--- G. DIRECTIONAL CHECK 2019-2023 ---\n")

print(
  directional_check,
  row.names = FALSE,
  digits = 12
)


cat("\n--- H. FORECAST TRAJECTORIES ---\n")

forecast_trajectories <- data.frame(
  Year = 2019:2023,
  Raw = as.numeric(raw_holdout),
  Baseline_ARIMA_200 =
    as.numeric(fc_baseline_200$mean),
  Baseline_ARIMA_013 =
    as.numeric(fc_baseline_013$mean),
  Baseline_ARIMA_211 =
    as.numeric(fc_baseline_211$mean),
  TestC_ARIMA_201 =
    as.numeric(fc_testC_201$mean),
  TestC_ARIMA_210 =
    as.numeric(fc_testC_210$mean)
)

print(
  forecast_trajectories,
  row.names = FALSE,
  digits = 12
)


# End of script