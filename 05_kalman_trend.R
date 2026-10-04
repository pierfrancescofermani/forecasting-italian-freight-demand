# =============================================================================
# 05_kalman_trend.R
#
# Kalman local-linear-trend extraction and calibration diagnostics.
#
# Input:
#   data/Total_Tonkm.xlsx
# =============================================================================

# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "KFAS",
  "forecast",
  "tseries",
  "FinTS",
  "strucchange",
  "ggplot2",
  "moments",
  "patchwork",
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

library(KFAS)
library(forecast)
library(tseries)
library(FinTS)
library(strucchange)
library(ggplot2)
library(moments)
library(patchwork)
library(gridExtra)

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
# CALIBRATION DIAGNOSTICS: LOCAL-LEVEL CONFIGURATIONS
# =============================================================================

get_break_years <- function(x) {
  bp <- strucchange::breakpoints(x ~ 1)
  idx <- bp$breakpoints
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0L) {
    return("")
  }

  paste(
    as.numeric(time(x)[idx]),
    collapse = ", "
  )
}

get_break_years_m <- function(x, m) {
  bp_full <- strucchange::breakpoints(x ~ 1)
  bp_selected <- strucchange::breakpoints(
    bp_full,
    breaks = m
  )

  idx <- bp_selected$breakpoints
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0L) {
    return("")
  }

  paste(
    as.numeric(time(x)[idx]),
    collapse = ", "
  )
}

# ---- A1 ----------------------------------------------------------------------

model_kf_d1_train <- SSModel(
  ts_kf_train ~ SSMtrend(
    degree = 1,
    Q = NA
  ),
  H = NA
)

fit_kf_d1_train <- fitSSM(
  model_kf_d1_train,
  inits = c(0.1, 0.1),
  method = "BFGS"
)

smoothed_kf_d1_train <- KFS(
  fit_kf_d1_train$model,
  filtering = "state",
  smoothing = "state"
)

ts_train_kf_d1 <- ts(
  smoothed_kf_d1_train$alphahat[, 1],
  start = 1990,
  frequency = 1
)

Q_A1 <- smoothed_kf_d1_train$model$Q[1, 1, 1]
H_A1 <- smoothed_kf_d1_train$model$H[1, 1, 1]
max_dev_A1 <- max(abs(ts_kf_train - ts_train_kf_d1), na.rm = TRUE)
se_mu_last_A1 <- sqrt(
  smoothed_kf_d1_train$V[1, 1, length(ts_kf_train)]
)
breaks_A1 <- get_break_years(ts_train_kf_d1)
breaks_A1_m1 <- get_break_years_m(ts_train_kf_d1, 1)
breaks_A1_m2 <- get_break_years_m(ts_train_kf_d1, 2)

# ---- A2 ----------------------------------------------------------------------

model_kf_d1_train_VAR <- SSModel(
  ts_kf_train ~ SSMtrend(
    degree = 1,
    Q = NA
  ),
  H = NA
)

fit_kf_d1_train_VAR <- try(
  fitSSM(
    model_kf_d1_train_VAR,
    inits = c(
      log(var(ts_kf_train)),
      log(var(ts_kf_train))
    ),
    method = "BFGS"
  ),
  silent = TRUE
)

if (inherits(fit_kf_d1_train_VAR, "try-error")) {
  A2_status <- "Fit not returned"
  Q_A2 <- NA_real_
  H_A2 <- NA_real_
  max_dev_A2 <- NA_real_
  breaks_A2 <- ""
  breaks_A2_m1 <- ""
  breaks_A2_m2 <- ""
} else {
  A2_status <- "Fit returned"

  smoothed_kf_d1_train_VAR <- KFS(
    fit_kf_d1_train_VAR$model,
    filtering = "state",
    smoothing = "state"
  )

  ts_train_kf_d1_A2 <- ts(
    smoothed_kf_d1_train_VAR$alphahat[, 1],
    start = 1990,
    frequency = 1
  )

  Q_A2 <- smoothed_kf_d1_train_VAR$model$Q[1, 1, 1]
  H_A2 <- smoothed_kf_d1_train_VAR$model$H[1, 1, 1]
  max_dev_A2 <- max(abs(ts_kf_train - ts_train_kf_d1_A2), na.rm = TRUE)
  breaks_A2 <- get_break_years(ts_train_kf_d1_A2)
  breaks_A2_m1 <- get_break_years_m(ts_train_kf_d1_A2, 1)
  breaks_A2_m2 <- get_break_years_m(ts_train_kf_d1_A2, 2)
}

# ---- A3 ----------------------------------------------------------------------

scale_fac_A3 <- 1e7
y_A3 <- ts_kf_train / scale_fac_A3

model_kf_d1_A3 <- SSModel(
  y_A3 ~ SSMtrend(
    degree = 1,
    Q = NA
  ),
  H = NA
)

init_A3 <- log(
  rep(
    var(y_A3, na.rm = TRUE),
    2
  )
)

fit_kf_d1_A3 <- fitSSM(
  model_kf_d1_A3,
  inits = init_A3,
  method = "BFGS"
)

kfs_kf_d1_A3 <- KFS(
  fit_kf_d1_A3$model,
  filtering = "state",
  smoothing = "state"
)

ts_train_kf_d1_A3 <- ts(
  drop(kfs_kf_d1_A3$alphahat) * scale_fac_A3,
  start = 1990,
  frequency = 1
)

Q_A3_scaled <- kfs_kf_d1_A3$model$Q[1, 1, 1]
H_A3_scaled <- kfs_kf_d1_A3$model$H[1, 1, 1]
max_dev_A3 <- max(abs(ts_kf_train - ts_train_kf_d1_A3), na.rm = TRUE)
breaks_A3 <- get_break_years(ts_train_kf_d1_A3)
breaks_A3_m1 <- get_break_years_m(ts_train_kf_d1_A3, 1)
breaks_A3_m2 <- get_break_years_m(ts_train_kf_d1_A3, 2)

# ---- A4 ----------------------------------------------------------------------

v_A4 <- var(ts_kf_train, na.rm = TRUE)
limit_A4 <- 1e7
scale_fac_A4 <- sqrt(
  v_A4 / (0.9 * limit_A4)
)

y_A4 <- ts_kf_train / scale_fac_A4

model_kf_d1_A4 <- SSModel(
  y_A4 ~ SSMtrend(
    degree = 1,
    Q = NA
  ),
  H = NA
)

init_A4 <- rep(
  log(var(y_A4, na.rm = TRUE)),
  2
)

fit_kf_d1_A4 <- fitSSM(
  model_kf_d1_A4,
  inits = init_A4,
  method = "BFGS"
)

kfs_kf_d1_A4 <- KFS(
  fit_kf_d1_A4$model,
  smoothing = "state"
)

ts_train_kf_d1_A4 <- ts(
  drop(kfs_kf_d1_A4$alphahat) * scale_fac_A4,
  start = 1990,
  frequency = 1
)

Q_A4_scaled <- kfs_kf_d1_A4$model$Q[1, 1, 1]
H_A4_scaled <- kfs_kf_d1_A4$model$H[1, 1, 1]
max_dev_A4 <- max(abs(ts_kf_train - ts_train_kf_d1_A4), na.rm = TRUE)
breaks_A4 <- get_break_years(ts_train_kf_d1_A4)
breaks_A4_m1 <- get_break_years_m(ts_train_kf_d1_A4, 1)
breaks_A4_m2 <- get_break_years_m(ts_train_kf_d1_A4, 2)

# ---- A5 ----------------------------------------------------------------------

v_A5 <- var(ts_kf_train, na.rm = TRUE)
limit_A5 <- 1e7
eps_A5 <- 1e-7

scale_fac_A5 <- sqrt(
  v_A5 / (limit_A5 * (1 - eps_A5))
)

y_A5 <- ts_kf_train / scale_fac_A5

model_kf_d1_A5 <- SSModel(
  y_A5 ~ SSMtrend(
    degree = 1,
    Q = NA
  ),
  H = NA
)

init_A5 <- rep(
  log(var(y_A5, na.rm = TRUE)),
  2
)

fit_kf_d1_A5 <- fitSSM(
  model_kf_d1_A5,
  inits = init_A5,
  method = "BFGS"
)

kfs_kf_d1_A5 <- KFS(
  fit_kf_d1_A5$model,
  smoothing = "state"
)

ts_train_kf_d1_A5 <- ts(
  drop(kfs_kf_d1_A5$alphahat) * scale_fac_A5,
  start = 1990,
  frequency = 1
)

Q_A5_scaled <- kfs_kf_d1_A5$model$Q[1, 1, 1]
H_A5_scaled <- kfs_kf_d1_A5$model$H[1, 1, 1]
max_dev_A5 <- max(abs(ts_kf_train - ts_train_kf_d1_A5), na.rm = TRUE)
breaks_A5 <- get_break_years(ts_train_kf_d1_A5)
breaks_A5_m1 <- get_break_years_m(ts_train_kf_d1_A5, 1)
breaks_A5_m2 <- get_break_years_m(ts_train_kf_d1_A5, 2)

local_level_diagnostics <- data.frame(
  Configuration = c("A1", "A2", "A3", "A4", "A5"),
  Status = c(
    "Fit returned",
    A2_status,
    "Fit returned",
    "Fit returned",
    "Fit returned"
  ),
  Q = c(
    Q_A1,
    Q_A2,
    Q_A3_scaled,
    Q_A4_scaled,
    Q_A5_scaled
  ),
  H = c(
    H_A1,
    H_A2,
    H_A3_scaled,
    H_A4_scaled,
    H_A5_scaled
  ),
  Max_abs_raw_minus_trend = c(
    max_dev_A1,
    max_dev_A2,
    max_dev_A3,
    max_dev_A4,
    max_dev_A5
  ),
  Break_years_BIC = c(
    breaks_A1,
    breaks_A2,
    breaks_A3,
    breaks_A4,
    breaks_A5
  ),
  Break_years_m1 = c(
    breaks_A1_m1,
    breaks_A2_m1,
    breaks_A3_m1,
    breaks_A4_m1,
    breaks_A5_m1
  ),
  Break_years_m2 = c(
    breaks_A1_m2,
    breaks_A2_m2,
    breaks_A3_m2,
    breaks_A4_m2,
    breaks_A5_m2
  )
)

print(local_level_diagnostics, row.names = FALSE)
cat("\nA1 final smoothed-level standard error:", se_mu_last_A1, "\n")

# =============================================================================
# CALIBRATION DIAGNOSTICS: LOCAL-LINEAR-TREND CONFIGURATIONS
# =============================================================================

# ---- B1 ----------------------------------------------------------------------

model_kf_d2_train <- SSModel(
  ts_kf_train ~ SSMtrend(
    degree = 2,
    Q = list(NA, NA)
  ),
  H = NA
)

fit_kf_d2_train <- fitSSM(
  model_kf_d2_train,
  inits = log(c(100, 100, 100)),
  method = "BFGS"
)

smoothed_kf_d2_train <- KFS(
  fit_kf_d2_train$model,
  filtering = "state",
  smoothing = "state"
)

ts_train_kf_d2_B1 <- ts(
  smoothed_kf_d2_train$alphahat[, 1],
  start = 1990,
  frequency = 1
)

Q_B1 <- smoothed_kf_d2_train$model$Q[, , 1]
H_B1 <- smoothed_kf_d2_train$model$H[, , 1]
max_dev_B1 <- max(abs(ts_kf_train - ts_train_kf_d2_B1), na.rm = TRUE)
breaks_B1 <- get_break_years(ts_train_kf_d2_B1)
breaks_B1_m1 <- get_break_years_m(ts_train_kf_d2_B1, 1)
breaks_B1_m2 <- get_break_years_m(ts_train_kf_d2_B1, 2)

# ---- B2 ----------------------------------------------------------------------

model_kf_d2_log_train <- SSModel(
  ts_kf_train ~ SSMtrend(
    degree = 2,
    Q = list(NA, NA)
  ),
  H = NA
)

fit_kf_d2_log_train <- try(
  fitSSM(
    model_kf_d2_log_train,
    inits = log(c(1000, 1000, 1000)),
    method = "BFGS"
  ),
  silent = TRUE
)

if (inherits(fit_kf_d2_log_train, "try-error")) {
  B2_status <- "Fit not returned"
  Q_B2 <- matrix(NA_real_, 2, 2)
  H_B2 <- matrix(NA_real_, 1, 1)
  max_dev_B2 <- NA_real_
  breaks_B2 <- ""
  breaks_B2_m1 <- ""
  breaks_B2_m2 <- ""
} else {
  B2_status <- "Fit returned"

  smoothed_kf_d2_log_train <- KFS(
    fit_kf_d2_log_train$model,
    filtering = "state",
    smoothing = "state"
  )

  ts_train_kf_d2_B2 <- ts(
    smoothed_kf_d2_log_train$alphahat[, 1],
    start = 1990,
    frequency = 1
  )

  Q_B2 <- smoothed_kf_d2_log_train$model$Q[, , 1]
  H_B2 <- smoothed_kf_d2_log_train$model$H[, , 1]
  max_dev_B2 <- max(abs(ts_kf_train - ts_train_kf_d2_B2), na.rm = TRUE)
  breaks_B2 <- get_break_years(ts_train_kf_d2_B2)
  breaks_B2_m1 <- get_break_years_m(ts_train_kf_d2_B2, 1)
  breaks_B2_m2 <- get_break_years_m(ts_train_kf_d2_B2, 2)
}

# ---- B3 ----------------------------------------------------------------------

model_kf_d2_train_b <- SSModel(
  ts_kf_train ~ SSMtrend(
    degree = 2,
    Q = list(NA, NA)
  ),
  H = NA
)

init_max_adm_b <- rep(
  log(1e7 - 1),
  3
)

fit_kf_d2_train_max_b <- fitSSM(
  model_kf_d2_train_b,
  inits = init_max_adm_b,
  method = "BFGS"
)

smoothed_kf_d2_train_max_b <- KFS(
  fit_kf_d2_train_max_b$model,
  filtering = "state",
  smoothing = "state"
)

ts_train_kf_d2_B3 <- ts(
  smoothed_kf_d2_train_max_b$alphahat[, 1],
  start = 1990,
  frequency = 1
)

Q_B3 <- smoothed_kf_d2_train_max_b$model$Q[, , 1]
H_B3 <- smoothed_kf_d2_train_max_b$model$H[, , 1]
max_dev_B3 <- max(abs(ts_kf_train - ts_train_kf_d2_B3), na.rm = TRUE)
breaks_B3 <- get_break_years(ts_train_kf_d2_B3)
breaks_B3_m1 <- get_break_years_m(ts_train_kf_d2_B3, 1)
breaks_B3_m2 <- get_break_years_m(ts_train_kf_d2_B3, 2)

# ---- B4 ----------------------------------------------------------------------

model_kf_d2_train_realinit <- SSModel(
  ts_kf_train ~ SSMtrend(
    degree = 2,
    Q = list(NA, NA)
  ),
  H = NA
)

var_train_d2_realinit <- var(
  ts_kf_train,
  na.rm = TRUE
)

v0_train_d2_realinit <- min(
  var_train_d2_realinit,
  1e7 - 1
)

init_d2_realinit_log <- rep(
  log(v0_train_d2_realinit),
  3
)

fit_kf_d2_train_realinit <- fitSSM(
  model_kf_d2_train_realinit,
  inits = init_d2_realinit_log,
  method = "BFGS"
)

kfs_kf_d2_train_realinit <- KFS(
  fit_kf_d2_train_realinit$model,
  filtering = "state",
  smoothing = "state"
)

ts_train_kf_d2_B4 <- ts(
  kfs_kf_d2_train_realinit$alphahat[, 1],
  start = 1990,
  frequency = 1
)

Q_B4 <- kfs_kf_d2_train_realinit$model$Q[, , 1]
H_B4 <- kfs_kf_d2_train_realinit$model$H[, , 1]
max_dev_B4 <- max(abs(ts_kf_train - ts_train_kf_d2_B4), na.rm = TRUE)
breaks_B4 <- get_break_years(ts_train_kf_d2_B4)
breaks_B4_m1 <- get_break_years_m(ts_train_kf_d2_B4, 1)
breaks_B4_m2 <- get_break_years_m(ts_train_kf_d2_B4, 2)

# ---- B5 ----------------------------------------------------------------------

scale_fac_B <- 1e7
y_B <- ts_kf_train / scale_fac_B

model_kf_d2_train_B <- SSModel(
  y_B ~ SSMtrend(
    degree = 2,
    Q = list(NA, NA)
  ),
  H = NA
)

init_B_log <- rep(
  log(var(y_B, na.rm = TRUE)),
  3
)

fit_kf_d2_train_B <- fitSSM(
  model_kf_d2_train_B,
  inits = init_B_log,
  method = "BFGS"
)

kfs_kf_d2_train_B <- KFS(
  fit_kf_d2_train_B$model,
  smoothing = "state"
)

ts_train_kf_d2_B5_scaled <- ts(
  kfs_kf_d2_train_B$alphahat[, 1],
  start = 1990,
  frequency = 1
)

ts_train_kf_d2_B5 <- ts_train_kf_d2_B5_scaled * scale_fac_B

Q_B5 <- kfs_kf_d2_train_B$model$Q[, , 1]
H_B5 <- kfs_kf_d2_train_B$model$H[, , 1]
max_dev_B5 <- max(abs(ts_kf_train - ts_train_kf_d2_B5), na.rm = TRUE)
breaks_B5 <- get_break_years(ts_train_kf_d2_B5)
breaks_B5_m1 <- get_break_years_m(ts_train_kf_d2_B5, 1)
breaks_B5_m2 <- get_break_years_m(ts_train_kf_d2_B5, 2)

local_linear_trend_diagnostics <- data.frame(
  Configuration = c("B1", "B2", "B3", "B4", "B5"),
  Status = c(
    "Fit returned",
    B2_status,
    "Fit returned",
    "Fit returned",
    "Fit returned"
  ),
  Q11 = c(
    Q_B1[1, 1],
    Q_B2[1, 1],
    Q_B3[1, 1],
    Q_B4[1, 1],
    Q_B5[1, 1]
  ),
  Q22 = c(
    Q_B1[2, 2],
    Q_B2[2, 2],
    Q_B3[2, 2],
    Q_B4[2, 2],
    Q_B5[2, 2]
  ),
  H = c(
    H_B1[1],
    H_B2[1],
    H_B3[1],
    H_B4[1],
    H_B5[1]
  ),
  Max_abs_raw_minus_trend = c(
    max_dev_B1,
    max_dev_B2,
    max_dev_B3,
    max_dev_B4,
    max_dev_B5
  ),
  Break_years_BIC = c(
    breaks_B1,
    breaks_B2,
    breaks_B3,
    breaks_B4,
    breaks_B5
  ),
  Break_years_m1 = c(
    breaks_B1_m1,
    breaks_B2_m1,
    breaks_B3_m1,
    breaks_B4_m1,
    breaks_B5_m1
  ),
  Break_years_m2 = c(
    breaks_B1_m2,
    breaks_B2_m2,
    breaks_B3_m2,
    breaks_B4_m2,
    breaks_B5_m2
  )
)

print(local_linear_trend_diagnostics, row.names = FALSE)


# =============================================================================
# BASELINE LOCAL-LINEAR-TREND EXTRACTION: TRAINING SAMPLE
# =============================================================================

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
  as.double(
    c(initial_level, initial_slope)
  ),
  ncol = 1
)

model_kf_train_manual2$P1 <- diag(
  as.double(
    c(P1_level, P1_slope)
  )
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
train_resid_kf_manual2 <- ts_kf_train - train_kf_manual2

ts_train_kf_manual2 <- ts(
  train_kf_manual2,
  start = 1990,
  frequency = 1
)

ts_train_resid_kf_manual2 <- ts(
  train_resid_kf_manual2,
  start = 1990,
  frequency = 1
)

cat(
  "\nBaseline training variances:\n",
  "Q11 =", smoothed_kf_train_manual2$model$Q[1, 1, 1], "\n",
  "Q22 =", smoothed_kf_train_manual2$model$Q[2, 2, 1], "\n",
  "H   =", smoothed_kf_train_manual2$model$H[1, 1, 1], "\n"
)

# =============================================================================
# BASELINE LOCAL-LINEAR-TREND EXTRACTION: FULL SAMPLE
# =============================================================================

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
  as.double(
    c(initial_level, initial_slope)
  ),
  ncol = 1
)

model_kf_full_manual2$P1 <- diag(
  as.double(
    c(P1_level, P1_slope)
  )
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
full_resid_kf_manual2 <- ts_kf_full - full_kf_manual2

ts_full_kf_manual2 <- ts(
  full_kf_manual2,
  start = 1990,
  frequency = 1
)

ts_full_resid_kf_manual2 <- ts(
  full_resid_kf_manual2,
  start = 1990,
  frequency = 1
)

cat(
  "\nBaseline full-sample variances:\n",
  "Q11 =", smoothed_kf_full_manual2$model$Q[1, 1, 1], "\n",
  "Q22 =", smoothed_kf_full_manual2$model$Q[2, 2, 1], "\n",
  "H   =", smoothed_kf_full_manual2$model$H[1, 1, 1], "\n"
)

# =============================================================================
# RAW SERIES AND BASELINE KALMAN TREND
# =============================================================================

df_combined_plot <- data.frame(
  Year = rep(1990:2018, 2),
  Value = c(
    as.numeric(ts_kf_train),
    as.numeric(ts_train_kf_manual2)
  ),
  Series = rep(
    c(
      "Raw series",
      "Kalman-smoothed trend"
    ),
    each = length(1990:2018)
  )
)

df_full_combined_plot <- data.frame(
  Year = rep(1990:2023, 2),
  Value = c(
    as.numeric(ts_kf_full),
    as.numeric(ts_full_kf_manual2)
  ),
  Series = rep(
    c(
      "Raw series",
      "Kalman-smoothed trend"
    ),
    each = length(1990:2023)
  )
)

p1 <- ggplot(
  df_combined_plot,
  aes(
    x = Year,
    y = Value,
    color = Series
  )
) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    name = "Series",
    values = c(
      "Raw series" = "grey25",
      "Kalman-smoothed trend" = "#1F78B4"
    )
  ) +
  labs(
    title = "Training sample",
    x = "Year",
    y = "Total freight activity (million tonne-km)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_blank()
  )

p2 <- ggplot(
  df_full_combined_plot,
  aes(
    x = Year,
    y = Value,
    color = Series
  )
) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    name = "Series",
    values = c(
      "Raw series" = "grey25",
      "Kalman-smoothed trend" = "#1F78B4"
    )
  ) +
  labs(
    title = "Full sample",
    x = "Year",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.subtitle = element_blank()
  )

kalman_trend_plot <- (
  p1 + p2 +
    plot_layout(
      guides = "collect",
      ncol = 2
    )
) &
  theme(
    legend.position = "top"
  )

print(kalman_trend_plot)

# =============================================================================
# BASELINE TREND DIAGNOSTICS
# =============================================================================

trend_diagnostics <- function(x, label) {

  descriptive <- data.frame(
    Series = label,
    Mean = mean(x),
    Median = median(x),
    Min = min(x),
    Max = max(x),
    Variance = var(x),
    Standard_Deviation = sd(x),
    Skewness = moments::skewness(as.numeric(x)),
    Kurtosis = moments::kurtosis(as.numeric(x))
  )

  shapiro_result <- shapiro.test(
    as.numeric(x)
  )

  ks_result <- suppressWarnings(
    ks.test(
      as.numeric(scale(x)),
      "pnorm"
    )
  )

  adf_result <- suppressWarnings(
    adf.test(x)
  )

  kpss_result <- suppressWarnings(
    kpss.test(x)
  )

  ljung_result <- Box.test(
    x,
    lag = 10,
    type = "Ljung-Box"
  )

  arch_result <- FinTS::ArchTest(
    as.numeric(x),
    lags = 5
  )

  bp_result <- strucchange::breakpoints(
    x ~ 1
  )

  bp_index <- bp_result$breakpoints
  bp_index <- bp_index[!is.na(bp_index)]

  bp_years <- if (length(bp_index) == 0L) {
    numeric(0)
  } else {
    as.numeric(
      time(x)[bp_index]
    )
  }

  cusum_result <- strucchange::efp(
    x ~ 1,
    type = "Rec-CUSUM"
  )

  cusum_test <- strucchange::sctest(
    cusum_result
  )

  cusum_values <- as.vector(
    cusum_result$process
  )

  cusum_bounds <- as.numeric(
    boundary(cusum_result)
  )

  cusum_df <- data.frame(
    Year = as.numeric(time(x)),
    CUSUM = cusum_values,
    Bound = cusum_bounds,
    Out_of_bounds = abs(cusum_values) > cusum_bounds
  )

  list(
    descriptive = descriptive,
    shapiro = shapiro_result,
    ks = ks_result,
    adf = adf_result,
    kpss = kpss_result,
    ljung_box = ljung_result,
    arch = arch_result,
    breakpoints = bp_result,
    break_years = bp_years,
    cusum = cusum_result,
    cusum_test = cusum_test,
    cusum_df = cusum_df
  )
}

diag_kf_train <- trend_diagnostics(
  ts_train_kf_manual2,
  "Training 1990-2018"
)

diag_kf_full <- trend_diagnostics(
  ts_full_kf_manual2,
  "Full 1990-2023"
)

print(diag_kf_train$descriptive, row.names = FALSE)
print(diag_kf_full$descriptive, row.names = FALSE)

print(diag_kf_train$shapiro)
print(diag_kf_train$ks)
print(diag_kf_train$adf)
print(diag_kf_train$kpss)
print(diag_kf_train$ljung_box)
print(diag_kf_train$arch)
print(summary(diag_kf_train$breakpoints))
print(diag_kf_train$cusum_test)
print(
  diag_kf_train$cusum_df[
    diag_kf_train$cusum_df$Out_of_bounds,
  ]
)

print(diag_kf_full$shapiro)
print(diag_kf_full$ks)
print(diag_kf_full$adf)
print(diag_kf_full$kpss)
print(diag_kf_full$ljung_box)
print(diag_kf_full$arch)
print(summary(diag_kf_full$breakpoints))
print(diag_kf_full$cusum_test)
print(
  diag_kf_full$cusum_df[
    diag_kf_full$cusum_df$Out_of_bounds,
  ]
)

trend_diagnostics_summary <- data.frame(
  Sample = c(
    "Training 1990-2018",
    "Full 1990-2023"
  ),
  Mean = c(
    mean(ts_train_kf_manual2),
    mean(ts_full_kf_manual2)
  ),
  SD = c(
    sd(ts_train_kf_manual2),
    sd(ts_full_kf_manual2)
  ),
  Skewness = c(
    moments::skewness(as.numeric(ts_train_kf_manual2)),
    moments::skewness(as.numeric(ts_full_kf_manual2))
  ),
  Kurtosis = c(
    moments::kurtosis(as.numeric(ts_train_kf_manual2)),
    moments::kurtosis(as.numeric(ts_full_kf_manual2))
  ),
  Shapiro_p = c(
    diag_kf_train$shapiro$p.value,
    diag_kf_full$shapiro$p.value
  ),
  KS_p = c(
    diag_kf_train$ks$p.value,
    diag_kf_full$ks$p.value
  ),
  ADF_p = c(
    diag_kf_train$adf$p.value,
    diag_kf_full$adf$p.value
  ),
  KPSS_p = c(
    diag_kf_train$kpss$p.value,
    diag_kf_full$kpss$p.value
  ),
  Ljung_Box_p = c(
    diag_kf_train$ljung_box$p.value,
    diag_kf_full$ljung_box$p.value
  ),
  ARCH_p = c(
    diag_kf_train$arch$p.value,
    diag_kf_full$arch$p.value
  ),
  Break_years = c(
    paste(
      diag_kf_train$break_years,
      collapse = ", "
    ),
    paste(
      diag_kf_full$break_years,
      collapse = ", "
    )
  ),
  CUSUM_statistic = c(
    unname(diag_kf_train$cusum_test$statistic),
    unname(diag_kf_full$cusum_test$statistic)
  ),
  CUSUM_p = c(
    diag_kf_train$cusum_test$p.value,
    diag_kf_full$cusum_test$p.value
  ),
  CUSUM_out_of_bounds_years = c(
    paste(
      diag_kf_train$cusum_df$Year[
        diag_kf_train$cusum_df$Out_of_bounds
      ],
      collapse = ", "
    ),
    paste(
      diag_kf_full$cusum_df$Year[
        diag_kf_full$cusum_df$Out_of_bounds
      ],
      collapse = ", "
    )
  )
)

print(
  trend_diagnostics_summary,
  row.names = FALSE
)

# ACF / PACF
gridExtra::grid.arrange(
  forecast::ggAcf(ts_train_kf_manual2) +
    ggtitle("ACF - training Kalman trend") +
    theme_minimal(),
  forecast::ggPacf(ts_train_kf_manual2) +
    ggtitle("PACF - training Kalman trend") +
    theme_minimal(),
  ncol = 2
)

gridExtra::grid.arrange(
  forecast::ggAcf(ts_full_kf_manual2) +
    ggtitle("ACF - full-sample Kalman trend") +
    theme_minimal(),
  forecast::ggPacf(ts_full_kf_manual2) +
    ggtitle("PACF - full-sample Kalman trend") +
    theme_minimal(),
  ncol = 2
)

# Bai-Perron 95% confidence intervals
ci_breaks_train_kf <- confint(
  diag_kf_train$breakpoints,
  breaks = 2,
  level = 0.95
)

ci_breaks_full_kf <- confint(
  diag_kf_full$breakpoints,
  breaks = 3,
  level = 0.95
)

print(ci_breaks_train_kf)
print(breakdates(ci_breaks_train_kf))

print(ci_breaks_full_kf)
print(breakdates(ci_breaks_full_kf))

# CUSUM plots
plot(diag_kf_train$cusum)
plot(diag_kf_full$cusum)

# =============================================================================
# BASELINE RESIDUAL DIAGNOSTICS
# =============================================================================

residual_diagnostics_baseline <- function(
  x,
  label,
  ljung_lag
) {

  shapiro_result <- shapiro.test(
    as.numeric(x)
  )

  adf_result <- suppressWarnings(
    adf.test(x)
  )

  kpss_result <- suppressWarnings(
    kpss.test(x)
  )

  ljung_result <- Box.test(
    x,
    lag = ljung_lag,
    type = "Ljung-Box"
  )

  arch_result <- FinTS::ArchTest(
    as.numeric(x)
  )

  data.frame(
    Sample = label,
    Variance = var(x),
    Shapiro_p = shapiro_result$p.value,
    ADF_p = adf_result$p.value,
    KPSS_p = kpss_result$p.value,
    Ljung_Box_lag = ljung_lag,
    Ljung_Box_statistic = unname(ljung_result$statistic),
    Ljung_Box_p = ljung_result$p.value,
    ARCH_p = arch_result$p.value
  )
}

residual_diagnostics_summary <- rbind(
  residual_diagnostics_baseline(
    ts_train_resid_kf_manual2,
    "Training 1990-2018",
    6
  ),
  residual_diagnostics_baseline(
    ts_full_resid_kf_manual2,
    "Full 1990-2023",
    7
  )
)

print(
  residual_diagnostics_summary,
  row.names = FALSE
)

# Standard residual diagnostic plots
checkresiduals(
  ts_train_resid_kf_manual2
)

checkresiduals(
  ts_full_resid_kf_manual2
)

acf(
  ts_train_resid_kf_manual2,
  main = "ACF of Kalman residuals - training"
)

pacf(
  ts_train_resid_kf_manual2,
  main = "PACF of Kalman residuals - training"
)

acf(
  ts_full_resid_kf_manual2,
  main = "ACF of Kalman residuals - full sample"
)

pacf(
  ts_full_resid_kf_manual2,
  main = "PACF of Kalman residuals - full sample"
)

# =============================================================================
# CONSTRAINED-ML KALMAN CALIBRATION
# =============================================================================

# ---- Training sample ----------------------------------------------------------

model_kf_train_testC <- SSModel(
  ts_kf_train ~ SSMtrend(
    degree = 2,
    Q = list(NA, NA)
  ),
  H = 9999999
)

fit_kf_train_testC <- fitSSM(
  model_kf_train_testC,
  inits = c(0.1, 0.1),
  method = "Nelder-Mead"
)

kfs_train_testC <- KFS(
  fit_kf_train_testC$model,
  filtering = "state",
  smoothing = "state"
)

train_kf_testC <- kfs_train_testC$alphahat[, 1]

ts_train_kf_testC <- ts(
  train_kf_testC,
  start = 1990,
  frequency = 1
)

Q_train_testC <- kfs_train_testC$model$Q[, , 1]
H_train_testC_matrix <- kfs_train_testC$model$H[, , 1]

Q11_train_testC <- Q_train_testC[1, 1]
Q22_train_testC <- Q_train_testC[2, 2]
H_train_testC <- H_train_testC_matrix[1]

bp_train_testC <- breakpoints(
  ts_train_kf_testC ~ 1
)

cat(
  "\nConstrained-ML training calibration:\n",
  "Q11 =", Q11_train_testC, "\n",
  "Q22 =", Q22_train_testC, "\n",
  "H   =", H_train_testC, "\n",
  "Convergence code =", fit_kf_train_testC$optim.out$convergence, "\n"
)

print(
  summary(bp_train_testC)
)

# ---- Full sample --------------------------------------------------------------

model_kf_full_testC <- SSModel(
  ts_kf_full ~ SSMtrend(
    degree = 2,
    Q = list(NA, NA)
  ),
  H = 9999999
)

fit_kf_full_testC <- fitSSM(
  model_kf_full_testC,
  inits = c(0.1, 0.1),
  method = "Nelder-Mead"
)

kfs_full_testC <- KFS(
  fit_kf_full_testC$model,
  filtering = "state",
  smoothing = "state"
)

full_kf_testC <- kfs_full_testC$alphahat[, 1]

ts_full_kf_testC <- ts(
  full_kf_testC,
  start = 1990,
  frequency = 1
)

Q_full_testC <- kfs_full_testC$model$Q[, , 1]
H_full_testC_matrix <- kfs_full_testC$model$H[, , 1]

Q11_full_testC <- Q_full_testC[1, 1]
Q22_full_testC <- Q_full_testC[2, 2]
H_full_testC <- H_full_testC_matrix[1]

bp_full_testC <- breakpoints(
  ts_full_kf_testC ~ 1
)

cat(
  "\nConstrained-ML full-sample calibration:\n",
  "Q11 =", Q11_full_testC, "\n",
  "Q22 =", Q22_full_testC, "\n",
  "H   =", H_full_testC, "\n",
  "Convergence code =", fit_kf_full_testC$optim.out$convergence, "\n"
)

print(
  summary(bp_full_testC)
)

# =============================================================================
# BASELINE AND CONSTRAINED-ML COMPARISON
# =============================================================================

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

baseline_train <- window(
  ts_train_kf_manual2,
  start = 1990,
  end = 2018
)

baseline_full <- window(
  ts_full_kf_manual2,
  start = 1990,
  end = 2023
)

testC_train <- window(
  ts_train_kf_testC,
  start = 1990,
  end = 2018
)

testC_full <- window(
  ts_full_kf_testC,
  start = 1990,
  end = 2023
)

smoothing_diagnostics <- function(
  raw,
  trend,
  label
) {

  d_raw <- diff(raw)
  d_trend <- diff(trend)

  d2_raw <- diff(d_raw)
  d2_trend <- diff(d_trend)

  direction_match <- sign(d_raw) == sign(d_trend)

  data.frame(
    Series = label,
    Var_raw_level = var(raw),
    Var_trend_level = var(trend),
    Trend_to_raw_level_variance = var(trend) / var(raw),
    Var_raw_D1 = var(d_raw),
    Var_trend_D1 = var(d_trend),
    Trend_to_raw_D1_variance = var(d_trend) / var(d_raw),
    Var_raw_D2 = var(d2_raw),
    Var_trend_D2 = var(d2_trend),
    Trend_to_raw_D2_variance = var(d2_trend) / var(d2_raw),
    Correlation_levels = cor(raw, trend),
    Correlation_changes = cor(d_raw, d_trend),
    Direction_matches = sum(direction_match),
    Direction_total = length(direction_match),
    Direction_agreement_pct = 100 * mean(direction_match),
    MAE_raw_vs_trend = mean(abs(raw - trend)),
    RMSE_raw_vs_trend = sqrt(mean((raw - trend)^2))
  )
}

smoothing_training <- rbind(
  smoothing_diagnostics(
    raw_train,
    baseline_train,
    "Baseline KF - training"
  ),
  smoothing_diagnostics(
    raw_train,
    testC_train,
    "Test C KF - training"
  )
)

smoothing_full <- rbind(
  smoothing_diagnostics(
    raw_full,
    baseline_full,
    "Baseline KF - full"
  ),
  smoothing_diagnostics(
    raw_full,
    testC_full,
    "Test C KF - full"
  )
)

print(
  smoothing_training,
  row.names = FALSE
)

print(
  smoothing_full,
  row.names = FALSE
)

resid_baseline_train <- raw_train - baseline_train
resid_testC_train <- raw_train - testC_train

resid_baseline_full <- raw_full - baseline_full
resid_testC_full <- raw_full - testC_full

residual_diagnostics_comparison <- function(
  resid,
  label
) {

  adf_res <- suppressWarnings(
    adf.test(resid)
  )

  kpss_res <- suppressWarnings(
    kpss.test(resid)
  )

  lb_res <- Box.test(
    resid,
    lag = 7,
    type = "Ljung-Box"
  )

  arch_res <- FinTS::ArchTest(
    as.numeric(resid),
    lags = 5
  )

  shapiro_res <- shapiro.test(
    as.numeric(resid)
  )

  data.frame(
    Series = label,
    Mean = mean(resid),
    Variance = var(resid),
    SD = sd(resid),
    Max_abs_residual = max(abs(resid)),
    ADF_p = adf_res$p.value,
    KPSS_p = kpss_res$p.value,
    Ljung_Box_p = lb_res$p.value,
    ARCH_p = arch_res$p.value,
    Shapiro_p = shapiro_res$p.value
  )
}

residual_comparison <- rbind(
  residual_diagnostics_comparison(
    resid_baseline_train,
    "Baseline KF - training"
  ),
  residual_diagnostics_comparison(
    resid_testC_train,
    "Test C KF - training"
  ),
  residual_diagnostics_comparison(
    resid_baseline_full,
    "Baseline KF - full"
  ),
  residual_diagnostics_comparison(
    resid_testC_full,
    "Test C KF - full"
  )
)

print(
  residual_comparison,
  row.names = FALSE
)

get_bic_breaks_compact <- function(
  series,
  label
) {

  bp_full <- strucchange::breakpoints(
    series ~ 1
  )

  bp_index <- bp_full$breakpoints
  bp_index <- bp_index[!is.na(bp_index)]

  bp_years <- if (length(bp_index) == 0L) {
    numeric(0)
  } else {
    as.numeric(
      time(series)[bp_index]
    )
  }

  data.frame(
    Series = label,
    BIC_breaks = length(bp_index),
    Break_years = if (length(bp_years) == 0L) {
      ""
    } else {
      paste(
        bp_years,
        collapse = ", "
      )
    }
  )
}

bp_comparison <- rbind(
  get_bic_breaks_compact(
    baseline_train,
    "Baseline KF - training"
  ),
  get_bic_breaks_compact(
    testC_train,
    "Test C KF - training"
  ),
  get_bic_breaks_compact(
    baseline_full,
    "Baseline KF - full"
  ),
  get_bic_breaks_compact(
    testC_full,
    "Test C KF - full"
  )
)

print(
  bp_comparison,
  row.names = FALSE
)

kalman_sensitivity_summary <- data.frame(
  Calibration = c(
    "Baseline",
    "Test C",
    "Baseline",
    "Test C"
  ),
  Sample = c(
    "Training",
    "Training",
    "Full",
    "Full"
  ),
  Directional_agreement_pct = c(
    smoothing_training$Direction_agreement_pct[
      smoothing_training$Series == "Baseline KF - training"
    ],
    smoothing_training$Direction_agreement_pct[
      smoothing_training$Series == "Test C KF - training"
    ],
    smoothing_full$Direction_agreement_pct[
      smoothing_full$Series == "Baseline KF - full"
    ],
    smoothing_full$Direction_agreement_pct[
      smoothing_full$Series == "Test C KF - full"
    ]
  ),
  Retained_D2_variance_pct = 100 * c(
    smoothing_training$Trend_to_raw_D2_variance[
      smoothing_training$Series == "Baseline KF - training"
    ],
    smoothing_training$Trend_to_raw_D2_variance[
      smoothing_training$Series == "Test C KF - training"
    ],
    smoothing_full$Trend_to_raw_D2_variance[
      smoothing_full$Series == "Baseline KF - full"
    ],
    smoothing_full$Trend_to_raw_D2_variance[
      smoothing_full$Series == "Test C KF - full"
    ]
  ),
  Residual_variance_million = c(
    residual_comparison$Variance[
      residual_comparison$Series == "Baseline KF - training"
    ],
    residual_comparison$Variance[
      residual_comparison$Series == "Test C KF - training"
    ],
    residual_comparison$Variance[
      residual_comparison$Series == "Baseline KF - full"
    ],
    residual_comparison$Variance[
      residual_comparison$Series == "Test C KF - full"
    ]
  ) / 1e6,
  Residual_Ljung_Box_p = c(
    residual_comparison$Ljung_Box_p[
      residual_comparison$Series == "Baseline KF - training"
    ],
    residual_comparison$Ljung_Box_p[
      residual_comparison$Series == "Test C KF - training"
    ],
    residual_comparison$Ljung_Box_p[
      residual_comparison$Series == "Baseline KF - full"
    ],
    residual_comparison$Ljung_Box_p[
      residual_comparison$Series == "Test C KF - full"
    ]
  ),
  BIC_break_years = c(
    bp_comparison$Break_years[
      bp_comparison$Series == "Baseline KF - training"
    ],
    bp_comparison$Break_years[
      bp_comparison$Series == "Test C KF - training"
    ],
    bp_comparison$Break_years[
      bp_comparison$Series == "Baseline KF - full"
    ],
    bp_comparison$Break_years[
      bp_comparison$Series == "Test C KF - full"
    ]
  )
)

print(
  kalman_sensitivity_summary,
  row.names = FALSE
)

sessionInfo()

# End of script
