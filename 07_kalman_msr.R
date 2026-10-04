# =============================================================================
# 07_kalman_msr.R
#
# MSR models on the baseline Kalman-smoothed trend.
#
# Input:
#   data/Total_Tonkm.xlsx
# =============================================================================

# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "KFAS",
  "MSwM",
  "nlme"
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
library(MSwM)
library(nlme)

# =============================================================================
# DATA 
# =============================================================================


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

ts_kf_full <- ts(
  Total_Tonkm$Total,
  start = Total_Tonkm$Year[1],
  frequency = 1
)

# =============================================================================
# BASELINE KALMAN-SMOOTHED TREND: TRAINING 1990-2018
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

train_resid_kf_manual2 <- ts_kf_train - train_kf_manual2

ts_train_kf_manual2 <- ts(
  train_kf_manual2,
  start = Total_Tonkm$Year[1],
  frequency = 1
)

ts_train_resid_kf_manual2 <- ts(
  train_resid_kf_manual2,
  start = Total_Tonkm$Year[1],
  frequency = 1
)

# =============================================================================
# BASELINE KALMAN-SMOOTHED TREND: FULL 1990-2023
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

full_resid_kf_manual2 <- ts_kf_full - full_kf_manual2

ts_full_kf_manual2 <- ts(
  full_kf_manual2,
  start = Total_Tonkm$Year[1],
  frequency = 1
)

ts_full_resid_kf_manual2 <- ts(
  full_resid_kf_manual2,
  start = Total_Tonkm$Year[1],
  frequency = 1
)

# =============================================================================
# MSR(2,0): TRAINING 1990-2018
# =============================================================================

lm_train <- lm(
  train_kf_manual2 ~ 1
)

msm_train_k2p0 <- msmFit(
  lm_train,
  k = 2,
  p = 0,
  sw = c(TRUE, TRUE),
  control = list(
    parallel = FALSE,
    maxiter = 5000,
    trace = TRUE
  )
)

summary(msm_train_k2p0)

ll_train <- msm_train_k2p0@Fit@logLikel
n_train <- nrow(msm_train_k2p0@model$model)

swi_train <- msm_train_k2p0@switch[
  -length(msm_train_k2p0@switch)
]

np_train <- msm_train_k2p0["k"] * sum(swi_train) +
  sum(!swi_train)

AIC_mswm_train <- 2 * ll_train + 2 * np_train
BIC_mswm_train <- 2 * ll_train + 2 * np_train * log(n_train)
logLik_train_reported <- -ll_train

# =============================================================================
# MSR(2,0): FULL 1990-2023
# =============================================================================

lm_full <- lm(
  ts_full_kf_manual2 ~ 1
)

msm_full_k2p0 <- msmFit(
  lm_full,
  k = 2,
  p = 0,
  sw = c(TRUE, TRUE),
  control = list(
    parallel = FALSE,
    maxiter = 5000,
    trace = TRUE
  )
)

summary(msm_full_k2p0)

ll_full <- msm_full_k2p0@Fit@logLikel
n_full <- nrow(msm_full_k2p0@model$model)

swi_full <- msm_full_k2p0@switch[
  -length(msm_full_k2p0@switch)
]

np_full <- msm_full_k2p0["k"] * sum(swi_full) +
  sum(!swi_full)

AIC_mswm_full <- 2 * ll_full + 2 * np_full
BIC_mswm_full <- 2 * ll_full + 2 * np_full * log(n_full)
logLik_full_reported <- -ll_full

# =============================================================================
# DETERMINISTIC REGIME PROJECTION: 2019-2023
# =============================================================================

P <- msm_train_k2p0@transMat

pi_t <- msm_train_k2p0@Fit@smoProb[
  nrow(msm_train_k2p0@Fit@smoProb),
]

PTregime_forecastk2 <- matrix(
  NA_real_,
  5,
  2
)

for (h in 1:5) {

  pi_t <- as.numeric(
    pi_t %*% t(P)
  )

  pi_t <- pi_t / sum(pi_t)

  PTregime_forecastk2[h, ] <- pi_t
}

years_forecastk2p0 <- 2019:2023

df_regime_forecastk2 <- data.frame(
  Year = years_forecastk2p0,
  PTRegime1_Probk2 = PTregime_forecastk2[, 1],
  PTRegime2_Probk2 = PTregime_forecastk2[, 2],
  PTPredicted_Regimek2 = apply(
    PTregime_forecastk2,
    1,
    which.max
  )
)

df_regime_forecastk2$PTProb_Predicted_Regime <- mapply(
  function(reg, p1, p2) {
    c(p1, p2)[reg]
  },
  df_regime_forecastk2$PTPredicted_Regimek2,
  df_regime_forecastk2$PTRegime1_Probk2,
  df_regime_forecastk2$PTRegime2_Probk2
)

# =============================================================================
# EXACT MSR REGIME-PATH PROBABILITIES: KALMAN-SMOOTHED SERIES
# =============================================================================

P_kf <- msm_train_k2p0@transMat

pi_2018_kf <- msm_train_k2p0@Fit@smoProb[
  nrow(msm_train_k2p0@Fit@smoProb),
]

future_paths_kf <- expand.grid(
  s2019 = 1:2,
  s2020 = 1:2,
  s2021 = 1:2,
  s2022 = 1:2,
  s2023 = 1:2
)

exact_future_prob_kf <- apply(
  future_paths_kf,
  1,
  function(path) {

    prob_path <- 0

    for (s2018 in 1:2) {

      prob_conditional <-
        pi_2018_kf[s2018] *
        P_kf[path[1], s2018]

      for (tt in 2:length(path)) {

        prob_conditional <-
          prob_conditional *
          P_kf[
            path[tt],
            path[tt - 1]
          ]
      }

      prob_path <- prob_path + prob_conditional
    }

    prob_path
  }
)

exact_paths_kf <- data.frame(
  Path = apply(
    future_paths_kf,
    1,
    function(x) {
      paste0(
        "R",
        x,
        collapse = "-"
      )
    }
  ),
  Probability = exact_future_prob_kf,
  Percent = 100 * exact_future_prob_kf
)

exact_paths_kf <- exact_paths_kf[
  order(
    exact_paths_kf$Probability,
    decreasing = TRUE
  ),
]

rownames(exact_paths_kf) <- NULL

class_probability_kf <- c(
  no_switch = 0,
  single_switch = 0,
  repeated_switching = 0
)

complete_paths_kf <- data.frame()

for (i in seq_len(nrow(future_paths_kf))) {

  future_path <- as.integer(
    future_paths_kf[i, ]
  )

  for (s2018 in 1:2) {

    complete_path <- c(
      s2018,
      future_path
    )

    prob_complete <- pi_2018_kf[s2018]

    for (tt in 2:length(complete_path)) {

      prob_complete <-
        prob_complete *
        P_kf[
          complete_path[tt],
          complete_path[tt - 1]
        ]
    }

    n_switches <- sum(
      diff(complete_path) != 0
    )

    if (n_switches == 0) {

      class_probability_kf["no_switch"] <-
        class_probability_kf["no_switch"] +
        prob_complete

    } else if (n_switches == 1) {

      class_probability_kf["single_switch"] <-
        class_probability_kf["single_switch"] +
        prob_complete

    } else {

      class_probability_kf["repeated_switching"] <-
        class_probability_kf["repeated_switching"] +
        prob_complete
    }

    complete_paths_kf <- rbind(
      complete_paths_kf,
      data.frame(
        Path = paste0(
          "R",
          complete_path,
          collapse = "-"
        ),
        Initial_2018 = s2018,
        Switches = n_switches,
        Probability = prob_complete,
        Percent = 100 * prob_complete
      )
    )
  }
}

exact_classes_kf <- data.frame(
  Path_class = c(
    "No switch",
    "Single switch",
    "Repeated switching"
  ),
  Number_of_transitions = c(
    "0",
    "1",
    ">=2"
  ),
  Probability = as.numeric(class_probability_kf),
  Percent = 100 * as.numeric(class_probability_kf)
)

complete_paths_kf <- complete_paths_kf[
  order(
    complete_paths_kf$Probability,
    decreasing = TRUE
  ),
]

rownames(complete_paths_kf) <- NULL

# =============================================================================
# TRANSITION PROBABILITIES AND 95% CONFIDENCE INTERVALS
# =============================================================================

transition_ci_mswm <- function(obj, label) {

  par0 <- c(
    log(obj@std[1]),
    log(obj@std[2]),
    qlogis(obj@transMat[1, 1]),
    qlogis(obj@transMat[1, 2]),
    obj@Coef[1, 1],
    obj@Coef[2, 1]
  )

  hess <- nlme::fdHess(
    pars = par0,
    fun = MSwM:::fopt.lm,
    object = obj
  )

  H <- hess$Hessian

  eig_H <- eigen(
    H,
    symmetric = TRUE
  )$values

  H_trans <- H[
    3:4,
    3:4
  ]

  V <- solve(H)

  eta_hat <- par0[3:4]

  se_eta <- sqrt(
    diag(V)
  )[3:4]

  z <- qnorm(0.975)

  lower_eta <- eta_hat - z * se_eta
  upper_eta <- eta_hat + z * se_eta

  p_hat <- plogis(eta_hat)
  p_lower <- plogis(lower_eta)
  p_upper <- plogis(upper_eta)

  p11 <- p_hat[1]
  p12 <- p_hat[2]

  p21 <- 1 - p11
  p22 <- 1 - p12

  p11_L <- p_lower[1]
  p11_U <- p_upper[1]

  p21_L <- 1 - p11_U
  p21_U <- 1 - p11_L

  p12_L <- p_lower[2]
  p12_U <- p_upper[2]

  p22_L <- 1 - p12_U
  p22_U <- 1 - p12_L

  transition_CI <- data.frame(
    Transition = c(
      "Regime 1 -> Regime 1",
      "Regime 1 -> Regime 2",
      "Regime 2 -> Regime 1",
      "Regime 2 -> Regime 2"
    ),
    Estimate = c(
      p11,
      p21,
      p12,
      p22
    ),
    CI_2.5 = c(
      p11_L,
      p21_L,
      p12_L,
      p22_L
    ),
    CI_97.5 = c(
      p11_U,
      p21_U,
      p12_U,
      p22_U
    )
  )

  cat("\n\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")

  cat("\nRegime means:\n")
  print(obj@Coef)

  cat("\nEstimated transition matrix:\n")
  print(obj@transMat)

  cat("\nHessian eigenvalues:\n")
  print(eig_H)

  cat("\nCondition number - transition block:\n")
  print(kappa(H_trans))

  cat("\nTransition parameters - logit scale:\n")
  print(
    data.frame(
      Parameter = c(
        "logit P(1|1)",
        "logit P(1|2)"
      ),
      Estimate = eta_hat,
      SE = se_eta
    )
  )

  cat("\nTransition probabilities with 95% CI:\n")
  print(
    transition_CI,
    digits = 6
  )

  invisible(transition_CI)
}

ci_kf_train <- transition_ci_mswm(
  msm_train_k2p0,
  "KALMAN TRAINING 1990-2018"
)

ci_kf_full <- transition_ci_mswm(
  msm_full_k2p0,
  "KALMAN FULL 1990-2023"
)

# =============================================================================
# SMOOTHED REGIME PROBABILITIES
# =============================================================================

extract_smoothed_probs <- function(
  obj,
  years,
  label
) {

  probs <- obj@Fit@smoProb[-1, ]

  means <- as.numeric(
    obj@Coef[, 1]
  )

  low_regime <- which.min(means)
  high_regime <- which.max(means)

  out <- data.frame(
    Year = years,
    Prob_Low = probs[, low_regime],
    Prob_High = probs[, high_regime]
  )

  out$Dominant_Regime <- ifelse(
    out$Prob_Low >= out$Prob_High,
    "Low",
    "High"
  )

  out$Dominant_Prob <- pmax(
    out$Prob_Low,
    out$Prob_High
  )

  cat("\n\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")

  cat("\nRegime means:\n")
  print(means)

  cat(
    "\nLow regime =",
    low_regime,
    "| mean =",
    means[low_regime],
    "\n"
  )

  cat(
    "High regime =",
    high_regime,
    "| mean =",
    means[high_regime],
    "\n"
  )

  cat("\n===== ALL SMOOTHED REGIME PROBABILITIES =====\n")
  print(
    out,
    digits = 6
  )

  invisible(out)
}

kf_prob_train <- extract_smoothed_probs(
  obj = msm_train_k2p0,
  years = 1990:2018,
  label = "KALMAN MSR - TRAINING 1990-2018"
)

kf_prob_full <- extract_smoothed_probs(
  obj = msm_full_k2p0,
  years = 1990:2023,
  label = "KALMAN MSR - FULL 1990-2023"
)

# =============================================================================
# ALTERNATIVE MSR SPECIFICATIONS
# =============================================================================

msr_control_kf <- list(
  parallel = FALSE,
  maxiter = 5000,
  tol = 1e-8,
  trace = TRUE
)

fit_msr_kf_safe <- function(
  y,
  k,
  p,
  sw,
  label
) {

  cat("\n\n")
  cat("====================================================================\n")
  cat(label, "\n")
  cat("====================================================================\n")

  base_lm <- lm(
    as.numeric(y) ~ 1
  )

  fit <- tryCatch(
    msmFit(
      base_lm,
      k = k,
      p = p,
      sw = sw,
      control = msr_control_kf
    ),
    error = function(e) {

      cat("\nFIT FAILED:\n")
      cat(
        conditionMessage(e),
        "\n"
      )

      return(NULL)
    }
  )

  if (is.null(fit)) {

    return(
      list(
        fit = NULL,
        fit_ok = FALSE,
        summary_ok = FALSE,
        std_ok = FALSE,
        trans_ok = FALSE,
        error = "Fit failed"
      )
    )
  }

  cat("\nMODEL OBJECT RETURNED.\n")

  summary_result <- tryCatch(
    summary(fit),
    error = function(e) e
  )

  if (inherits(summary_result, "error")) {

    cat("\nSUMMARY FAILED:\n")
    cat(
      conditionMessage(summary_result),
      "\n"
    )

    summary_ok <- FALSE
    summary_error <- conditionMessage(
      summary_result
    )

  } else {

    print(summary_result)

    summary_ok <- TRUE
    summary_error <- NA_character_
  }

  std_ok <- tryCatch(
    {

      std_values <- suppressWarnings(
        as.numeric(fit@std)
      )

      cat("\nResidual standard deviations:\n")
      print(fit@std)

      length(std_values) > 0 &&
        all(is.finite(std_values))
    },
    error = function(e) {

      cat("\nRESIDUAL-SD EXTRACTION FAILED:\n")
      cat(
        conditionMessage(e),
        "\n"
      )

      FALSE
    }
  )

  trans_ok <- tryCatch(
    {

      cat("\nTransition matrix:\n")
      print(fit@transMat)

      all(
        is.finite(
          as.numeric(fit@transMat)
        )
      )
    },
    error = function(e) {

      cat("\nTRANSITION-MATRIX EXTRACTION FAILED:\n")
      cat(
        conditionMessage(e),
        "\n"
      )

      FALSE
    }
  )

  cat("\nSTATUS\n")
  cat("Fit returned:          TRUE\n")
  cat("Summary available:    ", summary_ok, "\n")
  cat("Finite residual SDs:   ", std_ok, "\n")
  cat("Transition matrix OK: ", trans_ok, "\n")

  return(
    list(
      fit = fit,
      fit_ok = TRUE,
      summary_ok = summary_ok,
      std_ok = std_ok,
      trans_ok = trans_ok,
      error = summary_error
    )
  )
}

kfcheck_train_k2p0 <- fit_msr_kf_safe(
  y = ts_train_kf_manual2,
  k = 2,
  p = 0,
  sw = c(TRUE, TRUE),
  label = "KF TRAINING 1990-2018 - MSR(2,0)"
)

kfcheck_train_k2p1 <- fit_msr_kf_safe(
  y = ts_train_kf_manual2,
  k = 2,
  p = 1,
  sw = c(TRUE, TRUE, TRUE),
  label = "KF TRAINING 1990-2018 - MSR(2,1)"
)

kfcheck_train_k3p0 <- fit_msr_kf_safe(
  y = ts_train_kf_manual2,
  k = 3,
  p = 0,
  sw = c(TRUE, TRUE),
  label = "KF TRAINING 1990-2018 - MSR(3,0)"
)

kfcheck_train_k3p1 <- fit_msr_kf_safe(
  y = ts_train_kf_manual2,
  k = 3,
  p = 1,
  sw = c(TRUE, TRUE, TRUE),
  label = "KF TRAINING 1990-2018 - MSR(3,1)"
)

kfcheck_full_k2p0 <- fit_msr_kf_safe(
  y = ts_full_kf_manual2,
  k = 2,
  p = 0,
  sw = c(TRUE, TRUE),
  label = "KF FULL 1990-2023 - MSR(2,0)"
)

kfcheck_full_k2p1 <- fit_msr_kf_safe(
  y = ts_full_kf_manual2,
  k = 2,
  p = 1,
  sw = c(TRUE, TRUE, TRUE),
  label = "KF FULL 1990-2023 - MSR(2,1)"
)

kfcheck_full_k3p0 <- fit_msr_kf_safe(
  y = ts_full_kf_manual2,
  k = 3,
  p = 0,
  sw = c(TRUE, TRUE),
  label = "KF FULL 1990-2023 - MSR(3,0)"
)

kfcheck_full_k3p1 <- fit_msr_kf_safe(
  y = ts_full_kf_manual2,
  k = 3,
  p = 1,
  sw = c(TRUE, TRUE, TRUE),
  label = "KF FULL 1990-2023 - MSR(3,1)"
)

kf_fits <- list(
  kfcheck_train_k2p0,
  kfcheck_train_k2p1,
  kfcheck_train_k3p0,
  kfcheck_train_k3p1,
  kfcheck_full_k2p0,
  kfcheck_full_k2p1,
  kfcheck_full_k3p0,
  kfcheck_full_k3p1
)

kf_msr_status <- data.frame(
  Sample = c(
    rep(
      "KF training 1990-2018",
      4
    ),
    rep(
      "KF full 1990-2023",
      4
    )
  ),
  Specification = rep(
    c(
      "MSR(2,0)",
      "MSR(2,1)",
      "MSR(3,0)",
      "MSR(3,1)"
    ),
    2
  ),
  Fit_returned = vapply(
    kf_fits,
    function(x) {
      isTRUE(x$fit_ok)
    },
    logical(1)
  ),
  Summary_available = vapply(
    kf_fits,
    function(x) {
      isTRUE(x$summary_ok)
    },
    logical(1)
  ),
  Finite_residual_SDs = vapply(
    kf_fits,
    function(x) {
      isTRUE(x$std_ok)
    },
    logical(1)
  ),
  Transition_matrix = vapply(
    kf_fits,
    function(x) {
      isTRUE(x$trans_ok)
    },
    logical(1)
  )
)

# =============================================================================
# OUTPUT TO EXTRACT
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("OUTPUT TO EXTRACT\n")
cat("####################################################################\n")

cat("\n--- 1. MSR(2,0) TRAINING: MODEL SUMMARY ---\n")
print(
  summary(msm_train_k2p0)
)

cat("\nTraining log likelihood:\n")
print(logLik_train_reported)

cat("\nTraining AIC:\n")
print(AIC_mswm_train)

cat("\nTraining BIC:\n")
print(BIC_mswm_train)

cat("\nTraining regime means:\n")
print(
  msm_train_k2p0@Coef,
  digits = 12
)

cat("\nTraining regime residual SDs:\n")
print(
  msm_train_k2p0@std,
  digits = 12
)

cat("\nTraining transition matrix:\n")
print(
  msm_train_k2p0@transMat,
  digits = 12
)

cat("\n--- 2. DETERMINISTIC REGIME PROJECTION 2019-2023 ---\n")
print(
  df_regime_forecastk2,
  row.names = FALSE,
  digits = 12
)

cat("\n--- 3. EXACT REGIME-PATH PROBABILITIES ---\n")
print(
  exact_classes_kf,
  row.names = FALSE,
  digits = 12
)

cat("\nProbability mass by exact number of switches:\n")
print(
  aggregate(
    Probability ~ Switches,
    data = complete_paths_kf,
    FUN = sum
  ),
  row.names = FALSE,
  digits = 12
)

cat("\nFinal-paper path classes:\n")

cat(
  "No switch (0): ",
  sum(
    complete_paths_kf$Probability[
      complete_paths_kf$Switches == 0
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Single switch (1): ",
  sum(
    complete_paths_kf$Probability[
      complete_paths_kf$Switches == 1
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Switch and return (2): ",
  sum(
    complete_paths_kf$Probability[
      complete_paths_kf$Switches == 2
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "More than two switches (>2): ",
  sum(
    complete_paths_kf$Probability[
      complete_paths_kf$Switches > 2
    ]
  ),
  "\n",
  sep = ""
)

cat("\nCheck: sum of all 32 future-path probabilities:\n")
print(
  sum(exact_paths_kf$Probability),
  digits = 12
)

cat("\n--- 4. MSR(2,0) FULL: MODEL SUMMARY ---\n")
print(
  summary(msm_full_k2p0)
)

cat("\nFull-sample log likelihood:\n")
print(logLik_full_reported)

cat("\nFull-sample AIC:\n")
print(AIC_mswm_full)

cat("\nFull-sample BIC:\n")
print(BIC_mswm_full)

cat("\nFull-sample regime means:\n")
print(
  msm_full_k2p0@Coef,
  digits = 12
)

cat("\nFull-sample regime residual SDs:\n")
print(
  msm_full_k2p0@std,
  digits = 12
)

cat("\nFull-sample transition matrix:\n")
print(
  msm_full_k2p0@transMat,
  digits = 12
)

cat("\n--- 5. TRANSITION-PROBABILITY 95% CIs ---\n")

cat("\nTraining:\n")
print(
  ci_kf_train,
  row.names = FALSE,
  digits = 12
)

cat("\nFull sample:\n")
print(
  ci_kf_full,
  row.names = FALSE,
  digits = 12
)

cat("\n--- 6. SMOOTHED PROBABILITIES: TABLE YEARS ---\n")

table26_years <- c(
  1994,
  1995,
  1996,
  2010,
  2011,
  2012
)

cat("\nTraining:\n")
print(
  kf_prob_train[
    kf_prob_train$Year %in% table26_years,
  ],
  row.names = FALSE,
  digits = 12
)

cat("\nFull sample:\n")
print(
  kf_prob_full[
    kf_prob_full$Year %in% table26_years,
  ],
  row.names = FALSE,
  digits = 12
)

cat("\n--- 7. ALTERNATIVE MSR SPECIFICATION STATUS ---\n")
print(
  kf_msr_status,
  row.names = FALSE
)

# End of script
