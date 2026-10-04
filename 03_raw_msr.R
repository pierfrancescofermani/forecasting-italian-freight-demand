# =============================================================================
# 03_raw_msr.R
#
# Markov-switching regression on the raw freight-activity series.
#
# Input:
# Total_Tonkm.xlsx
# =============================================================================


# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "MSwM",
  "ggplot2",
  "dplyr",
  "tidyr",
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

library(MSwM)
library(ggplot2)
library(dplyr)
library(tidyr)
library(nlme)

# ---- Data --------------------------------------------------------------------

df_raw_train <- Total_Tonkm[
  Total_Tonkm$Year <= 2018,
]

df_raw_full <- Total_Tonkm[
  Total_Tonkm$Year <= 2023,
]

# =============================================================================
# MODEL ESTIMATION
# =============================================================================

msr_control <- list(
  parallel = FALSE,
  maxiter = 5000,
  tol = 1e-8,
  trace = TRUE
)

fit_msr_safe <- function(data, k, p, sw, label) {

  base_lm <- lm(
    Total ~ 1,
    data = data
  )

  fit <- tryCatch(
    msmFit(
      base_lm,
      k = k,
      p = p,
      sw = sw,
      control = msr_control
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(
      list(
        label = label,
        fit = NULL,
        fit_ok = FALSE,
        summary_ok = FALSE,
        std_ok = FALSE,
        trans_ok = FALSE
      )
    )
  }

  summary_result <- tryCatch(
    summary(fit),
    error = function(e) e
  )

  summary_ok <- !inherits(
    summary_result,
    "error"
  )

  std_ok <- tryCatch(
    {
      std_values <- suppressWarnings(
        as.numeric(fit@std)
      )
      length(std_values) > 0L &&
        all(is.finite(std_values))
    },
    error = function(e) FALSE
  )

  trans_ok <- tryCatch(
    {
      P <- fit@transMat
      all(is.finite(P))
    },
    error = function(e) FALSE
  )

  list(
    label = label,
    fit = fit,
    fit_ok = TRUE,
    summary_ok = summary_ok,
    std_ok = std_ok,
    trans_ok = trans_ok
  )
}

# ---- Training 1990-2018 -------------------------------------------------------

raw_train_k2p0 <- fit_msr_safe(
  data = df_raw_train,
  k = 2,
  p = 0,
  sw = c(TRUE, TRUE),
  label = "Raw training 1990-2018 - MSR(2,0)"
)

raw_train_k2p1 <- fit_msr_safe(
  data = df_raw_train,
  k = 2,
  p = 1,
  sw = c(TRUE, TRUE, TRUE),
  label = "Raw training 1990-2018 - MSR(2,1)"
)

raw_train_k3p0 <- fit_msr_safe(
  data = df_raw_train,
  k = 3,
  p = 0,
  sw = c(TRUE, TRUE),
  label = "Raw training 1990-2018 - MSR(3,0)"
)

raw_train_k3p1 <- fit_msr_safe(
  data = df_raw_train,
  k = 3,
  p = 1,
  sw = c(TRUE, TRUE, TRUE),
  label = "Raw training 1990-2018 - MSR(3,1)"
)

# ---- Full 1990-2023 -----------------------------------------------------------

raw_full_k2p0 <- fit_msr_safe(
  data = df_raw_full,
  k = 2,
  p = 0,
  sw = c(TRUE, TRUE),
  label = "Raw full 1990-2023 - MSR(2,0)"
)

raw_full_k2p1 <- fit_msr_safe(
  data = df_raw_full,
  k = 2,
  p = 1,
  sw = c(TRUE, TRUE, TRUE),
  label = "Raw full 1990-2023 - MSR(2,1)"
)

raw_full_k3p0 <- fit_msr_safe(
  data = df_raw_full,
  k = 3,
  p = 0,
  sw = c(TRUE, TRUE),
  label = "Raw full 1990-2023 - MSR(3,0)"
)

raw_full_k3p1 <- fit_msr_safe(
  data = df_raw_full,
  k = 3,
  p = 1,
  sw = c(TRUE, TRUE, TRUE),
  label = "Raw full 1990-2023 - MSR(3,1)"
)

fits_msr <- list(
  raw_train_k2p0,
  raw_train_k2p1,
  raw_train_k3p0,
  raw_train_k3p1,
  raw_full_k2p0,
  raw_full_k2p1,
  raw_full_k3p0,
  raw_full_k3p1
)

msr_status <- data.frame(
  Sample = c(
    rep("Raw training 1990-2018", 4),
    rep("Raw full 1990-2023", 4)
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
    fits_msr,
    function(x) isTRUE(x$fit_ok),
    logical(1)
  ),
  Summary_available = vapply(
    fits_msr,
    function(x) isTRUE(x$summary_ok),
    logical(1)
  ),
  Finite_standard_errors = vapply(
    fits_msr,
    function(x) isTRUE(x$std_ok),
    logical(1)
  ),
  Transition_matrix = vapply(
    fits_msr,
    function(x) isTRUE(x$trans_ok),
    logical(1)
  )
)

print(msr_status, row.names = FALSE)

stopifnot(
  isTRUE(raw_train_k2p0$fit_ok),
  isTRUE(raw_full_k2p0$fit_ok)
)

msr_k2_train <- raw_train_k2p0$fit
msr_k2_full <- raw_full_k2p0$fit

# =============================================================================
# BASELINE MSR(2,0): TRAINING AND FULL SAMPLE
# =============================================================================

cat("\nMSR(2,0) - training sample\n")
print(summary(msr_k2_train))

cat("\nTransition matrix - training sample\n")
print(msr_k2_train@transMat)

cat("\nRegime-specific coefficients - training sample\n")
print(msr_k2_train@Coef)

cat("\nRegime-specific residual standard deviations - training sample\n")
print(msr_k2_train@std)

cat("\nMSR(2,0) - full sample\n")
print(summary(msr_k2_full))

cat("\nTransition matrix - full sample\n")
print(msr_k2_full@transMat)

cat("\nRegime-specific coefficients - full sample\n")
print(msr_k2_full@Coef)

cat("\nRegime-specific residual standard deviations - full sample\n")
print(msr_k2_full@std)

stopifnot(
  isTRUE(
    all.equal(
      colSums(msr_k2_train@transMat),
      rep(1, 2),
      tolerance = 1e-8
    )
  ),
  isTRUE(
    all.equal(
      colSums(msr_k2_full@transMat),
      rep(1, 2),
      tolerance = 1e-8
    )
  )
)

# =============================================================================
# SMOOTHED REGIME PROBABILITIES
# =============================================================================

extract_smoothed_probs <- function(obj, years, label) {

  probs <- obj@Fit@smoProb[-1, , drop = FALSE]
  means <- as.numeric(obj@Coef[, 1])

  stopifnot(
    nrow(probs) == length(years)
  )

  low_regime <- which.min(means)
  high_regime <- which.max(means)

  out <- data.frame(
    Series = label,
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

  out
}

raw_prob_train <- extract_smoothed_probs(
  obj = msr_k2_train,
  years = 1990:2018,
  label = "Raw training 1990-2018"
)

raw_prob_full <- extract_smoothed_probs(
  obj = msr_k2_full,
  years = 1990:2023,
  label = "Raw full 1990-2023"
)

print(raw_prob_train, row.names = FALSE)
print(raw_prob_full, row.names = FALSE)

# =============================================================================
# REGIME CLASSIFICATION PLOTS
# =============================================================================

regime_class_train <- apply(
  msr_k2_train@Fit@smoProb[-1, , drop = FALSE],
  1,
  which.max
)

df_plot_train <- data.frame(
  Year = df_raw_train$Year,
  TonKm = df_raw_train$Total,
  Regime = factor(regime_class_train)
)

plot_msr_train <- ggplot(
  df_plot_train,
  aes(x = Year, y = TonKm)
) +
  geom_line(
    color = "grey25",
    linewidth = 0.8
  ) +
  geom_point(
    aes(
      color = Regime,
      shape = Regime
    ),
    size = 2.2,
    alpha = 0.9
  ) +
  scale_color_manual(
    name = "Estimated regime",
    values = c(
      "#1F78B4",
      "#D95F02"
    ),
    labels = c(
      "Regime 1",
      "Regime 2"
    )
  ) +
  scale_shape_manual(
    name = "Estimated regime",
    values = c(16, 17),
    labels = c(
      "Regime 1",
      "Regime 2"
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

print(plot_msr_train)

regime_class_full <- apply(
  msr_k2_full@Fit@smoProb[-1, , drop = FALSE],
  1,
  which.max
)

df_regime_plot <- data.frame(
  Year = df_raw_full$Year,
  TonKm = df_raw_full$Total,
  Regime = factor(regime_class_full)
)

plot_msr_full <- ggplot(
  df_regime_plot,
  aes(x = Year, y = TonKm)
) +
  geom_line(
    color = "grey25",
    linewidth = 0.8
  ) +
  geom_point(
    aes(
      color = Regime,
      shape = Regime
    ),
    size = 2.2,
    alpha = 0.9
  ) +
  scale_color_manual(
    name = "Estimated regime",
    values = c(
      "1" = "#1F78B4",
      "2" = "#D95F02"
    )
  ) +
  scale_shape_manual(
    name = "Estimated regime",
    values = c(
      "1" = 16,
      "2" = 17
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

print(plot_msr_full)

# =============================================================================
# DETERMINISTIC REGIME PROJECTION: 2019-2023
# =============================================================================

P_raw <- msr_k2_train@transMat

pi_2018_raw <- msr_k2_train@Fit@smoProb[
  nrow(msr_k2_train@Fit@smoProb),
]

regime_forecast_raw <- matrix(
  NA_real_,
  nrow = 5,
  ncol = 2
)

regime_forecast_raw[1, ] <-
  pi_2018_raw %*% t(P_raw)

for (h in 2:5) {
  regime_forecast_raw[h, ] <-
    regime_forecast_raw[h - 1, ] %*% t(P_raw)
}

df_regime_forecast_raw <- data.frame(
  Year = 2019:2023,
  Regime1_Prob = regime_forecast_raw[, 1],
  Regime2_Prob = regime_forecast_raw[, 2],
  Predicted_Regime = apply(
    regime_forecast_raw,
    1,
    which.max
  )
)

df_regime_forecast_raw$Observed <-
  Total_Tonkm$Total[
    Total_Tonkm$Year %in% 2019:2023
  ]

df_regime_forecast_raw$Prob_Predicted_Regime <-
  mapply(
    function(reg, p1, p2) {
      if (reg == 1) p1 else p2
    },
    df_regime_forecast_raw$Predicted_Regime,
    df_regime_forecast_raw$Regime1_Prob,
    df_regime_forecast_raw$Regime2_Prob
  )

print(
  df_regime_forecast_raw,
  row.names = FALSE
)

df_train_projection <- data.frame(
  Year = df_plot_train$Year,
  TonKm = df_plot_train$TonKm,
  Regime = factor(df_plot_train$Regime),
  Prob = 1
)

df_test_projection <- data.frame(
  Year = df_regime_forecast_raw$Year,
  TonKm = df_regime_forecast_raw$Observed,
  Regime = factor(
    df_regime_forecast_raw$Predicted_Regime
  ),
  Prob = df_regime_forecast_raw$Prob_Predicted_Regime
)

df_projection_plot <- rbind(
  df_train_projection,
  df_test_projection
)

plot_msr_projection <- ggplot(
  df_projection_plot,
  aes(x = Year, y = TonKm)
) +
  geom_line(
    color = "grey25",
    linewidth = 0.8
  ) +
  geom_point(
    aes(
      color = Regime,
      shape = Regime,
      alpha = Prob
    ),
    size = 2.2
  ) +
  scale_color_manual(
    name = "Regime",
    values = c(
      "#1F78B4",
      "#D95F02"
    ),
    labels = c(
      "Regime 1",
      "Regime 2"
    )
  ) +
  scale_shape_manual(
    name = "Regime",
    values = c(16, 17),
    labels = c(
      "Regime 1",
      "Regime 2"
    )
  ) +
  scale_alpha(
    name = "Probability",
    range = c(0.4, 1)
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

print(plot_msr_projection)

# =============================================================================
# EXACT REGIME-PATH PROBABILITIES: 2019-2023
# =============================================================================

future_paths_raw <- expand.grid(
  s2019 = 1:2,
  s2020 = 1:2,
  s2021 = 1:2,
  s2022 = 1:2,
  s2023 = 1:2
)

exact_future_prob_raw <- apply(
  future_paths_raw,
  1,
  function(path) {

    prob_path <- 0

    for (s2018 in 1:2) {

      prob_conditional <-
        pi_2018_raw[s2018] *
        P_raw[path[1], s2018]

      for (tt in 2:length(path)) {
        prob_conditional <-
          prob_conditional *
          P_raw[
            path[tt],
            path[tt - 1]
          ]
      }

      prob_path <- prob_path + prob_conditional
    }

    prob_path
  }
)

exact_paths_raw <- data.frame(
  Path = apply(
    future_paths_raw,
    1,
    function(x) {
      paste0(
        "R",
        x,
        collapse = "-"
      )
    }
  ),
  Probability = exact_future_prob_raw,
  Percent = 100 * exact_future_prob_raw
)

exact_paths_raw <- exact_paths_raw[
  order(
    exact_paths_raw$Probability,
    decreasing = TRUE
  ),
]

rownames(exact_paths_raw) <- NULL

complete_paths_raw <- data.frame()

for (i in seq_len(nrow(future_paths_raw))) {

  future_path <- as.integer(
    future_paths_raw[i, ]
  )

  for (s2018 in 1:2) {

    complete_path <- c(
      s2018,
      future_path
    )

    prob_complete <- pi_2018_raw[s2018]

    for (tt in 2:length(complete_path)) {
      prob_complete <-
        prob_complete *
        P_raw[
          complete_path[tt],
          complete_path[tt - 1]
        ]
    }

    n_switches <- sum(
      diff(complete_path) != 0
    )

    complete_paths_raw <- rbind(
      complete_paths_raw,
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

complete_paths_raw <- complete_paths_raw[
  order(
    complete_paths_raw$Probability,
    decreasing = TRUE
  ),
]

rownames(complete_paths_raw) <- NULL

exact_classes_raw <- data.frame(
  Path_class = c(
    "No switch",
    "Single switch",
    "Switch and return",
    "More than two switches"
  ),
  Number_of_transitions = c(
    "0",
    "1",
    "2",
    ">2"
  ),
  Probability = c(
    sum(
      complete_paths_raw$Probability[
        complete_paths_raw$Switches == 0
      ]
    ),
    sum(
      complete_paths_raw$Probability[
        complete_paths_raw$Switches == 1
      ]
    ),
    sum(
      complete_paths_raw$Probability[
        complete_paths_raw$Switches == 2
      ]
    ),
    sum(
      complete_paths_raw$Probability[
        complete_paths_raw$Switches > 2
      ]
    )
  )
)

exact_classes_raw$Percent <-
  100 * exact_classes_raw$Probability

print(
  exact_classes_raw,
  row.names = FALSE
)

print(exact_paths_raw)
print(complete_paths_raw)

stopifnot(
  abs(
    sum(exact_paths_raw$Probability) - 1
  ) < 1e-10,
  abs(
    sum(exact_classes_raw$Probability) - 1
  ) < 1e-10
)

# =============================================================================
# TRANSITION-PROBABILITY CONFIDENCE INTERVALS
# =============================================================================

transition_ci_model_indexed <- function(obj) {

  stopifnot(
    nrow(obj@transMat) == 2L,
    ncol(obj@transMat) == 2L
  )

  par0 <- c(
    log(obj@std[1]),
    log(obj@std[2]),
    qlogis(obj@transMat[1, 1]),
    qlogis(obj@transMat[1, 2]),
    obj@Coef[1, 1],
    obj@Coef[2, 1]
  )

  hess_res <- nlme::fdHess(
    pars = par0,
    fun = MSwM:::fopt.lm,
    object = obj
  )

  H <- hess_res$Hessian
  V <- solve(H)

  idx_trans <- 3:4

  eta_hat <- par0[idx_trans]
  se_eta <- sqrt(diag(V))[idx_trans]

  z <- qnorm(0.975)

  lower_eta <- eta_hat - z * se_eta
  upper_eta <- eta_hat + z * se_eta

  p_hat <- plogis(eta_hat)
  p_lower <- plogis(lower_eta)
  p_upper <- plogis(upper_eta)

  p11 <- p_hat[1]
  p12 <- p_hat[2]

  p11_L <- p_lower[1]
  p11_U <- p_upper[1]

  p12_L <- p_lower[2]
  p12_U <- p_upper[2]

  P_hat <- matrix(
    c(
      p11,
      1 - p11,
      p12,
      1 - p12
    ),
    nrow = 2,
    ncol = 2
  )

  P_lower <- matrix(
    c(
      p11_L,
      1 - p11_U,
      p12_L,
      1 - p12_U
    ),
    nrow = 2,
    ncol = 2
  )

  P_upper <- matrix(
    c(
      p11_U,
      1 - p11_L,
      p12_U,
      1 - p12_L
    ),
    nrow = 2,
    ncol = 2
  )

  list(
    estimate = P_hat,
    lower = P_lower,
    upper = P_upper,
    hessian = H,
    covariance = V,
    eigenvalues = eigen(
      H,
      symmetric = TRUE
    )$values,
    condition_number = kappa(H)
  )
}

transition_ci_low_high <- function(
  obj,
  ci_object,
  sample_label
) {

  means <- as.numeric(
    obj@Coef[, 1]
  )

  low <- which.min(means)
  high <- which.max(means)

  transitions <- data.frame(
    Series = "Raw",
    Sample = sample_label,
    Transition = c(
      "Low -> Low",
      "Low -> High",
      "High -> Low",
      "High -> High"
    ),
    Estimate = c(
      ci_object$estimate[low, low],
      ci_object$estimate[high, low],
      ci_object$estimate[low, high],
      ci_object$estimate[high, high]
    ),
    CI_2.5 = c(
      ci_object$lower[low, low],
      ci_object$lower[high, low],
      ci_object$lower[low, high],
      ci_object$lower[high, high]
    ),
    CI_97.5 = c(
      ci_object$upper[low, low],
      ci_object$upper[high, low],
      ci_object$upper[low, high],
      ci_object$upper[high, high]
    )
  )

  transitions
}

transition_ci_train_raw <-
  transition_ci_model_indexed(
    msr_k2_train
  )

transition_ci_full_raw <-
  transition_ci_model_indexed(
    msr_k2_full
  )

transition_CI_raw_train <-
  transition_ci_low_high(
    obj = msr_k2_train,
    ci_object = transition_ci_train_raw,
    sample_label = "Training 1990-2018"
  )

transition_CI_raw_full <-
  transition_ci_low_high(
    obj = msr_k2_full,
    ci_object = transition_ci_full_raw,
    sample_label = "Full 1990-2023"
  )

transition_CI_raw <- rbind(
  transition_CI_raw_train,
  transition_CI_raw_full
)

print(
  transition_CI_raw,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nTraining Hessian minimum eigenvalue:",
  min(transition_ci_train_raw$eigenvalues),
  "\nTraining Hessian condition number:",
  transition_ci_train_raw$condition_number,
  "\n"
)

cat(
  "\nFull-sample Hessian minimum eigenvalue:",
  min(transition_ci_full_raw$eigenvalues),
  "\nFull-sample Hessian condition number:",
  transition_ci_full_raw$condition_number,
  "\n"
)

# =============================================================================
# MSR(2,1) CLASSIFICATION STABILITY
# =============================================================================

if (
  isTRUE(raw_train_k2p1$fit_ok) &&
  isTRUE(raw_full_k2p1$fit_ok)
) {

  prob_train_p1 <-
    raw_train_k2p1$fit@Fit@smoProb[
      -1,
      ,
      drop = FALSE
    ]

  years_train_p1 <- 1991:2018

  stopifnot(
    nrow(prob_train_p1) ==
      length(years_train_p1)
  )

  class_train_p1 <- apply(
    prob_train_p1,
    1,
    which.max
  )

  train_p1 <- data.frame(
    Year = years_train_p1,
    Regime = class_train_p1,
    Prob_Regime1 = prob_train_p1[, 1],
    Prob_Regime2 = prob_train_p1[, 2],
    Dominant_Prob = apply(
      prob_train_p1,
      1,
      max
    )
  )

  prob_full_p1 <-
    raw_full_k2p1$fit@Fit@smoProb[
      -1,
      ,
      drop = FALSE
    ]

  years_full_p1 <- 1991:2023

  stopifnot(
    nrow(prob_full_p1) ==
      length(years_full_p1)
  )

  class_full_p1 <- apply(
    prob_full_p1,
    1,
    which.max
  )

  full_p1 <- data.frame(
    Year = years_full_p1,
    Regime = class_full_p1,
    Prob_Regime1 = prob_full_p1[, 1],
    Prob_Regime2 = prob_full_p1[, 2],
    Dominant_Prob = apply(
      prob_full_p1,
      1,
      max
    )
  )

  trans_train_p1 <-
    train_p1$Year[
      c(
        FALSE,
        diff(train_p1$Regime) != 0
      )
    ]

  trans_full_p1 <-
    full_p1$Year[
      c(
        FALSE,
        diff(full_p1$Regime) != 0
      )
    ]

  comparison_p1 <- merge(
    train_p1[
      ,
      c(
        "Year",
        "Regime",
        "Dominant_Prob"
      )
    ],
    full_p1[
      ,
      c(
        "Year",
        "Regime",
        "Dominant_Prob"
      )
    ],
    by = "Year",
    suffixes = c(
      "_Training",
      "_Full"
    )
  )

  cat(
    "\nMSR(2,1) training transition years:\n"
  )
  print(trans_train_p1)

  cat(
    "\nMSR(2,1) full-sample transition years:\n"
  )
  print(trans_full_p1)

  print(
    comparison_p1,
    row.names = FALSE
  )
}


# =============================================================================
# COMPACT AUDIT OUTPUT FOR REPLICATION CHECK
# =============================================================================

sink("03_raw_msr_AUDIT.txt")

cat("\n============================================================\n")
cat("1. SPECIFICATION STATUS\n")
cat("============================================================\n")
print(msr_status, row.names = FALSE)


cat("\n\n============================================================\n")
cat("2. BASELINE MSR(2,0) - TRAINING 1990-2018\n")
cat("============================================================\n")
print(summary(msr_k2_train))

cat("\nTransition matrix:\n")
print(msr_k2_train@transMat, digits = 10)

cat("\nRegime means:\n")
print(msr_k2_train@Coef, digits = 10)

cat("\nRegime residual SD:\n")
print(msr_k2_train@std, digits = 10)

means_train <- as.numeric(msr_k2_train@Coef[, 1])
cat("\nLOW regime = model regime", which.min(means_train),
    "| mean =", min(means_train), "\n")
cat("HIGH regime = model regime", which.max(means_train),
    "| mean =", max(means_train), "\n")


cat("\n\n============================================================\n")
cat("3. BASELINE MSR(2,0) - FULL 1990-2023\n")
cat("============================================================\n")
print(summary(msr_k2_full))

cat("\nTransition matrix:\n")
print(msr_k2_full@transMat, digits = 10)

cat("\nRegime means:\n")
print(msr_k2_full@Coef, digits = 10)

cat("\nRegime residual SD:\n")
print(msr_k2_full@std, digits = 10)

means_full <- as.numeric(msr_k2_full@Coef[, 1])
cat("\nLOW regime = model regime", which.min(means_full),
    "| mean =", min(means_full), "\n")
cat("HIGH regime = model regime", which.max(means_full),
    "| mean =", max(means_full), "\n")


cat("\n\n============================================================\n")
cat("4. BASELINE REGIME TRANSITION YEARS - LOW/HIGH ALIGNED\n")
cat("============================================================\n")

train_changes <- raw_prob_train$Year[
  c(FALSE, raw_prob_train$Dominant_Regime[-1] !=
      raw_prob_train$Dominant_Regime[-nrow(raw_prob_train)])
]

full_changes <- raw_prob_full$Year[
  c(FALSE, raw_prob_full$Dominant_Regime[-1] !=
      raw_prob_full$Dominant_Regime[-nrow(raw_prob_full)])
]

cat("\nTraining transition years:\n")
print(train_changes)

cat("\nFull-sample transition years:\n")
print(full_changes)


cat("\n\n============================================================\n")
cat("5. SMOOTHED PROBABILITIES AROUND TRANSITIONS\n")
cat("============================================================\n")

focus_years <- c(1994:1996, 2010:2012)

cat("\nTRAINING:\n")
print(
  raw_prob_train[
    raw_prob_train$Year %in% focus_years,
    c("Year", "Prob_Low", "Prob_High",
      "Dominant_Regime", "Dominant_Prob")
  ],
  row.names = FALSE,
  digits = 8
)

cat("\nFULL:\n")
print(
  raw_prob_full[
    raw_prob_full$Year %in% focus_years,
    c("Year", "Prob_Low", "Prob_High",
      "Dominant_Regime", "Dominant_Prob")
  ],
  row.names = FALSE,
  digits = 8
)


cat("\n\n============================================================\n")
cat("6. DETERMINISTIC REGIME PROJECTION 2019-2023\n")
cat("============================================================\n")
print(df_regime_forecast_raw, row.names = FALSE, digits = 8)


cat("\n\n============================================================\n")
cat("7. EXACT REGIME-PATH CLASSES 2019-2023\n")
cat("============================================================\n")
print(exact_classes_raw, row.names = FALSE, digits = 10)

cat("\nSum exact path probabilities:\n")
print(sum(exact_paths_raw$Probability), digits = 12)


cat("\n\n============================================================\n")
cat("8. TRANSITION PROBABILITIES + 95% CI - LOW/HIGH ALIGNED\n")
cat("============================================================\n")
print(transition_CI_raw, row.names = FALSE, digits = 8)


cat("\n\n============================================================\n")
cat("9. MSR(2,1) STABILITY\n")
cat("============================================================\n")

if (exists("trans_train_p1")) {
  cat("\nTraining transition years MSR(2,1):\n")
  print(trans_train_p1)
}

if (exists("trans_full_p1")) {
  cat("\nFull transition years MSR(2,1):\n")
  print(trans_full_p1)
}

if (exists("comparison_p1")) {
  cat("\nClassification agreement on common years:\n")
  print(
    table(
      comparison_p1$Regime_Training ==
        comparison_p1$Regime_Full
    )
  )
}


cat("\n\n============================================================\n")
cat("END AUDIT\n")
cat("============================================================\n")

sink()

sessionInfo()

# End of script
