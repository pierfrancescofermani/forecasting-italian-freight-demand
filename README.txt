REPLICATION FILES FOR
"Forecasting Italian Freight Demand in the Presence of Regimes"
================================================================

This repository contains the data and R scripts used for the empirical
analysis in the paper "Forecasting Italian Freight Demand in the Presence
of Regimes".

The analysis uses annual Italian freight activity over 1990-2023.
The main estimation sample is 1990-2018, while 2019-2023 is used as the
out-of-sample holdout period.


DATA
====

The repository includes two Excel files:

- data/Total_Tonkm.xlsx
  Main freight-activity dataset used throughout the analysis.
  Variables:
  - Year
  - Total

- data/Tkm_Gdp.xlsx
  Freight activity and GDP dataset used only for the GDP-augmented
  ARIMAX robustness benchmark.
  Variables:
  - Year
  - Tkm
  - GdP

The R scripts assume that the relevant Excel file has already been imported into the R session using the object name shown below.

library(readxl)

Total_Tonkm <- read_excel("data/Total_Tonkm.xlsx")
Tkm_Gdp <- read_excel("data/Tkm_Gdp.xlsx")

Only Total_Tonkm is required for scripts 01–07 and 09.
Tkm_Gdp is required for script 08.


R SCRIPTS
=========

01_raw_series_diagnostics.R
Purpose:
  Descriptive statistics and structural diagnostics for the raw freight
  series, including ACF/PACF, stationarity, serial dependence, ARCH,
  Bai-Perron breakpoints, and CUSUM.

Paper:
  Section 3; Figures 1, 6-7; Tables 1-2; raw-series breakpoint
  confidence intervals.


02_raw_arima_models.R
Purpose:
  ARIMA model selection and residual diagnostics on the raw training
  series, out-of-sample forecasts for 2019-2023, and differencing
  sensitivity analysis.

Paper:
  Section 5.1.1; Table 3; Appendix B, Tables 10-12 and Figures 8-9.


03_raw_msr.R
Purpose:
  Markov-switching analysis of the raw series, including training- and
  full-sample estimation, regime classification, deterministic regime
  projection, exact regime-path probabilities, and transition-probability
  uncertainty.

Paper:
  Sections 5.1.2-5.1.5; Figure 2; Appendix B, Figures 10-11 and
  related MSR tables.


04_loess_screening.R
Purpose:
  LOESS structural screening of the training series and diagnostics for
  the selected LOESS trend and residuals.

Paper:
  Section 5.2.1; Figure 3; Appendix B, Tables 15-16.


05_kalman_trend.R
Purpose:
  Kalman local-linear-trend extraction and calibration diagnostics,
  including automatic calibration tests, the baseline calibration,
  training/full-sample diagnostics, and the constrained-ML Test C
  comparison.

Paper:
  Section 5.2.2; Figure 4; Appendix A; Appendix B, Kalman trend and
  residual diagnostics.


06_kalman_arima.R
Purpose:
  ARIMA estimation and forecasting on the baseline Kalman-smoothed
  training trend.

Paper:
  Section 5.3.1; Table 5; Figure 5; Appendix B, Tables 20-21.


07_kalman_msr.R
Purpose:
  Markov-switching analysis on the baseline Kalman-smoothed trend,
  including deterministic forecasts, exact regime-path probabilities,
  transition uncertainty, and training/full-sample structural comparison.

Paper:
  Sections 5.3.2-5.3.5; Appendix B, Figure 13 and related MSR tables.


08_gdp_arimax_benchmark.R
Purpose:
  Auxiliary GDP-augmented ARIMAX robustness analysis and comparison
  with the univariate ARIMA benchmark.

Paper:
  Section 5.1.1; Table 4.


09_kalman_sensitivity.R
Purpose:
  Downstream sensitivity analysis comparing the baseline Kalman
  calibration with constrained-ML Test C for the extracted trend,
  MSR results, and ARIMA forecasts.

Paper:
  Section 5.4; Tables 6-8.




RUNNING THE ANALYSIS
====================
1. Set the repository root as the working directory.
2. Import the required dataset as shown above.
3. Run the desired R script.

The scripts are self-contained once the corresponding data object has been loaded and can therefore be run individually.
Each script checks whether its required R packages are installed and stops with an informative message if any dependency is missing.

The scripts follow the empirical sequence used in the paper:
01 → 02/03 → 04 → 05 → 06/07 → 09

Script 08_gdp_arimax_benchmark.R is an auxiliary robustness analysis and is not part of the main sequential modelling chain.


ADDITIONAL OUTPUT
=================
03_raw_msr.R also creates:
03_raw_msr_AUDIT.txt
This file provides a compact summary of the main MSR replication outputs and model checks.