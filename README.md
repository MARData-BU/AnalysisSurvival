# MARData_surv_analysis

`MARData_surv_analysis` provides a standardized, publication-ready workflow for clinical survival analysis in R. It automates Kaplan-Meier curve generation, hazard ratio calculations via Cox proportional hazards models, proportional hazards assumption testing, and custom plot formatting.

# MARData_surv_cutpoint

`MARData_surv_cutpoint` provides a standardized, publication-ready workflow for clinical survival analysis in R, determining the optimal cutpoint for continuous variables. It automates Kaplan-Meier curve generation, hazard ratio calculations via Cox proportional hazards models, proportional hazards assumption testing, and custom plot formatting.

---

## Features

- **Automated KM Plots:** Generates publication-ready Kaplan-Meier curves (`.png` output) powered by `survminer` and `ggplot2` with customized axis breaks and in-plot statistical annotations.
- **Unified Cox Modeling:** Performs both univariate and multivariable (adjusted) Cox proportional hazards regressions in a single function call.
- **Proportional Hazards Diagnostic:** Automatically evaluates the proportional hazards assumption using Schoenfeld residuals (`cox.zph()`) for each variable and globally.
- **Interaction Testing:** Optional automated Likelihood Ratio Tests (LRT) for specified interaction terms.
- **Structured Summary Output:** Returns a clean, consolidated `data.frame` containing:
  - Patient sample size ($N$) per stratum
  - Median survival times with 95% Confidence Intervals (or any other CI specified)
  - Hazard Ratios (HR / aHR) with $95\%$ CIs (or any other CI specified) and $p$-values
  - Overall log-rank test $p$-value
  - Proportional hazards assumption test $p$-values
  - Interaction $p$-values (if specified)

## Examples:

# Univariate model (strata > 1)
MARData_surv_analysis(data = data, 
  outcome = "PFS", 
  time_var = "pfs_time", 
  event_var = "pfs_event", 
  test_var = "Sex", 
  ref_level = "Male", 
  analysis_name = "Sex", 
  filename = file.path(resultsDir, "Male_vs_Female.png"), 
  covariates = NULL, 
  interaction_var = NULL,
  conf_int = 0.90)

# Multivariate model (without interaction)
MARData_surv_analysis(data = data, 
  outcome = "PFS", 
  time_var = "pfs_time", 
  event_var = "pfs_event", 
  test_var = "Sex", 
  ref_level = "Male", 
  analysis_name = "Sex_multi_Age", 
  filename = file.path(resultsDir, "Male_vs_Female.png"), 
  covariates = "Age", 
  interaction_var = NULL,
  conf_int = 0.90)

# Multivariate model (with interaction)
MARData_surv_analysis(data = data, 
  outcome = "PFS", 
  time_var = "pfs_time", 
  event_var = "pfs_event", 
  test_var = "Sex", 
  ref_level = "Male", 
  analysis_name = "Sex_interact_Age", 
  filename = file.path(resultsDir, "Male_vs_Female.png"), 
  covariates = "Age", 
  interaction_var = "Age",
  conf_int = 0.90)

# Univariate model (strata == 1)
MARData_surv_analysis(data = data[which(data$Sex == "Female"),], 
  outcome = "PFS", 
  time_var = "pfs_time", 
  event_var = "pfs_event", 
  test_var = "Sex", 
  analysis_name = "Sex", 
  filename = file.path(resultsDir, "Female_distribution.png"), 
  covariates = NULL, 
  interaction_var = NULL,
  conf_int = 0.90)

# Optimal cutpoint
MARData_surv_cutpoiont(data = data, 
  outcome = "PFS", 
  time_var = "pfs_time", 
  event_var = "pfs_event", 
  test_var = "Age", 
  analysis_name = "Age_optimal", 
  filename = file.path(resultsDir, "Age_optimal_categorization.png"), 
  covariates = NULL, 
  interaction_var = NULL,
  conf_int = 0.90)
