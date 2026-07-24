# AnalysisSurvival

`AnalysisSurvival` provides a standardized, publication-ready workflow for clinical survival analysis in R. It automates Kaplan-Meier curve generation, hazard ratio calculations via Cox proportional hazards models, proportional hazards assumption testing, and custom plot formatting.

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
