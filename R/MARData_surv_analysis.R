#' Run Automated Survival Analysis, KM Plotting, and Cox Modeling
#'
#' @description
#' Generates Kaplan-Meier survival curves, fits Cox proportional hazards 
#' models (univariate or multivariable), checks proportional hazards 
#' assumptions via Schoenfeld residuals, and outputs summary statistics.
#'
#' @param data A data.frame containing clinical/survival data.
#' @param outcome Character string describing the outcome (e.g., "Overall Survival").
#' @param time_var Character string specifying column name for time.
#' @param event_var Character string specifying column name for event/status (0/1 - where 0 is censored and 1 event - or 1/2 - where 1 is censored and 2 is event -).
#' @param test_var Character string specifying primary grouping variable(s) to test.
#' @param covariates Optional character vector of covariate column names for multivariable adjustment.
#' @param interaction_var Optional character vector for interaction term testing. The variable must be within the "covariates" variables as well. 
#' @param ref_level Optional character string specifying reference level for `test_var`.
#' @param analysis_name Title label for the analysis and plot.
#' @param filename File path/name for saving the KM plot (default "KM.png"). Must have an available extension among png, pdf, jpeg, jpg, tiff, bmp or svg. An unrecognized or missing extension will fall back to ".png".
#' @param conf_int Numeric confidence level (default 0.95).
#' @param custom_colors Vector of hex color codes or color names for plot strata. Default colours are "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#E6AB02" and "#66A61E".
#' @param custom_linetypes Vector of linetypes for plot strata as in ggsurvplot. Default is 1 (solid). If the same linetype is applied to all strata, the specific linetype can only be specified once.
#' @param break_time_by Step size for x-axis time breaks. Automatically calculated if NULL (default).
#' @param plot Logical; whether to draw and export plot (default TRUE).
#' @param width Numeric; width of the plot to save, in inches (default 8).
#' @param height Numeric; height of the plot to save, in inches (default 6).
#' @param res Numeric; resolution of the plot to save (default 300).
#'
#' @return A data.frame containing sample sizes, median survival with CIs, 
#'   hazard ratios (HR/aHR) with CIs, log-rank p-values, proportional hazards test 
#'   p-values, and interaction p-values.
#' @export

MARData_surv_analysis <- function(data, outcome = NA, time_var, event_var, test_var,
                              covariates = c(), interaction_var = NULL,
                              ref_level = NULL, analysis_name = NA, filename = "KM.png",
                              conf_int = 0.95, custom_colors = default_colors,
                              custom_linetypes = 1, break_time_by = NULL,
                              plot = TRUE, width = 8, height = 6, res = 300) {

  data[[time_var]]  <- as.numeric(data[[time_var]])
  data[[event_var]] <- as.numeric(data[[event_var]])
  is_combined_test_var <- length(test_var) > 1

  if (length(test_var) > 1) {
    for (v in test_var) {
      if (is.numeric(data[[v]]) && length(unique(stats::na.omit(data[[v]]))) > 10) {
        warning("'", v, "' looks continuous (more than 10 unique values) and is being combined ",
                "as a categorical grouping variable; consider discretizing it first if that is not intended.")
      }
    }
    combo_name <- paste(test_var, collapse = "_x_")
    while (combo_name %in% names(data)) combo_name <- paste0(combo_name, "_")

    labeled_parts <- lapply(test_var, function(v) paste(v, as.character(data[[v]])))
    combo_vals <- do.call(paste, c(labeled_parts, list(sep = " + ")))
    # A row with a missing value in ANY of the component variables must stay NA in the combined column too
    any_na <- Reduce(`|`, lapply(test_var, function(v) is.na(data[[v]])))
    combo_vals[any_na] <- NA

    data[[combo_name]] <- droplevels(factor(combo_vals))
    test_var <- combo_name
  }

  # A numeric test_var is treated as continuous: no KM curves, no median survival, and no plot - only a Cox model with test_var as a continuous predictor.
  test_var_is_numeric <- is.numeric(data[[test_var]])

  if (test_var_is_numeric) {
    if (!is.null(ref_level)) {
      warning("ref_level is ignored because '", test_var, "' is numeric (continuous); ",
              "reference levels only apply to categorical variables.")
    }
  } else {
    data <- set_ref_level(data, test_var, ref_level)
  }

  for (v in covariates) {
    if (is.character(data[[v]]) || is.factor(data[[v]])) data[[v]] <- as.factor(data[[v]])
  }

  model_cols <- c(time_var, event_var, test_var, covariates)
  data <- data[complete.cases(data[, model_cols, drop = FALSE]), ]

  if (nrow(data) == 0) {
    message(paste("Skipping analysis:", analysis_name, "- Data is empty (or no complete cases)."))
    return(NULL)
  }

  # ============================================================
  # Numeric test_var: Cox model only. 
  # ============================================================
  if (test_var_is_numeric) {
    message("'", test_var, "' is numeric: fitting a Cox model only - ",
            "Kaplan-Meier curves, median survival, and the plot do not apply to a continuous variable.")

    df_hrs <- get_cox_hrs(data, time_var, event_var, test_var, covariates, interaction_var,
                          ref_level = NULL, conf_level = conf_int)
    if (is.null(df_hrs)) return(NULL)

    hr_col <- grep("^a?HR_.*CI$", colnames(df_hrs), value = TRUE)[1]

    summary_results <- data.frame(
      Scenario            = analysis_name,
      Outcome             = outcome,
      Variable_Term       = df_hrs$Term,
      Group               = df_hrs$Group,
      N                   = nrow(data),
      Median_Survival     = NA_character_,
      HR_CI               = df_hrs[[hr_col]],
      P_Value             = df_hrs$P_Value,
      Overall_LogRank_P   = "N/A",
      Variable_PH_p_value = df_hrs$Variable_PH_p_value,
      Global_PH_p_value   = df_hrs$Global_PH_p_value,
      P_for_Interaction   = df_hrs$P_for_Interaction,
      Adjusted_For        = if (length(covariates) > 0) paste(covariates, collapse = ", ") else "Unadjusted (univariate)",
      stringsAsFactors    = FALSE
    )
    colnames(summary_results)[7] <- hr_col
    return(summary_results)
  }

  # ============================================================
  # Categorical test_var (factor/character): existing KM + Cox + plot
  # ============================================================
  km_formula <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, event_var, test_var))
  fit_km <- survfit(km_formula, data = data, conf.type = "log-log", conf.int = conf_int)
  fit_km$call$formula <- km_formula

  df_medians <- get_km_medians(km_formula, data, conf_level = conf_int)
  df_hrs <- get_cox_hrs(data, time_var, event_var, test_var, covariates, interaction_var,
                        ref_level = ref_level, conf_level = conf_int)
  if (is.null(df_hrs)) return(NULL)

  df_stats <- merge(df_medians, df_hrs, by = "Group", all.y = TRUE)

  n_strata <- length(fit_km$strata)
  if (n_strata > 1) {
    surv_diff <- survdiff(km_formula, data = data)
    overall_p <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
    overall_p_str <- if (overall_p < 0.001) "<0.001" else sprintf("%.4f", overall_p)
  } else {
    overall_p_str <- "N/A"
  }

  med_col <- grep("^Median_Survival", colnames(df_stats), value = TRUE)[1]
  hr_col  <- grep("^a?HR_.*CI$", colnames(df_stats), value = TRUE)[1]

  summary_results <- data.frame(
    Scenario           = analysis_name,
    Outcome            = outcome,
    Variable_Term      = df_stats$Term,
    Group              = df_stats$Group,
    N                  = df_stats$N,
    Median_Survival    = df_stats[[med_col]],
    HR_CI              = df_stats[[hr_col]],
    P_Value            = df_stats$P_Value,
    Overall_LogRank_P  = overall_p_str,
    Variable_PH_p_value = df_stats$Variable_PH_p_value,
    Global_PH_p_value  = df_stats$Global_PH_p_value,
    P_for_Interaction  = df_stats$P_for_Interaction,
    Adjusted_For       = if (length(covariates) > 0) paste(covariates, collapse = ", ") else "Unadjusted (univariate)",
    stringsAsFactors   = FALSE
  )
  colnames(summary_results)[7] <- hr_col

  if (plot) {
    n_curves <- max(n_strata, 1)
    if (n_curves > length(custom_colors)) {
      custom_colors <- grDevices::colorRampPalette(custom_colors)(n_curves)
    }
    if (length(custom_linetypes) == 1) {
      custom_linetypes <- rep(custom_linetypes, n_curves)
    } else if (length(custom_linetypes) < n_curves) {
      message("The number of strata is ", n_curves, " but the number of linetypes set is ", length(custom_linetypes))
      return(NULL)
    }

    p <- ggsurvplot(
      fit_km, data = data, risk.table = TRUE, legend.title = "",
      legend.labs = if (is_combined_test_var) levels(droplevels(as.factor(data[[test_var]]))) else NULL,
      break.time.by = if (is.null(break_time_by)) choose_break_time(max(data[[time_var]], na.rm = TRUE)) else break_time_by,
      fontsize = 3, title = analysis_name,
      ggtheme = theme_publish(), xlab = "Time (months)", ylab = paste(outcome, "(%)"),
      palette = custom_colors, linetype = custom_linetypes, legend = c(0.7, 0.9),
      linewidth = 1, surv.median.line = "hv", risk.table.height = 0.15,
      tables.theme = clean_theme(), break.y.by = 0.1, surv.scale = "percent", pval = FALSE
    )

    annot_text <- generate_surv_annotation(km_formula, data, df_hrs, test_var, conf_level = conf_int)
    if (!is.null(annot_text)) {
      p$plot <- p$plot + ggplot2::annotate(
        "text", x = 1, y = 0.15, label = annot_text,
        hjust = 0, vjust = 0, size = 3.5, fontface = "italic"
      )
    }
    
    save_survival_plot(p, filename, width = width, height = height, res = res)

  }

  return(summary_results)
}

