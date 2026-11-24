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
#' @param covariates Optional character vector of covariate column names for multivariable adjustment. If `test_var` has length > 1 (combined grouping), any covariate that duplicates one of `test_var`'s components is dropped automatically (with a warning), since it is already fully captured by the combined grouping.
#' @param interaction_var Optional character vector for interaction term testing. The variable must be within the "covariates" variables as well. If `test_var` has length > 1 (combined grouping), any interaction_var that duplicates one of `test_var`'s components is dropped automatically (with a warning): keeping it would otherwise give a degenerate interaction p-value of 1.
#' @param ref_level Optional character string specifying reference level for `test_var`.
#' @param analysis_name Title label for the analysis and plot.
#' @param imgname File path/name for saving the KM plot (default "KM.png"). Must have an available extension among png, pdf, jpeg, jpg, tiff, bmp or svg. An unrecognized or missing extension will fall back to ".png".
#' @param filename File path/name for saving the summary table as CSV (default "table.csv"). An extension other than ".csv" (or none) will have ".csv" appended.
#' @param conf_int Numeric confidence level (default 0.95).
#' @param conf_type Character; the confidence interval transformation passed to survival::survfit()'s `conf.type` for both median survival and the plotted KM curves - one of "log-log" (default), "log" (survival::survfit()'s own default when unspecified), "plain", or "logit". Different conf_type values give the same point estimates (median, HR) but different CIs - if you are comparing against a manual survival::survfit() call that did not set conf.type, that call used "log", not this function's default of "log-log".
#' @param custom_colors Vector of hex color codes or color names for plot strata. Default colours are "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#E6AB02" and "#66A61E".
#' @param custom_linetypes Vector of linetypes for plot strata as in ggsurvplot. Default is 1 (solid). If the same linetype is applied to all strata, the specific linetype can only be specified once.
#' @param break_time_by Step size for x-axis time breaks. Automatically calculated if NULL (default).
#' @param plot Logical; whether to draw and export plot (default TRUE).
#' @param width Numeric; width of the plot to save, in inches (default 8).
#' @param height Numeric; height of the plot to save, in inches (default 6).
#' @param res Numeric; resolution of the plot to save (default 300).
#'
#' @param return_details If TRUE, return a list containing the raw summary, cached
#'   plot, annotation grob, group key and presentation arguments. Default FALSE
#'   preserves the data.frame return type and machine-readable column names.
#' @param group_labels Optional display labels, in fitted-stratum order or named
#'   by exact original group values. Never changes factor levels or contrasts.
#' @param variable_labels Optional named mapping from internal variable names to display names.
#' @param legend_font_size Legend text size in points.
#' @param legend_ncol Number of legend columns; NULL or 0 selects automatically.
#' @param legend_wrap_width Maximum approximate characters per legend line.
#' @param risk_table_labels Use "group_id" (default) or short "label" row labels.
#' @param annotation_position Normalised c(x, y) position of the annotation's
#'   top-left corner, measured from the top-left of the ENTIRE exported figure.
#' @param annotation_font_size Statistics-box font size in points.
#' @return A data.frame containing sample sizes (N) and number of events (N_events) per group
#'   (for a numeric `test_var`, the total N and total events of the model), median survival with CIs, 
#'   hazard ratios (HR/aHR) with CIs, log-rank p-values, proportional hazards test 
#'   p-values, P_for_Interaction (test between test_var and an external interaction_var),
#'   and, when `test_var` has length > 1, P_for_TestVar_Interaction - a likelihood-ratio
#'   test comparing an additive model to a fully crossed one, testing for an interaction
#'   among the test_var components themselves (N/A when `test_var` has length 1).
#' @export

MARData_surv_analysis <- function(data, outcome = NA, time_var, event_var, test_var,
                              covariates = c(), interaction_var = NULL,
                              ref_level = NULL, analysis_name = NA, imgname = "KM.png", filename = "table.csv",
                              conf_int = 0.95, conf_type = "log-log", custom_colors = default_colors,
                              custom_linetypes = 1, break_time_by = NULL,
                              plot = TRUE, width = 8, height = 6, res = 300,
                              return_details = FALSE, group_labels = NULL, variable_labels = NULL,
                              legend_font_size = 9, legend_ncol = NULL, legend_wrap_width = 28,
                              risk_table_labels = "group_id",
                              annotation_position = c(x = 0.14, y = 0.42),
                              annotation_font_size = 9) {

  original_test_vars <- test_var

  data[[time_var]]  <- as.numeric(data[[time_var]])
  data[[event_var]] <- as.numeric(data[[event_var]])
  is_combined_test_var <- length(test_var) > 1
  testvar_interaction_p <- "N/A"

  if (is_combined_test_var) {
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

    # Safeguard: a variable that is one of the combined test_var components is
    # fully redundant if it is ALSO listed in covariates/interaction_var - the
    # combined grouping already encodes its value for every row. Leaving it in
    # as well creates perfect collinearity: its own coefficient becomes
    # unidentifiable (NA), and any interaction test involving it compares two
    # models with an identical likelihood, so the LRT degenerates to p = 1
    # instead of erroring. Drop it from covariates/interaction_var and warn,
    # rather than silently returning that degenerate result.
    dup_covariates <- intersect(covariates, test_var)
    if (length(dup_covariates) > 0) {
      warning("Variable(s) also present in the combined test_var, so removed from 'covariates' ",
              "(redundant - already captured by the combined grouping): ",
              paste(dup_covariates, collapse = ", "))
      covariates <- setdiff(covariates, dup_covariates)
    }

    dup_interaction <- intersect(interaction_var, test_var)
    if (length(dup_interaction) > 0) {
      warning("Variable(s) also present in the combined test_var, so removed from 'interaction_var' ",
              "(redundant - already captured by the combined grouping, which would otherwise give a ",
              "degenerate interaction p-value of 1): ", paste(dup_interaction, collapse = ", "))
      interaction_var <- setdiff(interaction_var, dup_interaction)
      if (length(interaction_var) == 0) interaction_var <- NULL
    }

    # Formal interaction test between the test_var components themselves
    # (distinct from P_for_Interaction, which tests the combined test_var
    # against an external interaction_var): compares an additive model to a
    # fully crossed one via LRT. Uses the ORIGINAL component columns, not the
    # combined factor - a saturated combined factor and "a * b" span the same
    # design space (same overall fit), but only "a * b" decomposes that into
    # a testable interaction term rather than K-1 cell-vs-reference contrasts.
    add_formula  <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time_var, event_var,
                                        paste(c(test_var, covariates), collapse = " + ")))
    full_formula <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time_var, event_var,
                                        paste(c(paste(test_var, collapse = " * "), covariates), collapse = " + ")))
    fit_add  <- fit_coxph_checked(add_formula, data)
    fit_full <- fit_coxph_checked(full_formula, data)

    if (is.null(fit_add$fit) || is.null(fit_full$fit)) {
      warning("Could not fit the additive/full models needed to test for an interaction among ",
              "the test_var components; P_for_TestVar_Interaction will be N/A.")
    } else {
      if (!fit_add$converged || !fit_full$converged) {
        warning("Cox model(s) for the test_var interaction test did not fully converge; ",
                "P_for_TestVar_Interaction may be unreliable - interpret with caution.")
      }
      lrt <- tryCatch(stats::anova(fit_add$fit, fit_full$fit, test = "LRT"), error = function(e) NULL)
      p_int <- if (!is.null(lrt)) lrt$`Pr(>|Chi|)`[2] else NA
      if (!is.na(p_int)) {
        testvar_interaction_p <- if (p_int < 0.001) "<0.001" else sprintf("%.4f", p_int)
      }
    }

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
  data <- data[stats::complete.cases(data[, model_cols, drop = FALSE]), ]

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
      N_events            = sum(survival::Surv(data[[time_var]], data[[event_var]])[, "status"]),
      Median_Survival     = NA_character_,
      HR_CI               = df_hrs[[hr_col]],
      P_Value             = df_hrs$P_Value,
      Overall_LogRank_P   = "N/A",
      Variable_PH_p_value = df_hrs$Variable_PH_p_value,
      Global_PH_p_value   = df_hrs$Global_PH_p_value,
      P_for_Interaction   = df_hrs$P_for_Interaction,
      P_for_TestVar_Interaction = testvar_interaction_p,
      Adjusted_For        = if (length(covariates) > 0) paste(covariates, collapse = ", ") else "Unadjusted (univariate)",
      stringsAsFactors    = FALSE
    )
    colnames(summary_results)[colnames(summary_results) == "HR_CI"] <- hr_col
    save_summary_table(summary_results, filename)
    if (isTRUE(return_details)) return(list(
      summary = summary_results, plot = NULL, annotation = NULL, group_key = NULL,
      plot_args = NULL, annotation_args = NULL, main_term = test_var,
      component_vars = original_test_vars, conf_int = conf_int, width = width, height = height,
      annotation_position = annotation_position
    ))
    return(summary_results)
  }

  # ============================================================
  # Categorical test_var (factor/character): existing KM + Cox + plot
  # ============================================================
  km_formula <- stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time_var, event_var, test_var))
  fit_km <- survival::survfit(km_formula, data = data, conf.type = conf_type, conf.int = conf_int)
  fit_km$call$formula <- km_formula

  df_medians <- get_km_medians(km_formula, data, conf_level = conf_int, conf_type = conf_type)

  # Events per group, read from the same survfit table that provides N/median
  # (so it follows survival's own event coding: 0/1 or 1/2). Attached by
  # matching on Group - not merged - so the row order of df_medians is kept.
  km_tab <- summary(fit_km)$table
  if (is.null(dim(km_tab))) {
    km_tab <- matrix(km_tab, nrow = 1, dimnames = list("Overall", names(km_tab)))
  }
  km_groups <- gsub("^[^=]*=", "", rownames(km_tab))
  df_medians$N_events <- as.vector(km_tab[, "events"])[match(df_medians$Group, km_groups)]

  df_hrs <- get_cox_hrs(data, time_var, event_var, test_var, covariates, interaction_var,
                        ref_level = ref_level, conf_level = conf_int)
  if (is.null(df_hrs)) return(NULL)

  df_stats <- merge(df_medians, df_hrs, by = "Group", all.y = TRUE)

  n_strata <- length(fit_km$strata)
  if (n_strata > 1) {
    surv_diff <- survival::survdiff(km_formula, data = data)
    overall_p <- stats::pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
    overall_p_str <- if (overall_p < 0.001) "<0.001" else sprintf("%.4f", overall_p)
  } else {
    overall_p_str <- "N/A"
  }

  med_col <- grep("^Median_Survival", colnames(df_stats), value = TRUE)[1]
  hr_col  <- grep("^a?HR_.*CI$", colnames(df_stats), value = TRUE)[1]

  covariate_rows <- !is.na(df_stats$Term) & df_stats$Term != test_var
  df_stats[covariate_rows, c("N", "N_events", med_col)] <- NA

  summary_results <- data.frame(
    Scenario           = analysis_name,
    Outcome            = outcome,
    Variable_Term      = df_stats$Term,
    Group              = df_stats$Group,
    N                  = df_stats$N,
    N_events           = df_stats$N_events,
    Median_Survival    = df_stats[[med_col]],
    HR_CI              = df_stats[[hr_col]],
    P_Value            = df_stats$P_Value,
    Overall_LogRank_P  = overall_p_str,
    Variable_PH_p_value = df_stats$Variable_PH_p_value,
    Global_PH_p_value  = df_stats$Global_PH_p_value,
    P_for_Interaction  = df_stats$P_for_Interaction,
    P_for_TestVar_Interaction = testvar_interaction_p,
    Adjusted_For       = if (length(covariates) > 0) paste(covariates, collapse = ", ") else "Unadjusted (univariate)",
    stringsAsFactors   = FALSE
  )
  colnames(summary_results)[colnames(summary_results) == "HR_CI"] <- hr_col

  p <- NULL
  key <- NULL
  annotation <- NULL
  plot_args <- NULL
  annotation_args <- NULL
  if (isTRUE(plot) || isTRUE(return_details)) {
    key <- make_survival_group_key(fit_km, data, test_var, original_test_vars,
                                  group_labels, variable_labels)
    plot_args <- list(
      fit_km = fit_km, data = data, group_key = key,
      analysis_name = analysis_name, outcome = outcome,
      custom_colors = custom_colors, custom_linetypes = custom_linetypes,
      break_time_by = break_time_by, time_var = time_var,
      legend_font_size = legend_font_size, legend_ncol = legend_ncol,
      legend_wrap_width = legend_wrap_width, risk_table_labels = risk_table_labels,
      width = width, height = height
    )
    p <- do.call(build_survival_plot, plot_args)
    annotation_args <- list(
      cox_df = df_hrs, test_var = test_var, group_key = key, overall_p = overall_p_str,
      conf_int = conf_int, component_interaction_p = testvar_interaction_p,
      font_size = annotation_font_size
    )
    annotation <- do.call(make_survival_annotation, annotation_args)
    annotation_position <- normalise_annotation_position(annotation_position, annotation, width, height)
    if (isTRUE(plot)) save_survival_plot(p, imgname, width = width, height = height, res = res,
                                       annotation = annotation, annotation_position = annotation_position)
  }

  save_summary_table(summary_results, filename)
  if (isTRUE(return_details)) return(list(
    summary = summary_results, plot = p, annotation = annotation, group_key = key,
    plot_args = plot_args, annotation_args = annotation_args, main_term = test_var,
    component_vars = original_test_vars, conf_int = conf_int, width = width, height = height,
    annotation_position = annotation_position
  ))
  return(summary_results)
}
