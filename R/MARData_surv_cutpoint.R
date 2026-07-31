#' Run Automated Survival Analysis, KM Plotting, and Cox Modeling
#'
#' @description
#' Dichotomizes a numerical variable using surv_cutpoint and surv_categorize,
#' generates Kaplan-Meier survival curves, fits Cox proportional hazards 
#' models (univariate or multivariable), checks proportional hazards 
#' assumptions via Schoenfeld residuals, and outputs summary statistics.
#'
#' @param data A data.frame containing clinical/survival data.
#' @param outcome Character string describing the outcome (e.g., "Overall Survival").
#' @param time_var Character string specifying column name for time.
#' @param event_var Character string specifying column name for event/status (0/1 - where 0 is censored and 1 event - or 1/2 - where 1 is censored and 2 is event -).
#' @param test_var Character string specifying primary grouping numerical variable to test.
#' @param covariates Optional character vector of covariate column names for multivariable adjustment.
#' @param interaction_var Optional character vector for interaction term testing. The variable must be within the "covariates" variables as well. 
#' @param ref_level Optional character string specifying reference level for `test_var`.
#' @param analysis_name Title label for the analysis and plot.
#' @param filename File path/name for saving the KM plot (default "KM.png"). Must have an available extension among png, pdf, jpeg, jpg, tiff, bmp or svg. An unrecognized or missing extension will fall back to ".png".
#' @param conf_int Numeric confidence level (default 0.95).
#' @param custom_colors Vector of hex color codes or color names for plot strata. Default colours are "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#E6AB02" and "#66A61E".
#' @param custom_linetypes Vector of linetypes for plot strata as in ggsurvplot. Default is 1 (solid). If the same linetype is applied to all strata, the specific linetype can only be specified once.
#' @param break_time_by Step size for x-axis time breaks. Automatically calculated if NULL (default).
#' @param minprop Numeric the minimal proportion of observations per group (default 0.1).
#' @param plot Logical; whether to draw and export plot (default TRUE).
#' @param width Numeric; width of the plot to save, in inches (default 8).
#' @param height Numeric; height of the plot to save, in inches (default 6).
#' @param res Numeric; resolution of the plot to save (default 300).
#'
#' @return A data.frame containing sample sizes, median survival with CIs, 
#'   hazard ratios (HR/aHR) with CIs, log-rank p-values, proportional hazards test 
#'   p-values, and interaction p-values.
#' @export

MARData_surv_cutpoint <- function(data, outcome = NA, time_var, event_var, test_var,
                                       covariates = c(), interaction_var = NULL, ref_level = NULL,
                                       analysis_name = NA, filename = paste0(test_var, "_optimal_cutpoint.png"),
                                       conf_int = 0.95, custom_colors = default_colors,
                                       custom_linetypes = 1, break_time_by = NULL, minprop = 0.1,
                                       plot = TRUE, width = 8, height = 6, res = 300) {

  # Check numeric-ness on the ORIGINAL column, before any coercion. Coercing
  # first and checking is.numeric() after (as in the draft) can never catch
  # a genuinely non-numeric test_var: as.numeric() on a character/factor
  # column still leaves the column typed "numeric" (just full of NA), so
  # is.numeric() would return TRUE regardless and silently hand surv_cutpoint()
  # a garbage all-NA variable instead of stopping here with a clear message.
  test_var_is_numeric <- is.numeric(data[[test_var]])
  if (!test_var_is_numeric) {
    warning(test_var, " is not numerical, so the optimal cutpoint cannot be computed.")
    return(NULL)
  }

  data[[time_var]]  <- as.numeric(data[[time_var]])
  data[[event_var]] <- as.numeric(data[[event_var]])

  model_cols <- c(time_var, event_var, test_var, covariates)
  data <- data[complete.cases(data[, model_cols, drop = FALSE]), ]

  if (nrow(data) == 0) {
    message(paste("Skipping analysis:", analysis_name, "- Data is empty (or no complete cases)."))
    return(NULL)
  }

  categ_col <- paste0(test_var, "_categ")
  while (categ_col %in% names(data)) categ_col <- paste0(categ_col, "_")

  cut <- surv_cutpoint(data, time = time_var, event = event_var, variables = test_var, minprop = minprop)
  cutpoint <- cut$cutpoint$cutpoint
  cat_df <- surv_categorize(cut)

  data[[categ_col]] <- factor(cat_df[[test_var]], levels = c("low", "high"))

  # ============================================================
  # KM + Cox pipeline on the categorized column - identical to the
  # categorical branch of analyze_survival(), just called on categ_col
  # (the low/high grouping) instead of test_var (the raw numeric column).
  # ============================================================
  data <- set_ref_level(data, categ_col, ref_level)

  km_formula <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, event_var, categ_col))
  fit_km <- survfit(km_formula, data = data, conf.type = "log-log", conf.int = conf_int)
  fit_km$call$formula <- km_formula

  df_medians <- get_km_medians(km_formula, data, conf_level = conf_int)
  df_hrs <- get_cox_hrs(data, time_var, event_var, categ_col, covariates, interaction_var,
                        ref_level = NULL, conf_level = conf_int)
  if (is.null(df_hrs)) return(NULL)

  # Both df_medians and df_hrs now key their "Group" on the SAME low/high
  # factor, so this merge lines up cleanly - all.y = TRUE (matching
  # analyze_survival()) keeps any covariate rows too, which have no KM
  # counterpart and legitimately get NA for N/median.
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
    Optimal_Cutpoint   = round(cutpoint, 2),
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
      break.time.by = if (is.null(break_time_by)) choose_break_time(max(data[[time_var]], na.rm = TRUE)) else break_time_by,
      fontsize = 3, title = analysis_name,
      ggtheme = theme_publish(), xlab = "Time (months)", ylab = paste(outcome, "(%)"),
      palette = custom_colors, linetype = custom_linetypes, legend = c(0.7, 0.9),
      linewidth = 1, surv.median.line = "hv", risk.table.height = 0.15,
      tables.theme = clean_theme(), break.y.by = 0.1, surv.scale = "percent", pval = FALSE
    )

    annot_text <- generate_surv_annotation(km_formula, data, df_hrs, categ_col, conf_level = conf_int)
    if (!is.null(annot_text)) {
      p$plot <- p$plot + ggplot2::annotate(
        "text", x = 1, y = 0.15, label = annot_text,
        hjust = 0, vjust = 0, size = 3.5, fontface = "italic"
      )
    }

    # Same filename-extension -> device logic as analyze_survival().
    ext <- tolower(tools::file_ext(filename))
    raster_devices <- list(png = png, jpeg = jpeg, jpg = jpeg, tiff = tiff, bmp = bmp)
    if (!(ext %in% c(names(raster_devices), "pdf", "svg"))) {
      filename <- paste0(filename, ".png")
      ext <- "png"
    }

    if (ext == "pdf") {
      pdf(filename, width = width, height = height)
    } else if (ext == "svg") {
      svg(filename, width = width, height = height)
    } else {
      raster_devices[[ext]](filename, width = width, height = height, units = "in", res = res)
    }
    on.exit(dev.off(), add = TRUE)
    print(p)
  }

  return(summary_results)
}
