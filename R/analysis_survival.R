#' Run Survival Analysis and Generate Plots
#'
#' @param data Data frame containing survival data.
#' @param outcome Name of outcome variable (e.g. "PFS").
#' @param time_var Time variable name. 
#' @param event_var Event variable name. Must be either 0/1 (being 0 censored, and 1 event) or 1/2 (being 1 censored and 2 event).
#' @param test_var Predictor variable name.
#' @param covariates Optional covariates vector (can be >1). 
#' @param interaction_var Optional interaction variable  (can be >1).
#' @param ref_level Reference level for test_var. If no level is specified, the reference level will be selected by alphanumerical order. 
#' @param analysis_name Analysis title label.
#' @param filename File path to save output plot.
#' @param conf_int Confidence interval level (default 0.95).
#' @param custom_colors Palette color vector. Default colours are #1B9E77, #D95F02, #7570B3, #E7298A, #E6AB02 and #66A61E.
#' @param custom_linetypes Linetypes for curves as per ggsurvplot function (default is 1, solid).
#' @param break_time_by Step size for x-axis time breaks. If this argument is not specified, it is automatically set. 
#' @param plot Logical; whether to generate and save plot.
#' @return A summary data frame of median survival and hazard ratios.
#' @export

analyze_survival <- function(data, outcome = NA, time_var, event_var, test_var,
                             covariates = c(), interaction_var = NULL,
                             ref_level = NULL, analysis_name = NA, filename = "KM.png",
                             conf_int = 0.95, custom_colors = default_colors,
                             custom_linetypes = 1, break_time_by = NULL, plot = TRUE) {
  
  data[[time_var]]  <- as.numeric(data[[time_var]])
  data[[event_var]] <- as.numeric(data[[event_var]])
  data <- set_ref_level(data, test_var, ref_level)
  for (v in covariates) {
    if (is.character(data[[v]]) || is.factor(data[[v]])) data[[v]] <- as.factor(data[[v]])
  }
  
  model_cols <- c(time_var, event_var, test_var, covariates)
  data <- data[complete.cases(data[, model_cols, drop = FALSE]), ]
  
  if (nrow(data) == 0) {
    message(paste("Skipping analysis:", analysis_name, "- Data is empty (or no complete cases)."))
    return(NULL)
  }

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
  
  if(plot == TRUE){
  
    p <- ggsurvplot(
      fit_km, data = data, risk.table = TRUE, legend.title = "",
      break.time.by = if (is.null(break_time_by)) choose_break_time(max(data[[time_var]], na.rm = TRUE)) else break_time_by,
      fontsize = 3, title = analysis_name,
      ggtheme = theme_classic2(), xlab = "Time (months)", ylab = paste(outcome, "(%)"),
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
    
    if (!grepl("\\.png$", filename, ignore.case = TRUE)) filename <- paste0(filename, ".png")
    png(filename, width = 8, height = 6, units = "in", res = 300)
    on.exit(dev.off(), add = TRUE)
    print(p)
  }
  
  return(summary_results)
  
}
