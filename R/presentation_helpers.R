# Presentation helpers. Model names, factor levels and contrasts are NOT renamed.

survival_display_name <- function(x) {
  trimws(gsub("[[:space:]]+", " ", gsub("[._]+", " ", as.character(x))))
}

survival_wrap <- function(x, width = 28L) {
  width <- max(8L, as.integer(width))
  vapply(as.character(x), function(s) {
    if (is.na(s)) return("")
    # Also break a single very long token; strwrap() alone does not do this.
    words <- strsplit(s, "[[:space:]]+")[[1L]]
    words <- unlist(lapply(words, function(w) {
      if (nchar(w) <= width) return(w)
      starts <- seq.int(1L, nchar(w), by = width)
      substring(w, starts, pmin(starts + width - 1L, nchar(w)))
    }), use.names = FALSE)
    paste(strwrap(paste(words, collapse = " "), width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

make_survival_group_key <- function(fit, data, test_var, component_vars,
                                    group_labels = NULL, variable_labels = NULL) {
  raw <- if (is.null(fit$strata)) "Overall" else sub("^[^=]*=", "", names(fit$strata))
  key <- data.frame(id = paste0("G", seq_along(raw)), raw = raw,
                    label = raw, definition = raw, stringsAsFactors = FALSE)
  for (i in seq_along(raw)) {
    row <- if (is.null(fit$strata)) 1L else match(raw[i], as.character(data[[test_var]]))
    if (is.na(row)) next
    values <- vapply(component_vars, function(v) as.character(data[[v]][row]), character(1))
    vlabels <- vapply(component_vars, function(v) {
      if (!is.null(variable_labels) && v %in% names(variable_labels)) variable_labels[[v]]
      else survival_display_name(v)
    }, character(1))
    key$definition[i] <- paste(paste0(vlabels, " = ", values), collapse = "; ")
    key$label[i] <- if (length(component_vars) > 1L) paste(values, collapse = " / ") else values[1L]
  }
  # Codebook-style names stay in the key; very long default value labels are shortened.
  key$label <- ifelse(nchar(key$label) > 56L, paste0(substr(key$label, 1L, 53L), "..."), key$label)
  if (!is.null(group_labels)) {
    if (is.null(names(group_labels))) {
      if (length(group_labels) != nrow(key)) stop("Supply one display label per fitted group.")
      key$label <- as.character(group_labels)
    } else {
      idx <- match(key$raw, names(group_labels))
      use <- !is.na(idx)
      key$label[use] <- as.character(group_labels[idx[use]])
    }
  }
  if (anyNA(key$label) || any(!nzchar(trimws(key$label)))) stop("Group labels cannot be empty.")
  key
}

build_survival_plot <- function(fit_km, data, group_key, analysis_name, outcome,
                                custom_colors = default_colors, custom_linetypes = 1,
                                break_time_by = NULL, time_var,
                                legend_font_size = 9, legend_ncol = NULL,
                                legend_wrap_width = 28, risk_table_labels = "group_id",
                                width = 9, height = 7) {
  n <- nrow(group_key)
  if (n < 1L) stop("There are no fitted groups to plot.")
  risk_table_labels <- match.arg(risk_table_labels, c("group_id", "label"))
  if (!length(custom_colors)) stop("At least one plot colour is required.")
  if (length(custom_colors) < n) custom_colors <- grDevices::colorRampPalette(custom_colors)(n)
  if (length(custom_linetypes) == 1L) custom_linetypes <- rep(custom_linetypes, n)
  if (length(custom_linetypes) < n) stop("Supply one linetype, or one per fitted group.")
  legend_labs <- survival_wrap(paste0(group_key$id, ": ", group_key$label), legend_wrap_width)
  max_chars <- max(nchar(unlist(strsplit(legend_labs, "\n", fixed = TRUE))))
  if (is.null(legend_ncol) || legend_ncol == 0) {
    legend_ncol <- max(1L, min(n, 3L, floor(width * 9 / (max_chars + 7))))
  }
  legend_ncol <- max(1L, min(n, as.integer(legend_ncol)))
  risk_text <- paste0(group_key$id, ": ", group_key$label)
  risk_text <- ifelse(nchar(risk_text) > 32L, paste0(substr(risk_text, 1L, 29L), "..."), risk_text)
  risk_labels <- if (risk_table_labels == "group_id") group_key$id else survival_wrap(risk_text, 18)
  risk_lines <- max(lengths(strsplit(risk_labels, "\n", fixed = TRUE)))
  table_height <- min(0.45, max(0.18, (0.55 + 0.23 * n * risk_lines) / height))
  title <- if (length(analysis_name) && !is.na(analysis_name)) survival_wrap(analysis_name, 70) else NULL
  p <- survminer::ggsurvplot(
    fit_km, data = data, risk.table = TRUE, legend.title = "",
    legend.labs = legend_labs, legend = "top",
    break.time.by = if (is.null(break_time_by)) choose_break_time(max(data[[time_var]], na.rm = TRUE)) else break_time_by,
    fontsize = 3, title = title, ggtheme = ggplot2::theme_classic(base_size = 11),
    xlab = "Time (months)", ylab = paste(outcome, "(%)"),
    palette = custom_colors, linetype = custom_linetypes,
    size = 0.85, surv.median.line = "hv", risk.table.height = table_height,
    risk.table.title = "Number at risk", risk.table.y.text.col = FALSE,
    tables.theme = survminer::theme_cleantable(),
    break.y.by = 0.1, surv.scale = "percent", pval = FALSE
  )
  guide <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::guide_legend(ncol = legend_ncol, theme = ggplot2::theme(legend.byrow = TRUE))
  } else ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE)
  p$plot <- p$plot + ggplot2::guides(colour = guide, fill = guide, linetype = guide) +
    ggplot2::theme(
      legend.position = "top", legend.justification = "center", legend.box = "horizontal",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_font_size),
      legend.key.width = grid::unit(0.8, "lines"),
      legend.key.height = grid::unit(0.8, "lines"),
      legend.spacing.x = grid::unit(3, "pt"),
      legend.margin = ggplot2::margin(3, 3, 3, 3),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      plot.margin = ggplot2::margin(8, 12, 5, 8)
    )
  lookup <- stats::setNames(risk_labels, legend_labs)
  # survminer uses y = rev(strata) plus reversed labels internally. Override
  # BOTH together so shortening labels never associates a count with the wrong group.
  strata <- NULL # avoid a package-check note; aes() resolves the data column.
  p$table <- suppressMessages(p$table + ggplot2::aes(y = strata) +
    ggplot2::scale_y_discrete(limits = rev(legend_labs), breaks = rev(legend_labs),
                            labels = unname(lookup[rev(legend_labs)]))) +
    ggplot2::labs(y = NULL) + ggplot2::theme(
    legend.position = "none", axis.text.y = ggplot2::element_text(size = 9, colour = "black"),
    plot.title = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(4, 12, 5, 8)
  )
  p
}

make_survival_annotation <- function(cox_df, test_var, group_key, overall_p,
                                     conf_int = 0.95, component_interaction_p = "N/A",
                                     font_size = 9) {
  rows <- cox_df[cox_df$Term == test_var, , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  hr_col <- grep("^a?HR_.*CI$", names(rows), value = TRUE)[1L]
  if (is.na(hr_col)) return(NULL)
  idx <- match(group_key$raw, rows$Group)
  keep <- !is.na(idx)
  if (!any(keep)) return(NULL)
  tab <- data.frame(Group = group_key$id[keep], Estimate = rows[[hr_col]][idx[keep]],
                    P = rows$P_Value[idx[keep]], stringsAsFactors = FALSE)
  names(tab) <- c("Group", paste0(if (startsWith(hr_col, "aHR")) "aHR" else "HR",
                                " (", format(100 * conf_int, trim = TRUE), "% CI)"), "p-value")
  foot <- if (!identical(overall_p, "N/A")) paste0("Log-rank p ", if (startsWith(overall_p, "<")) "" else "= ", overall_p) else
    "Log-rank: not applicable"
  ip <- unique(rows$P_for_Interaction[!is.na(rows$P_for_Interaction) & rows$P_for_Interaction != "N/A"])
  if (length(ip)) foot <- c(foot, paste0("External interaction p: ", paste(ip, collapse = " / ")))
  if (!identical(component_interaction_p, "N/A")) foot <- c(foot, paste0("Component interaction p: ", component_interaction_p))
  foot <- survival_wrap(foot, 52)
  font_size <- max(6, min(18, as.numeric(font_size)))
  # Measure on a real raster device, not the user's current device.
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 1600, height = 1600, res = 160)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  theme <- gridExtra::ttheme_minimal(
    base_size = font_size, padding = grid::unit(c(2.1, 2.8), "mm"),
    core = list(fg_params = list(hjust = 0, x = 0.05), bg_params = list(fill = "white", col = NA)),
    colhead = list(fg_params = list(hjust = 0, x = 0.05, fontface = "bold"),
                   bg_params = list(fill = "#F0F2F4", col = NA))
  )
  g <- gridExtra::tableGrob(tab, rows = NULL, theme = theme)
  footer <- grid::textGrob(paste(foot, collapse = "\n"), x = grid::unit(5, "pt"), hjust = 0,
                           gp = grid::gpar(fontsize = font_size - 0.5, lineheight = 1.15))
  footer_width <- grid::convertWidth(grid::grobWidth(footer) + grid::unit(10, "pt"), "inches", valueOnly = TRUE)
  table_width <- grid::convertWidth(sum(g$widths), "inches", valueOnly = TRUE)
  if (footer_width > table_width) g$widths[length(g$widths)] <- g$widths[length(g$widths)] +
    grid::unit(footer_width - table_width, "inches")
  g <- gtable::gtable_add_rows(g, grid::grobHeight(footer) + grid::unit(10, "pt"))
  g <- gtable::gtable_add_grob(g, footer, t = nrow(g), l = 1, r = ncol(g), clip = "on")
  g <- gtable::gtable_add_grob(g, grid::rectGrob(gp = grid::gpar(fill = "white", col = "#737B84", lwd = 0.6)),
                              t = 1, l = 1, b = nrow(g), r = ncol(g), z = -Inf)
  list(grob = g,
       width = grid::convertWidth(sum(g$widths), "inches", valueOnly = TRUE),
       height = grid::convertHeight(sum(g$heights), "inches", valueOnly = TRUE),
       alt = paste(c(apply(tab, 1, paste, collapse = ": "), foot), collapse = "; "))
}

normalise_annotation_position <- function(position, annotation, width, height) {
  if (!is.numeric(position) || length(position) != 2L || any(!is.finite(position))) {
    stop("annotation_position must be two finite numbers: x and y from the top-left, in [0, 1].")
  }
  if (is.null(annotation)) return(stats::setNames(pmin(1, pmax(0, position)), c("x", "y")))
  if (annotation$width > width || annotation$height > height) {
    stop("The statistics box is larger than the figure. Increase figure dimensions or reduce annotation font size.")
  }
  stats::setNames(pmax(0, pmin(as.numeric(position), c(1 - annotation$width / width, 1 - annotation$height / height))), c("x", "y"))
}

draw_survival_annotation <- function(annotation, position, width, height) {
  if (is.null(annotation)) return(invisible(NULL))
  position <- normalise_annotation_position(position, annotation, width, height)
  grid::upViewport(0)
  grid::pushViewport(grid::viewport(
    x = grid::unit(position[1L], "npc"), y = grid::unit(1 - position[2L], "npc"),
    width = grid::unit(annotation$width, "inches"), height = grid::unit(annotation$height, "inches"),
    just = c("left", "top"), clip = "off"
  ))
  on.exit(grid::popViewport(), add = TRUE)
  grid::grid.draw(annotation$grob)
  invisible(position)
}

save_survival_annotation <- function(annotation, filename, res = 180) {
  stopifnot(!is.null(annotation))
  grDevices::png(filename, width = ceiling(annotation$width * res),
                 height = ceiling(annotation$height * res), res = res, bg = "transparent")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(annotation$grob)
  invisible(filename)
}

#' Readable result-column names and definitions
#' @param columns Character vector of raw result column names.
#' @param conf_int Confidence level used for the analysis.
#' @export
survival_column_dictionary <- function(columns, conf_int = 0.95) {
  labels <- c(Scenario = "Analysis", Outcome = "Outcome", Variable_Term = "Variable / term",
              Group = "Group", N = "N", N_events = "Events",
              Median_Survival = paste0("Median survival (", format(100 * conf_int, trim = TRUE), "% CI)"),
              P_Value = "Cox p-value", Overall_LogRank_P = "Overall log-rank p-value",
              Variable_PH_p_value = "Variable PH p-value", Global_PH_p_value = "Global PH p-value",
              P_for_Interaction = "External interaction p-value",
              P_for_TestVar_Interaction = "Group-component interaction p-value",
              Adjusted_For = "Adjusted for")
  descriptions <- c(
    Scenario = "Name assigned to this analysis.",
    Outcome = "Outcome being analysed; the event indicator determines which observations are events.",
    Variable_Term = "Predictor associated with this row of the Cox model. Includes the primary grouping variable and, when present, adjustment covariates.",
    Group = "Category or combined category represented by this row. G1, G2, etc. link primary groups to the plot and group key. A continuous predictor is interpreted per one-unit increase.",
    N = "Number of analysed participants in this primary group, after complete-case exclusions. For a continuous primary predictor, total complete-case model N is repeated on each row. Categorical-model covariate-only rows have no separate group N.",
    N_events = "Number of observed events in this analysed primary group. Censored observations are not events. For a continuous primary predictor, the total model event count is repeated. Covariate-only rows of a categorical model have no separate event count.",
    Median_Survival = "Kaplan-Meier median survival and its confidence interval, in months (the app assumes time is in months). NR means the median was not reached. Not applicable to a continuous primary predictor or covariate-only rows. An unavailable interval bound may be NA.",
    P_Value = "Two-sided Wald p-value for this Cox coefficient (null: HR = 1). Ref marks the reference category. This is not the overall log-rank p-value and is not an omnibus test of a multi-level variable.",
    Overall_LogRank_P = "Unadjusted overall log-rank test comparing Kaplan-Meier curves across the primary groups; repeated across rows. Not applicable with one group or a continuous primary predictor.",
    Variable_PH_p_value = "Term-level proportional-hazards diagnostic from cox.zph, based on scaled Schoenfeld residuals. Small p-values suggest departure from proportional hazards; they do not test whether HR = 1. Repeated across levels of the same term.",
    Global_PH_p_value = "Global proportional-hazards diagnostic across the Cox model. Repeated across rows; a small p-value suggests a possible violation in at least part of the model.",
    P_for_Interaction = "Likelihood-ratio test comparing the main-effects Cox model with a model that adds an interaction between the primary test variable and a selected external covariate. Multiple values separated by / follow the selected interaction-variable order. This is distinct from interaction among combined grouping components.",
    P_for_TestVar_Interaction = "Likelihood-ratio test comparing additive and fully crossed models of the variables used to create a combined group. With more than two components this jointly tests the added interaction terms. N/A for a single primary test variable.",
    Adjusted_For = "Other covariates included in the main Cox model. Unadjusted (univariate) means no additional adjustment covariates were included. Kaplan-Meier curves and log-rank tests remain unadjusted."
  )
  out <- data.frame(key = columns, label = unname(labels[columns]),
                    description = unname(descriptions[columns]), stringsAsFactors = FALSE)
  for (i in seq_along(columns)) {
    nm <- columns[i]
    if (grepl("^a?HR_.*CI$", nm)) {
      # Use the actual confidence level; legacy raw column names round it.
      ci <- format(100 * conf_int, trim = TRUE)
      adjusted <- startsWith(nm, "aHR")
      out$label[i] <- paste0(if (adjusted) "Adjusted HR" else "HR", " (", ci, "% CI)")
      out$description[i] <- paste0(if (adjusted) "Adjusted hazard ratio" else "Hazard ratio",
        " and ", ci, "% confidence interval from the Cox model. For categorical terms, relative to their model reference category; for a continuous term, per one-unit increase. 1.00 (Ref) marks the primary reference group. A hazard ratio is not a probability or risk ratio.")
    }
    if (is.na(out$label[i])) out$label[i] <- survival_display_name(nm)
    if (is.na(out$description[i])) out$description[i] <- paste("Result column:", out$label[i])
  }
  out
}

#' Create a presentation copy without changing the machine-readable results
#' @param results Raw results from MARData_surv_analysis().
#' @param conf_int Confidence level used for the analysis.
#' @param group_key Optional display key from return_details = TRUE.
#' @param main_term Optional primary model term, used to avoid relabelling covariate rows.
#' @param variable_labels Optional named vector mapping internal variable names to display names.
#' @export
format_survival_results <- function(results, conf_int = 0.95, group_key = NULL,
                                    main_term = NULL, variable_labels = NULL) {
  x <- results
  if (!is.null(group_key) && !is.null(main_term)) {
    idx <- match(as.character(x$Group), group_key$raw)
    use <- !is.na(idx) & !is.na(x$Variable_Term) & x$Variable_Term == main_term
    x$Group[use] <- paste0(group_key$id[idx[use]], ": ", group_key$label[idx[use]])
  }
  friendly <- function(v) {
    if (!is.null(variable_labels) && v %in% names(variable_labels)) variable_labels[[v]]
    else survival_display_name(v)
  }
  if ("Variable_Term" %in% names(x)) x$Variable_Term <- vapply(as.character(x$Variable_Term), friendly, character(1))
  if ("Adjusted_For" %in% names(x)) x$Adjusted_For <- vapply(as.character(x$Adjusted_For), function(s) {
    if (is.na(s) || s == "Unadjusted (univariate)") return(s)
    paste(vapply(strsplit(s, ", ", fixed = TRUE)[[1L]], friendly, character(1)), collapse = ", ")
  }, character(1))
  names(x) <- survival_column_dictionary(names(results), conf_int)$label
  x
}

# Rebuild presentation only: no survfit(), coxph(), survdiff() or interaction tests.
update_survival_display <- function(details, group_labels = NULL, legend_font_size = NULL,
                                     legend_ncol = NULL, legend_wrap_width = NULL,
                                     risk_table_labels = NULL, annotation_font_size = NULL) {
  if (is.null(details$plot_args)) return(details)
  if (!is.null(group_labels)) {
    if (length(group_labels) != nrow(details$group_key) || anyNA(group_labels) ||
        any(!nzchar(trimws(group_labels)))) stop("Supply a non-empty label for every group.")
    details$group_key$label <- as.character(group_labels)
  }
  args <- details$plot_args
  args$group_key <- details$group_key
  for (nm in c("legend_font_size", "legend_ncol", "legend_wrap_width", "risk_table_labels")) {
    value <- get(nm)
    if (!is.null(value)) args[[nm]] <- value
  }
  details$plot <- do.call(build_survival_plot, args)
  details$plot_args <- args
  if (!is.null(annotation_font_size)) details$annotation_args$font_size <- annotation_font_size
  details$annotation_args$group_key <- details$group_key
  details$annotation <- do.call(make_survival_annotation, details$annotation_args)
  details
}
