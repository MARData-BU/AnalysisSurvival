# HELPER FUNCIONS AND VARIABLES

default_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#E6AB02", "#66A61E")

format_level_list <- function(levels) {
  levels <- as.character(levels)
  if (length(levels) <= 1) return(levels)
  paste0(paste(levels[-length(levels)], collapse = ", "), " and ", levels[length(levels)])
}

set_ref_level <- function(data, test_var, ref_level = NULL) {
  data[[test_var]] <- as.factor(data[[test_var]])
  if (!is.null(ref_level)) {
    lv <- levels(data[[test_var]])
    if (!(ref_level %in% lv)) {
      stop(sprintf("Level %s is not in variable %s; the available levels are %s",
                   ref_level, test_var, format_level_list(lv)), call. = FALSE)
    }
    data[[test_var]] <- stats::relevel(data[[test_var]], ref = ref_level)
  }
  data
}

choose_break_time <- function(max_time, n_breaks = 6) {
  if (!is.finite(max_time) || max_time <= 0) return(5)
  step <- diff(pretty(c(0, max_time), n = n_breaks))[1]
  if (!is.finite(step) || step <= 0) return(5)
  step
}

save_survival_plot <- function(p, filename, width = 8, height = 6, res = 300,
                               annotation = NULL,
                               annotation_position = c(x = 0.14, y = 0.42)) {
  ext <- tolower(tools::file_ext(filename))
  raster_devices <- list(png = grDevices::png, jpeg = grDevices::jpeg, jpg = grDevices::jpeg, tiff = grDevices::tiff, bmp = grDevices::bmp)
  if (!(ext %in% c(names(raster_devices), "pdf", "svg"))) {
    filename <- paste0(filename, ".png")
    ext <- "png"
  }
 
  if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = height)
  } else if (ext == "svg") {
    grDevices::svg(filename, width = width, height = height)
  } else {
    raster_devices[[ext]](filename, width = width, height = height, units = "in", res = res)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  if (!is.null(annotation)) draw_survival_annotation(annotation, annotation_position, width, height)
 
  invisible(filename)
}

save_summary_table <- function(df, filename) {
  ext <- tolower(tools::file_ext(filename))
  if (ext != "csv") filename <- paste0(filename, ".csv")
  utils::write.table(df, file = filename, sep = ",", row.names = FALSE, col.names = TRUE)
  invisible(filename)
}
