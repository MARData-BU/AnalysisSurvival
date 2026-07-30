get_km_medians <- function(surv_formula, data, conf_level = 0.95) {
  fit_km <- survfit(surv_formula, data = data, conf.type = "log-log", conf.int = conf_level)
  fit_km$call$formula <- surv_formula
  
  km_tab <- summary(fit_km)$table
  if (is.null(dim(km_tab))) {
    km_tab <- matrix(km_tab, nrow = 1, dimnames = list("Overall", names(km_tab)))
  }

  group_names <- gsub("^[^=]*=", "", rownames(km_tab))
  strata_counts <- if (!is.null(fit_km$strata)) fit_km$n else nrow(data)
  
  lcl_col <- grep("LCL$", colnames(km_tab), value = TRUE)
  ucl_col <- grep("UCL$", colnames(km_tab), value = TRUE)
  
  formatted_medians <- apply(km_tab, 1, function(row) {
    if (is.na(row["median"])) return("NR")
    sprintf("%.2f (%.2f-%.2f)", row["median"], row[lcl_col], row[ucl_col])
  })
  
  ci_label <- paste0("Median_Survival_", round(conf_level * 100), "CI")
  
  res <- data.frame(
    Group = group_names,
    N = as.vector(strata_counts),
    Median_CI = formatted_medians,
    stringsAsFactors = FALSE
  )
  colnames(res)[3] <- ci_label
  return(res)
}

get_cox_hrs <- function(data, time_var, event_var, test_var,
                        covariates = c(), interaction_var = NULL,
                        ref_level = NULL, conf_level = 0.95) {
  
  # A numeric test_var is treated as continuous 
  test_var_is_numeric <- is.numeric(data[[test_var]])
  
  if (test_var_is_numeric) {
    if (!is.null(ref_level)) {
      warning("ref_level is ignored because '", test_var, "' is numeric (continuous); ",
              "reference levels only apply to categorical variables.")
    }
  } else {
    data <- set_ref_level(data, test_var, ref_level)
  }
  
  is_adjusted <- length(covariates) > 0
  ci_label <- paste0(if (is_adjusted) "aHR_" else "HR_", round(conf_level * 100), "CI")
  
  if (!test_var_is_numeric) {
    km_formula_test <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, event_var, test_var))
    fit_km <- survfit(km_formula_test, data = data)
    n_strata <- length(fit_km$strata)
    
    if (n_strata <= 1) {
      res <- data.frame(
        Term = test_var, Group = "Overall", HR_CI = "1.00 (Ref)",
        P_Value = "N/A", Variable_PH_p_value = "N/A", Global_PH_p_value = "N/A",
        P_for_Interaction = "N/A", stringsAsFactors = FALSE
      )
      colnames(res)[3] <- ci_label
      return(res)
    }
  }
  
  # Full model: test_var + covariates
  all_vars <- unique(c(test_var, covariates))
  main_formula <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, event_var, paste(all_vars, collapse = " + ")))
  cox_main <- tryCatch(coxph(main_formula, data = data), error = function(e) NULL)
  if (is.null(cox_main)) {
    warning("Cox model failed to converge.")
    return(NULL)
  }
  
  ph_test <- tryCatch(cox.zph(cox_main), error = function(e) NULL)
  global_ph_p <- "N/A"
  ph_by_term  <- character(0)
  
  if (!is.null(ph_test)) {
    ph_tab <- ph_test$table
    if (any(ph_tab[, "p"] < 0.05, na.rm = TRUE)) {
      warning("The Proportional Hazards assumption is not met for at least one term. Please review cox.zph() output.")
    }
    if ("GLOBAL" %in% rownames(ph_tab)) {
      g_p <- ph_tab["GLOBAL", "p"]
      global_ph_p <- ifelse(g_p < 0.001, "<0.001", sprintf("%.3f", g_p))
    }
    term_rows <- setdiff(rownames(ph_tab), "GLOBAL")
    ph_by_term <- ifelse(ph_tab[term_rows, "p"] < 0.001, "<0.001", sprintf("%.3f", ph_tab[term_rows, "p"]))
    names(ph_by_term) <- term_rows
  }
  
  cox_sum    <- summary(cox_main)
  coef_mat   <- cox_sum$coefficients
  term_names <- rownames(coef_mat)
  hrs        <- as.vector(coef_mat[, "exp(coef)"])
  p_vals     <- as.vector(coef_mat[, "Pr(>|z|)"])
  cis        <- cox_sum$conf.int[, c(3, 4), drop = FALSE]
  
  formatted_hrs <- sprintf("%.2f (%.2f-%.2f)", hrs, cis[, 1], cis[, 2])
  formatted_p   <- ifelse(p_vals < 0.001, "<0.001", sprintf("%.3f", p_vals))

  vars_by_length <- all_vars[order(-nchar(all_vars))]
  owner_var <- vapply(term_names, function(tn) {
    if (tn %in% all_vars) return(tn)
    hit <- vars_by_length[startsWith(tn, vars_by_length)]
    if (length(hit) > 0) hit[1] else tn
  }, character(1))
  
  var_ph_p_col <- sapply(owner_var, function(v) {
    if (v %in% names(ph_by_term)) ph_by_term[[v]] else "N/A"
  })
  
  # --- Optional single interaction term (test_var * interaction_var) ---
  interaction_col <- rep("N/A", length(term_names))
  if (!is.null(interaction_var) && length(interaction_var) > 0) {
    valid_int_vars <- intersect(interaction_var, all_vars)
    test_var_pattern <- paste0("^", test_var)
    
    for (int_v in valid_int_vars) {
      other_covars <- setdiff(covariates, int_v)
      int_term <- paste0(test_var, " * ", int_v)
      rhs <- paste(c(int_term, other_covars), collapse = " + ")
      int_formula <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, event_var, rhs))
      cox_int <- tryCatch(coxph(int_formula, data = data), error = function(e) NULL)
      
      if (!is.null(cox_int)) {
        lrt <- anova(cox_main, cox_int, test = "LRT")
        p_int <- lrt$`Pr(>|Chi|)`[2]
        if (!is.na(p_int)) {
          p_int_str <- if (p_int < 0.001) "<0.001" else sprintf("%.4f", p_int)
          matching_terms <- grepl(paste0("^", int_v), term_names) | grepl(test_var_pattern, term_names)
          for (idx in which(matching_terms)) {
            interaction_col[idx] <- if (interaction_col[idx] == "N/A") p_int_str else paste(interaction_col[idx], p_int_str, sep = " / ")
          }
        }
      }
    }
  }
  
  is_test_var_term <- owner_var == test_var
  clean_groups <- ifelse(term_names == owner_var, owner_var, substring(term_names, nchar(owner_var) + 1))
  
  res_non_ref <- data.frame(
    Term = owner_var, Group = clean_groups, HR_CI = formatted_hrs,
    P_Value = formatted_p, Variable_PH_p_value = var_ph_p_col,
    Global_PH_p_value = global_ph_p, P_for_Interaction = interaction_col,
    stringsAsFactors = FALSE
  )
  
  if (test_var_is_numeric) {
    res <- res_non_ref
    colnames(res)[3] <- ci_label
    return(res)
  }
  
  # --- Explicit reference row for test_var (categorical only) ---
  all_levels <- levels(as.factor(data[[test_var]]))
  ref_level_label <- setdiff(all_levels, clean_groups[is_test_var_term])[1]
  ref_ph_p  <- if (any(is_test_var_term)) var_ph_p_col[which(is_test_var_term)[1]] else "N/A"
  ref_int_p <- if (any(is_test_var_term)) interaction_col[which(is_test_var_term)[1]] else "N/A"
  
  res_ref <- data.frame(
    Term = test_var, Group = ref_level_label, HR_CI = "1.00 (Ref)", P_Value = "Ref",
    Variable_PH_p_value = ref_ph_p, Global_PH_p_value = global_ph_p,
    P_for_Interaction = ref_int_p, stringsAsFactors = FALSE
  )
  
  res <- rbind(res_ref, res_non_ref)
  colnames(res)[3] <- ci_label
  return(res)
}

generate_surv_annotation <- function(km_formula, data, cox_df, test_var, conf_level = 0.95) {
  fit_km <- survfit(km_formula, data = data)
  n_strata <- length(fit_km$strata)
  if (n_strata == 0) return(NULL)
  
  surv_diff <- survdiff(km_formula, data = data)
  p_val <- pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1, lower.tail = FALSE)
  p_str <- if (p_val < 0.001) "Log-rank p < 0.001" else sprintf("Log-rank p = %.4f", p_val)
  
  hr_col <- grep("^a?HR_.*CI$", colnames(cox_df), value = TRUE)[1]
  test_rows <- cox_df[cox_df$Term == test_var, ]
  if (nrow(test_rows) == 0) return(NULL)
  
  lines <- sapply(seq_len(nrow(test_rows)), function(i) {
    if (test_rows[[hr_col]][i] == "1.00 (Ref)") {
      paste0(test_rows$Group[i], ": Ref")
    } else {
      sprintf("%s: HR = %s (p = %s)", test_rows$Group[i], test_rows[[hr_col]][i], test_rows$P_Value[i])
    }
  })
  
  lines <- c(lines, p_str)
  
  p_int_val <- unique(test_rows$P_for_Interaction[test_rows$P_for_Interaction != "N/A"])
  if (length(p_int_val) > 0) {
    lines <- c(lines, paste0("P for Int = ", paste(p_int_val, collapse = " / ")))
  }
  
  paste(lines, collapse = "\n")
}
