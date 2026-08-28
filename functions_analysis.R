# =============================================================================
# functions_analysis.R
# Statistical models: LMER, ART, screening, significance tables
# =============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(emmeans)
  library(ARTool)
})

recode_hormone_labels <- function(df, exp_short, e2_label, hormone_levels) {
  df %>%
    dplyr::mutate(
      HORMONE = dplyr::case_match(
        as.character(HORMONE),
        "NONE" ~ exp_short,
        "Estradiol" ~ e2_label,
        .default = as.character(HORMONE)
      ),
      HORMONE = factor(HORMONE, levels = hormone_levels)
    )
}

make_log2fc_long <- function(df, cytokine_cols,
                             pbs_level         = PBS_LEVEL,
                             exclude_exposures = EXCLUDE_EXPOSURES,
                             pseudocount       = PSEUDOCOUNT) {
  d <- df %>%
    dplyr::filter(!EXPOSURE %in% exclude_exposures) %>%
    dplyr::mutate(
      SEX         = as.character(SEX),
      CELLTYPE    = factor(CELLTYPE),
      HORMONE     = factor(HORMONE),
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  if ("TIMEPOINT" %in% names(d)) {
    d <- d %>% dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }

  cytokine_cols <- intersect(cytokine_cols, names(d))
  if (length(cytokine_cols) == 0) {
    return(tibble::tibble(PATIENTCODE = character(), SEX = character(), CELLTYPE = factor(), HORMONE = factor(), TIMEPOINT = factor(), EXPOSURE = factor(), CYTOKINE = character(), VALUE = numeric(), log2_val = numeric(), log2_pbs = numeric(), log2FC = numeric()))
  }

  long <- d %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cytokine_cols), names_to = "CYTOKINE", values_to = "VALUE") %>%
    dplyr::filter(!is.na(VALUE)) %>%
    dplyr::mutate(log2_val = log2(as.numeric(VALUE) + pseudocount))

  join_keys <- c("PATIENTCODE", "SEX", "CELLTYPE", "HORMONE", "CYTOKINE")
  if ("TIMEPOINT" %in% names(long)) join_keys <- c(join_keys, "TIMEPOINT")

  pbs <- long %>%
    dplyr::filter(EXPOSURE == pbs_level) %>%
    dplyr::select(PATIENTCODE, SEX, CELLTYPE, HORMONE, dplyr::any_of("TIMEPOINT"), CYTOKINE, log2_pbs = log2_val)

  long %>%
    dplyr::filter(EXPOSURE != pbs_level) %>%
    dplyr::left_join(pbs, by = join_keys) %>%
    dplyr::mutate(log2FC = log2_val - log2_pbs)
}

exposure_lmer_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                   response_columns = NULL,
                                   ctrl_level = "PBS_Control") {
  if (group != "All" && "SEX" %in% names(data)) data <- data %>% dplyr::filter(SEX == group)
  if ("TIMEPOINT" %in% names(data)) {
    data <- data %>% dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }

  stopifnot("EXPOSURE" %in% names(data), "PATIENTCODE" %in% names(data))
  data <- data %>% dplyr::mutate(EXPOSURE = factor(EXPOSURE), PATIENTCODE = factor(PATIENTCODE))
  if (ctrl_level %in% levels(data$EXPOSURE)) data$EXPOSURE <- relevel(data$EXPOSURE, ref = ctrl_level)
  if (is.null(response_columns)) response_columns <- names(data)[sapply(data, is.numeric)]

  results_list <- list()
  for (resp in response_columns) {
    if (!resp %in% names(data)) next
    df <- data %>% dplyr::filter(!is.na(.data[[resp]]))
    if (length(unique(df$EXPOSURE)) < 2) next
    tp_term <- if ("TIMEPOINT" %in% names(df) && dplyr::n_distinct(df$TIMEPOINT[!is.na(df$TIMEPOINT)]) > 1) {
      " + TIMEPOINT"
    } else ""
    model <- try(lme4::lmer(as.formula(paste0(resp, " ~ EXPOSURE", tp_term, " + (1|PATIENTCODE)")), data = df), silent = TRUE)
    if (inherits(model, "try-error")) next
    emm      <- emmeans::emmeans(model, ~ EXPOSURE, weights = "equal")
    ctrl_idx <- which(levels(data$EXPOSURE) == ctrl_level)
    pairwise <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ctrl_idx, adjust = adjust_method)
    pairwise_df          <- as.data.frame(summary(pairwise))
    pairwise_df$response <- resp
    results_list[[resp]] <- pairwise_df
  }
  dplyr::bind_rows(results_list)
}

interaction_lmer_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                      response_columns = NULL,
                                      ctrl_level = "PBS_Control") {
  if (group != "All" && "SEX" %in% names(data)) data <- data %>% dplyr::filter(SEX == group)
  if ("TIMEPOINT" %in% names(data)) {
    data <- data %>% dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }

  stopifnot("EXPOSURE" %in% names(data), "PATIENTCODE" %in% names(data), "SEX" %in% names(data))
  data <- data %>% dplyr::mutate(EXPOSURE = factor(EXPOSURE), PATIENTCODE = factor(PATIENTCODE), SEX = factor(SEX))
  if (ctrl_level %in% levels(data$EXPOSURE)) data$EXPOSURE <- relevel(data$EXPOSURE, ref = ctrl_level)
  if (is.null(response_columns)) response_columns <- names(data)[sapply(data, is.numeric)]

  results_list <- list()
  for (resp in response_columns) {
    if (!resp %in% names(data)) next
    df <- data %>% dplyr::filter(!is.na(.data[[resp]]))
    if (length(unique(df$EXPOSURE)) < 2 || length(unique(df$SEX)) < 2) next
    tp_term <- if ("TIMEPOINT" %in% names(df) && dplyr::n_distinct(df$TIMEPOINT[!is.na(df$TIMEPOINT)]) > 1) {
      " + TIMEPOINT"
    } else ""
    model <- try(lme4::lmer(as.formula(paste0(resp, " ~ EXPOSURE * SEX", tp_term, " + (1|PATIENTCODE)")), data = df), silent = TRUE)
    if (inherits(model, "try-error")) next
    emm_exp  <- emmeans::emmeans(model, ~ EXPOSURE, weights = "equal")
    ctrl_idx <- which(levels(data$EXPOSURE) == ctrl_level)
    pw_exp   <- emmeans::contrast(emm_exp, "trt.vs.ctrl", ref = ctrl_idx, adjust = adjust_method)
    exp_df   <- as.data.frame(summary(pw_exp)); exp_df$type <- "Exposure_vs_Control"; exp_df$response <- resp
    emm_int  <- emmeans::emmeans(model, ~ SEX | EXPOSURE)
    pw_int   <- emmeans::contrast(emm_int, "pairwise", simple = "SEX", combine = TRUE)
    int_df   <- as.data.frame(summary(pw_int)); int_df$type <- "Sex_within_Exposure"; int_df$response <- resp
    results_list[[resp]] <- dplyr::bind_rows(exp_df, int_df)
  }
  dplyr::bind_rows(results_list)
}

exposure_art_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                  response_columns = NULL,
                                  ctrl_level = "PBS_Control") {
  filtered_data <- switch(group,
                          "M"   = data %>% dplyr::filter(EXPOSURE != "Untreated_Control", SEX == "M"),
                          "F"   = data %>% dplyr::filter(EXPOSURE != "Untreated_Control", SEX == "F"),
                          "All" = data %>% dplyr::filter(EXPOSURE != "Untreated_Control"),
                          stop("Invalid group. Choose 'All', 'M', or 'F'.")) %>%
    dplyr::mutate(SEX = factor(SEX), EXPOSURE = factor(EXPOSURE), PATIENTCODE = factor(PATIENTCODE))
  if ("TIMEPOINT" %in% names(filtered_data)) {
    filtered_data <- filtered_data %>% dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }
  if (is.null(response_columns)) response_columns <- names(filtered_data)[sapply(filtered_data, is.numeric)]

  results_list <- list()
  for (response in response_columns) {
    if (!response %in% names(filtered_data)) next
    tp_term <- if ("TIMEPOINT" %in% names(filtered_data) && dplyr::n_distinct(filtered_data$TIMEPOINT[!is.na(filtered_data$TIMEPOINT)]) > 1) {
      " + TIMEPOINT"
    } else ""
    formula <- if (group == "All") {
      as.formula(paste0(response, " ~ EXPOSURE * SEX", tp_term, " + (1|PATIENTCODE)"))
    } else {
      as.formula(paste0(response, " ~ EXPOSURE", tp_term, " + (1|PATIENTCODE)"))
    }
    m.art <- try(ARTool::art(formula, data = filtered_data), silent = TRUE)
    if (inherits(m.art, "try-error")) next
    anova_res   <- anova(m.art)
    exposure_p  <- { r <- anova_res[grepl("^EXPOSURE$", anova_res[[1]], ignore.case = TRUE), ]; if (nrow(r) > 0) r[["Pr(>F)"]][1] else NA }
    interaction_p <- if (group == "All") {
      r <- anova_res[grepl("EXPOSURE:SEX", anova_res[[1]], ignore.case = TRUE), ]; if (nrow(r) > 0) r[["Pr(>F)"]][1] else NA
    } else NA
    pairwise_con <- try(ARTool::art.con(m.art, "EXPOSURE", adjust = adjust_method), silent = TRUE)
    if (inherits(pairwise_con, "try-error")) next
    ctrl_idx <- which(levels(filtered_data$EXPOSURE) == ctrl_level)
    pairwise <- emmeans::contrast(pairwise_con, "trt.vs.ctrl", ref = ctrl_idx, adjust = adjust_method)
    res <- as.data.frame(pairwise); res$response <- response; res$exposure_p <- exposure_p; res$interaction_p <- interaction_p
    results_list[[response]] <- res
  }
  dplyr::bind_rows(results_list)
}

run_lmer_chunk <- function(label, celltype_filter, hormone_filter,
                           sala_full = NULL,
                           llod_table = cytokine_llod,
                           zero_co = ZERO_CUTOFF,
                           alpha_q = ALPHA_Q,
                           trend_a = TREND_ALPHA) {
  cat("\n── LMER:", label, "──\n")
  if (is.null(sala_full)) stop("sala_full must be provided as a data frame containing the full SALA dataset")
  d_sub <- sala_full %>% dplyr::filter(CELLTYPE == celltype_filter, HORMONE == hormone_filter)
  if (nrow(d_sub) == 0) {
    warning(sprintf("[%s] No rows for CELLTYPE=%s HORMONE=%s; skipping.", label, celltype_filter, hormone_filter))
    return(invisible(list(sig_table = tibble::tibble(), lmer_interaction = tibble::tibble(), lmer_plot_data = tibble::tibble())))
  }
  imp_res <- impute_lod_sqrt2(input_data = d_sub, cols = intersect(names(d_sub), llod_table$Analyte), cytokine_llod = llod_table, zero_cutoff = zero_co)
  d_imp <- imp_res$data
  valid_cyts <- intersect(imp_res$valid_cytokines, names(d_imp))
  if (length(valid_cyts) == 0) {
    warning(sprintf("[%s] No valid cytokines after imputation/filters; skipping.", label))
    return(invisible(list(sig_table = tibble::tibble(), lmer_interaction = tibble::tibble(), lmer_plot_data = tibble::tibble())))
  }
  d_filt <- d_imp %>% dplyr::filter(EXPOSURE != "Untreated_Control") %>% dplyr::mutate(EXPOSURE = factor(EXPOSURE), SEX = factor(SEX), PATIENTCODE = factor(PATIENTCODE))
  if ("TIMEPOINT" %in% names(d_filt)) d_filt <- d_filt %>% dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  if (!(PBS_LEVEL %in% unique(as.character(d_filt$EXPOSURE)))) {
    warning(sprintf("[%s] PBS level '%s' absent in this stratum; skipping.", label, PBS_LEVEL))
    return(invisible(list(sig_table = tibble::tibble(), lmer_interaction = tibble::tibble(), lmer_plot_data = tibble::tibble())))
  }
  msd_sum <- summarize_to_wide(d_filt, measure_vars = valid_cyts)
  if (!("PBS_Control" %in% names(msd_sum))) {
    warning(sprintf("[%s] summarize_to_wide() missing PBS_Control; skipping.", label))
    return(invisible(list(sig_table = tibble::tibble(), lmer_interaction = tibble::tibble(), lmer_plot_data = tibble::tibble())))
  }
  pbs_ctl <- msd_sum %>% dplyr::select(Measurement, `PBS_Control`) %>% dplyr::arrange(Measurement)
  lmer_All <- exposure_lmer_pairwise(d_filt, "All", "fdr", valid_cyts)
  lmer_F   <- exposure_lmer_pairwise(d_filt, "F",   "fdr", valid_cyts)
  lmer_M   <- exposure_lmer_pairwise(d_filt, "M",   "fdr", valid_cyts)
  lmer_int <- interaction_lmer_pairwise(d_filt, "All", "fdr", valid_cyts)
  if (nrow(lmer_All) == 0 && nrow(lmer_F) == 0 && nrow(lmer_M) == 0) {
    warning(sprintf("[%s] No successful LMER results; skipping table export.", label))
    return(invisible(list(sig_table = tibble::tibble(), lmer_interaction = lmer_int, lmer_plot_data = tibble::tibble())))
  }
  lmer_plot_data <- lmer_results_to_plot_format(lmer_All, lmer_F, lmer_M, celltype = celltype_filter, hormone = hormone_filter, alpha = alpha_q)
  out_csv <- paste0("MSD_SALA_", label, "_lmer.csv")
  sig_tbl <- build_lmer_sig_table(lmer_All = lmer_All, lmer_F = lmer_F, lmer_M = lmer_M, msd_summary = msd_sum, pbs_control = pbs_ctl, cytokine_llod = llod_table, out_filename = out_csv, alpha = alpha_q, trend_alpha = trend_a)
  cat("  Saved:", out_csv, "\n")
  invisible(list(sig_table = sig_tbl, lmer_interaction = lmer_int, lmer_plot_data = lmer_plot_data))
}

build_lmer_sig_table <- function(lmer_All, lmer_F, lmer_M,
                                 msd_summary, pbs_control,
                                 cytokine_llod,
                                 out_filename,
                                 alpha = 0.1,
                                 trend_alpha = 0.2,
                                 epsilon = 1e-5) {
  fmt <- function(df, grp) {
    if (is.null(df) || nrow(df) == 0) return(tibble::tibble())
    if (!all(c("response", "contrast", "p.value") %in% names(df))) return(tibble::tibble())
    col <- paste0("response_", grp)
    df %>% dplyr::mutate(!!col := paste0(response, "_", grp), EXPOSURE = stringr::str_replace(contrast, " - PBS_Control", ""))
  }
  combined <- dplyr::bind_rows(fmt(lmer_All, "All"), fmt(lmer_F, "F"), fmt(lmer_M, "M"))
  if (nrow(combined) == 0) { save_table(tibble::tibble(), out_filename); return(invisible(tibble::tibble())) }
  stars_df <- combined %>% dplyr::group_by(response) %>% dplyr::mutate(q = p.adjust(p.value, method = "fdr")) %>% dplyr::ungroup() %>% dplyr::mutate(stars = dplyr::case_when(q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ ""))
  star_cols <- intersect(c("response_All", "response_F", "response_M"), names(stars_df))
  if (length(star_cols) == 0) { save_table(tibble::tibble(), out_filename); return(invisible(tibble::tibble())) }
  stars_long <- stars_df %>% tidyr::pivot_longer(cols = dplyr::all_of(star_cols), names_to = "Group", values_to = "Measurement", values_drop_na = TRUE) %>% dplyr::select(Measurement, EXPOSURE, stars) %>% tidyr::pivot_wider(names_from = EXPOSURE, values_from = stars) %>% tidyr::pivot_longer(cols = -Measurement, names_to = "Exposure", values_to = "Stars") %>% dplyr::mutate(Exposure = stringr::str_trim(gsub(" - PBS_Control$", "", Exposure)))
  if (!("PBS_Control" %in% names(msd_summary))) { save_table(tibble::tibble(), out_filename); return(invisible(tibble::tibble())) }
  wide_no_pbs <- msd_summary %>% dplyr::select(-`PBS_Control`) %>% dplyr::arrange(Measurement) %>% dplyr::mutate(Analyte_Base = gsub("_.*", "", Measurement)) %>% dplyr::left_join(cytokine_llod, by = c("Analyte_Base" = "Analyte"))
  if (!"LLOD" %in% colnames(wide_no_pbs)) stop("LLOD column not found. Check cytokine_llod data.")
  long_vals <- wide_no_pbs %>% tidyr::pivot_longer(cols = -c(Measurement, Analyte_Base, LLOD), names_to = "Exposure", values_to = "Value_full") %>% dplyr::mutate(Value_num = as.numeric(sub("^(\\d*\\.?\\d+).*", "\\1", Value_full)))
  long_vals <- long_vals %>% dplyr::left_join(stars_long, by = c("Measurement", "Exposure")) %>% dplyr::mutate(Value_final = ifelse(!is.na(Stars) & !is.na(Value_num) & Value_num >= LLOD, paste0(Value_full, " ", Stars), Value_full))
  wide_final <- long_vals %>% dplyr::select(Measurement, Exposure, Value_final, LLOD, Analyte_Base) %>% tidyr::pivot_wider(names_from = Exposure, values_from = Value_final)
  pbs_numeric <- as.numeric(sub(" ±.*", "", pbs_control$`PBS_Control`))
  sig_table <- wide_final %>% dplyr::left_join(pbs_control, by = "Measurement") %>% dplyr::mutate(dplyr::across(dplyr::where(is.character) & !dplyr::all_of(c("Measurement", "PBS_Control")), ~ ifelse(pbs_numeric < epsilon, gsub("[*#~]+", "", .), .))) %>% dplyr::select(Measurement, `PBS_Control`, dplyr::everything(), -Analyte_Base, -LLOD)
  save_table(sig_table, out_filename)
  invisible(sig_table)
}

lmer_results_to_plot_format <- function(lmer_All, lmer_F, lmer_M,
                                        celltype, hormone,
                                        alpha = ALPHA_Q) {
  parse_rows <- function(df, sex_label) {
    if (is.null(df) || nrow(df) == 0) return(tibble::tibble())
    required <- c("contrast", "response", "estimate", "SE", "p.value")
    if (!all(required %in% names(df))) return(tibble::tibble())
    df %>% dplyr::mutate(SEX = sex_label, CELLTYPE = celltype, HORMONE = hormone, EXPOSURE = stringr::str_replace(contrast, " - PBS_Control$", ""), CYTOKINE = response) %>% dplyr::select(SEX, CELLTYPE, HORMONE, EXPOSURE, CYTOKINE, estimate, SE, p.value)
  }
  combined <- dplyr::bind_rows(parse_rows(lmer_All, "All"), parse_rows(lmer_F, "F"), parse_rows(lmer_M, "M"))
  if (nrow(combined) == 0) {
    return(tibble::tibble(SEX = character(), CELLTYPE = character(), HORMONE = character(), EXPOSURE = character(), CYTOKINE = character(), estimate = numeric(), SE = numeric(), p.value = numeric(), q = numeric(), sig = logical()))
  }
  combined %>% dplyr::group_by(CYTOKINE) %>% dplyr::mutate(q = p.adjust(p.value, method = "fdr")) %>% dplyr::ungroup() %>% dplyr::mutate(sig = q < alpha)
}

impute_lod_sqrt2 <- function(input_data, cols = NULL, cytokine_llod, zero_cutoff = ZERO_CUTOFF) {
  if (is.null(cols)) cols <- intersect(names(input_data), cytokine_llod$Analyte)
  llod_map <- stats::setNames(cytokine_llod$LLOD, cytokine_llod$Analyte)
  valid_cytokines <- cols[vapply(cols, function(col) {
    vals <- suppressWarnings(as.numeric(input_data[[col]]))
    if (all(is.na(vals))) return(FALSE)
    if ("SEX" %in% names(input_data)) {
      max_zero <- input_data %>% dplyr::mutate(.val = vals) %>% dplyr::group_by(SEX) %>% dplyr::summarise(p_zero = mean(is.na(.val) | .val <= 0), .groups = "drop") %>% dplyr::summarise(max_p = max(p_zero, na.rm = TRUE)) %>% dplyr::pull(max_p)
      return(is.finite(max_zero) && max_zero <= zero_cutoff)
    }
    mean(is.na(vals) | vals <= 0) <= zero_cutoff
  }, logical(1))]
  skipped_cytokines <- setdiff(cols, valid_cytokines)
  out <- input_data
  for (col in valid_cytokines) {
    llod_val <- llod_map[[col]]
    if (!is.null(llod_val) && is.finite(llod_val)) {
      vals <- suppressWarnings(as.numeric(out[[col]]))
      vals[!is.na(vals) & vals < llod_val] <- llod_val / sqrt(2)
      out[[col]] <- vals
    }
  }
  list(data = out, valid_cytokines = valid_cytokines, skipped_cytokines = skipped_cytokines)
}

summarize_to_wide <- function(data, measure_vars) {
  measure_vars <- intersect(measure_vars, names(data))
  if (length(measure_vars) == 0) return(tibble::tibble(Measurement = character(), PBS_Control = character()))
  long <- data %>% dplyr::select(dplyr::any_of(c("EXPOSURE", "SEX", measure_vars))) %>% tidyr::pivot_longer(cols = dplyr::all_of(measure_vars), names_to = "CYTOKINE", values_to = "VALUE")
  by_sex <- long %>% dplyr::group_by(CYTOKINE, SEX, EXPOSURE) %>% dplyr::summarise(mu = mean(VALUE, na.rm = TRUE), sd = stats::sd(VALUE, na.rm = TRUE), .groups = "drop") %>% dplyr::mutate(Measurement = paste0(CYTOKINE, "_", SEX))
  all_rows <- long %>% dplyr::group_by(CYTOKINE, EXPOSURE) %>% dplyr::summarise(mu = mean(VALUE, na.rm = TRUE), sd = stats::sd(VALUE, na.rm = TRUE), .groups = "drop") %>% dplyr::mutate(Measurement = paste0(CYTOKINE, "_All"))
  dplyr::bind_rows(by_sex, all_rows) %>%
    dplyr::mutate(mu = ifelse(is.nan(mu), NA_real_, mu), sd = ifelse(is.na(sd) | is.nan(sd), 0, sd), Value = sprintf("%.3f ± %.3f", mu, sd)) %>%
    dplyr::select(Measurement, EXPOSURE, Value) %>%
    tidyr::pivot_wider(names_from = EXPOSURE, values_from = Value) %>%
    dplyr::arrange(Measurement)
}
