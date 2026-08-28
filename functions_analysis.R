# =============================================================================
# functions_analysis.R
# Statistical models: LMER, ART, screening, significance tables
# =============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(emmeans)
  library(ARTool)
})

# ── Hormone label recoding ────────────────────────────────────────────────────

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

# ── log2FC long format ────────────────────────────────────────────────────────

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
    d <- d %>%
      dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }
  
  long <- d %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cytokine_cols),
                        names_to = "CYTOKINE", values_to = "VALUE") %>%
    dplyr::filter(!is.na(VALUE)) %>%
    dplyr::mutate(log2_val = log2(as.numeric(VALUE) + pseudocount))
  
  join_keys <- c("PATIENTCODE", "SEX", "CELLTYPE", "HORMONE", "CYTOKINE")
  if ("TIMEPOINT" %in% names(long)) {
    join_keys <- c(join_keys, "TIMEPOINT")
  }
  
  pbs <- long %>%
    dplyr::filter(EXPOSURE == pbs_level) %>%
    dplyr::select(PATIENTCODE, SEX, CELLTYPE, HORMONE,
                  dplyr::any_of("TIMEPOINT"),
                  CYTOKINE,
                  log2_pbs = log2_val)
  
  long %>%
    dplyr::filter(EXPOSURE != pbs_level) %>%
    dplyr::left_join(pbs, by = join_keys) %>%
    dplyr::mutate(log2FC = log2_val - log2_pbs)
}

# ── LMER screening: one exposure vs PBS, per CELLTYPE × HORMONE × SEX ────────
# NOTE: screen_one_exposure_lmer_log2() is SUPERSEDED and no longer called in
# 01_timecourse_analysis.Rmd.  The single source of statistical truth is now the
# exposure_lmer_pairwise() pipeline run inside run_lmer_chunk(), whose output
# is returned as lmer_plot_data and used for both significance tables AND plots.
# This function is kept here for reference only.

screen_one_exposure_lmer_log2 <- function(df, cytokine_cols, target_exposure,
                                          ctrl_level        = "PBS_Control",
                                          exclude_exposures = character(0),
                                          pseudocount       = 1e-6,
                                          group             = c("All","F","M"),
                                          alpha             = 0.1,
                                          emmeans_weights   = c("equal","proportional")) {
  group           <- match.arg(group)
  emmeans_weights <- match.arg(emmeans_weights)
  
  d0 <- df %>%
    dplyr::filter(!EXPOSURE %in% exclude_exposures,
                  EXPOSURE %in% c(ctrl_level, target_exposure)) %>%
    dplyr::mutate(
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE),
      CELLTYPE    = factor(CELLTYPE),
      HORMONE     = factor(HORMONE)
    )
  
  if (group %in% c("F","M"))
    d0 <- d0 %>% dplyr::filter(as.character(SEX) == group)
  if ("TIMEPOINT" %in% names(d0)) {
    d0 <- d0 %>%
      dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }
  
  if (ctrl_level %in% levels(d0$EXPOSURE))
    d0$EXPOSURE <- relevel(d0$EXPOSURE, ref = ctrl_level)
  
  # Empty result template, used whenever nothing survives screening
  # (e.g. too few samples in this SEX/CELLTYPE/HORMONE stratum). Returning
  # this instead of an empty bind_rows() output avoids a downstream
  # `dplyr::group_by()` error on a 0-row/0-column tibble.
  empty_result <- tibble::tibble(
    SEX      = character(),
    CELLTYPE = character(),
    HORMONE  = character(),
    EXPOSURE = character(),
    CYTOKINE = character(),
    estimate = numeric(),
    SE       = numeric(),
    p.value  = numeric(),
    q        = numeric(),
    sig      = logical()
  )
  
  get_contrast <- function(fit) {
    emm  <- emmeans::emmeans(fit, ~ EXPOSURE, weights = emmeans_weights)
    levs <- levels(emmeans::summary(emm)$EXPOSURE)
    v    <- if (identical(levs, c(ctrl_level, target_exposure))) c(-1, 1) else c(1, -1)
    contrast_list        <- list(v)
    names(contrast_list) <- paste0(target_exposure, " - ", ctrl_level)
    emmeans::contrast(emm, method = contrast_list, adjust = "none")
  }
  
  combos        <- d0 %>% dplyr::distinct(CELLTYPE, HORMONE)
  cytokine_cols <- intersect(cytokine_cols, names(d0))
  out           <- list()
  
  for (i in seq_len(nrow(combos))) {
    ct   <- combos$CELLTYPE[i]
    ho   <- combos$HORMONE[i]
    dsub <- d0 %>% dplyr::filter(CELLTYPE == ct, HORMONE == ho)
    if (!all(c(ctrl_level, target_exposure) %in% unique(dsub$EXPOSURE))) next
    
    for (cyt in cytokine_cols) {
      dat <- dsub %>% dplyr::filter(!is.na(.data[[cyt]]))
      if (nrow(dat) < 3) next
      if (!all(c(ctrl_level, target_exposure) %in% unique(dat$EXPOSURE))) next
      
      dat  <- dat %>% dplyr::mutate(resp = log2(.data[[cyt]] + pseudocount))
      tp_term <- if ("TIMEPOINT" %in% names(dat) &&
                     dplyr::n_distinct(dat$TIMEPOINT[!is.na(dat$TIMEPOINT)]) > 1) {
        " + TIMEPOINT"
      } else {
        ""
      }
      form <- as.formula(paste0("resp ~ EXPOSURE", tp_term, " + (1 | PATIENTCODE)"))
      
      fit <- try(lme4::lmer(form, data = dat), silent = TRUE)
      if (inherits(fit, "try-error")) next
      # Check for convergence / singular-fit warnings stored in the fit object
      if (length(lme4::isSingular(fit)) > 0 && lme4::isSingular(fit)) {
        warning("Singular fit for ", cyt, " in CELLTYPE=", ct, " HORMONE=", ho,
                "; estimates may be unreliable.")
      }
      
      con <- try(get_contrast(fit), silent = TRUE)
      if (inherits(con, "try-error")) next
      
      s <- as.data.frame(summary(con))
      out[[paste(group, ct, ho, cyt, sep = "|")]] <- tibble::tibble(
        SEX      = ifelse(group == "All", "All", group),
        CELLTYPE = as.character(ct),
        HORMONE  = as.character(ho),
        EXPOSURE = target_exposure,
        CYTOKINE = cyt,
        estimate = s$estimate[1],
        SE       = s$SE[1],
        p.value  = s$p.value[1]
      )
    }
  }
  
  result <- dplyr::bind_rows(out)
  
  if (nrow(result) == 0) {
    warning("screen_one_exposure_lmer_log2(): no cytokines survived screening for ",
            "group='", group, "', target_exposure='", target_exposure,
            "'. Returning an empty result.")
    return(empty_result)
  }
  
  result %>%
    dplyr::group_by(SEX, CELLTYPE, HORMONE, EXPOSURE) %>%
    dplyr::mutate(q   = p.adjust(p.value, method = "fdr"),
                  sig = q < alpha) %>%
    dplyr::ungroup()
}

# ── LMER pairwise: exposure vs PBS ────────────────────────────────────────────

exposure_lmer_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                   response_columns = NULL,
                                   ctrl_level = "PBS_Control",
                                   by_timepoint = FALSE) {
  if (group != "All" && "SEX" %in% names(data))
    data <- data %>% dplyr::filter(SEX == group)
  if ("TIMEPOINT" %in% names(data)) {
    data <- data %>%
      dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }
  
  stopifnot("EXPOSURE" %in% names(data), "PATIENTCODE" %in% names(data))
  
  data <- data %>%
    dplyr::mutate(
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  
  if (ctrl_level %in% levels(data$EXPOSURE))
    data$EXPOSURE <- relevel(data$EXPOSURE, ref = ctrl_level)
  
  if (is.null(response_columns))
    response_columns <- names(data)[sapply(data, is.numeric)]
  
  results_list <- list()
  
  for (resp in response_columns) {
    if (!resp %in% names(data)) next
    df <- data %>% dplyr::filter(!is.na(.data[[resp]]))
    if (length(unique(df$EXPOSURE)) < 2) next
    
    tp_term <- if ("TIMEPOINT" %in% names(df) &&
                   dplyr::n_distinct(df$TIMEPOINT[!is.na(df$TIMEPOINT)]) > 1) {
      " + TIMEPOINT"
    } else {
      ""
    }
    model <- try(
      lme4::lmer(as.formula(paste0(resp, " ~ EXPOSURE", tp_term, " + (1|PATIENTCODE)")),
                 data = df),
      silent = TRUE)
    if (inherits(model, "try-error")) { warning("Model failed for ", resp); next }
    
    has_timepoint <- by_timepoint &&
      "TIMEPOINT" %in% names(df) &&
      dplyr::n_distinct(df$TIMEPOINT[!is.na(df$TIMEPOINT)]) > 1
    emm_spec <- if (has_timepoint) {
      stats::as.formula("~ EXPOSURE | TIMEPOINT")
    } else {
      stats::as.formula("~ EXPOSURE")
    }
    emm      <- emmeans::emmeans(model, emm_spec, weights = "equal")
    ctrl_idx <- which(levels(data$EXPOSURE) == ctrl_level)
    pairwise <- emmeans::contrast(emm, method = "trt.vs.ctrl",
                                  ref = ctrl_idx, adjust = adjust_method)
    pairwise_df          <- as.data.frame(summary(pairwise))
    if (by_timepoint &&
        "TIMEPOINT" %in% names(df) &&
        !"TIMEPOINT" %in% names(pairwise_df)) {
      if (has_timepoint) {
        stop("TIMEPOINT missing from timepoint-specific contrast output for ", resp)
      }
      tp_values <- unique(as.character(stats::na.omit(df$TIMEPOINT)))
      if (length(tp_values) == 1) {
        pairwise_df$TIMEPOINT <- tp_values
      }
    }
    pairwise_df$response <- resp
    results_list[[resp]] <- pairwise_df
  }
  
  dplyr::bind_rows(results_list)
}

# ── LMER pairwise: sex × exposure interaction ─────────────────────────────────

interaction_lmer_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                      response_columns = NULL,
                                      ctrl_level = "PBS_Control") {
  if (group != "All" && "SEX" %in% names(data))
    data <- data %>% dplyr::filter(SEX == group)
  if ("TIMEPOINT" %in% names(data)) {
    data <- data %>%
      dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }
  
  stopifnot("EXPOSURE" %in% names(data), "PATIENTCODE" %in% names(data),
            "SEX"      %in% names(data))
  
  data <- data %>%
    dplyr::mutate(
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE),
      SEX         = factor(SEX)
    )
  
  if (ctrl_level %in% levels(data$EXPOSURE))
    data$EXPOSURE <- relevel(data$EXPOSURE, ref = ctrl_level)
  
  if (is.null(response_columns))
    response_columns <- names(data)[sapply(data, is.numeric)]
  
  results_list <- list()
  
  for (resp in response_columns) {
    if (!resp %in% names(data)) next
    df <- data %>% dplyr::filter(!is.na(.data[[resp]]))
    if (length(unique(df$EXPOSURE)) < 2 || length(unique(df$SEX)) < 2) next
    
    tp_term <- if ("TIMEPOINT" %in% names(df) &&
                   dplyr::n_distinct(df$TIMEPOINT[!is.na(df$TIMEPOINT)]) > 1) {
      " + TIMEPOINT"
    } else {
      ""
    }
    model <- try(
      lme4::lmer(
        as.formula(paste0(resp, " ~ EXPOSURE * SEX", tp_term, " + (1|PATIENTCODE)")),
        data = df),
      silent = TRUE)
    if (inherits(model, "try-error")) {
      warning("Interaction model failed for ", resp); next
    }
    
    emm_exp  <- emmeans::emmeans(model, ~ EXPOSURE, weights = "equal")
    ctrl_idx <- which(levels(data$EXPOSURE) == ctrl_level)
    pw_exp   <- emmeans::contrast(emm_exp, "trt.vs.ctrl",
                                  ref = ctrl_idx, adjust = adjust_method)
    exp_df   <- as.data.frame(summary(pw_exp))
    exp_df$type     <- "Exposure_vs_Control"
    exp_df$response <- resp
    
    emm_int  <- emmeans::emmeans(model, ~ SEX | EXPOSURE)
    pw_int   <- emmeans::contrast(emm_int, "pairwise",
                                  simple = "SEX", combine = TRUE)
    int_df   <- as.data.frame(summary(pw_int))
    int_df$type     <- "Sex_within_Exposure"
    int_df$response <- resp
    
    results_list[[resp]] <- dplyr::bind_rows(exp_df, int_df)
  }
  
  dplyr::bind_rows(results_list)
}

# ── LMER pairwise: sex × exposure interaction ─────────────────────────────────
timepoint_lmer_within_exposure_sex <- function(
    data,
    response,
    min_donors = 4,
    time_levels = c("4", "144")
) {
  required <- c(response, "TIMEPOINT", "PATIENTCODE", "CELLTYPE", "HORMONE", "EXPOSURE", "SEX")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("timepoint_lmer_within_exposure_sex(): missing columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  df <- data %>%
    dplyr::filter(
      !is.na(.data[[response]]),
      !is.na(TIMEPOINT),
      !is.na(PATIENTCODE),
      !is.na(CELLTYPE),
      !is.na(HORMONE),
      !is.na(EXPOSURE),
      !is.na(SEX)
    ) %>%
    dplyr::mutate(
      TIMEPOINT = factor(as.character(TIMEPOINT), levels = time_levels),
      PATIENTCODE = factor(PATIENTCODE),
      CELLTYPE = as.character(CELLTYPE),
      HORMONE  = as.character(HORMONE),
      EXPOSURE = as.character(EXPOSURE),
      SEX      = as.character(SEX)
    ) %>%
    dplyr::filter(TIMEPOINT %in% time_levels)
  
  # Split by strata: CELLTYPE × HORMONE × EXPOSURE × SEX
  strata <- df %>%
    dplyr::group_by(CELLTYPE, HORMONE, EXPOSURE, SEX) %>%
    dplyr::group_split(.keep = TRUE)
  
  one_group <- function(sdf) {
    key <- sdf %>%
      dplyr::slice(1) %>%
      dplyr::select(CELLTYPE, HORMONE, EXPOSURE, SEX)
    
    # donors that have BOTH timepoints
    donor_tp <- sdf %>%
      dplyr::group_by(PATIENTCODE) %>%
      dplyr::summarise(n_tp = dplyr::n_distinct(TIMEPOINT), .groups = "drop")
    
    keep_donors <- donor_tp %>%
      dplyr::filter(n_tp == length(time_levels)) %>%
      dplyr::pull(PATIENTCODE)
    
    sdf2 <- sdf %>% dplyr::filter(PATIENTCODE %in% keep_donors)
    
    n_donors <- dplyr::n_distinct(sdf2$PATIENTCODE)
    n_rows   <- nrow(sdf2)
    n_tp     <- dplyr::n_distinct(sdf2$TIMEPOINT)
    
    # defaults
    out <- key %>%
      dplyr::mutate(
        response = response,
        n_donors = n_donors,
        n_rows   = n_rows,
        estimate_144_vs_4 = NA_real_,
        p_value = NA_real_,
        converged = NA,
        note = NA_character_
      )
    
    if (n_donors < min_donors || n_tp < 2) {
      out$note <- paste0("insufficient data: donors=", n_donors, ", timepoints=", n_tp)
      return(out)
    }
    
    fml <- stats::as.formula(paste0("`", response, "` ~ TIMEPOINT + (1|PATIENTCODE)"))
    
    fit <- tryCatch(
      lme4::lmer(fml, data = sdf2, REML = FALSE),
      error = function(e) e
    )
    
    if (inherits(fit, "error")) {
      out$note <- paste("lmer_error:", fit$message)
      return(out)
    }
    
    # p-values via lmerTest if available; otherwise fallback to Wald z approximation
    tidy_fix <- tryCatch(
      broom.mixed::tidy(fit, effects = "fixed"),
      error = function(e) NULL
    )
    
    term_name <- paste0("TIMEPOINT", time_levels[2])
    if (is.null(tidy_fix) || !term_name %in% tidy_fix$term) {
      out$note <- "TIMEPOINT term not estimable"
      return(out)
    }
    
    est <- tidy_fix$estimate[tidy_fix$term == term_name][1]
    
    pval <- NA_real_
    if ("p.value" %in% names(tidy_fix)) {
      pval <- tidy_fix$p.value[tidy_fix$term == term_name][1]
    } else {
      # fallback using estimate/std.error
      se <- tidy_fix$std.error[tidy_fix$term == term_name][1]
      if (is.finite(se) && se > 0) {
        z <- est / se
        pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      }
    }
    
    out$estimate_144_vs_4 <- est
    out$p_value <- pval
    out$converged <- is.null(fit@optinfo$conv$lme4$messages)
    out$note <- ""
    out
  }
  
  res <- dplyr::bind_rows(lapply(strata, one_group)) %>%
    dplyr::mutate(
      p_adj_BH = stats::p.adjust(p_value, method = "BH")
    ) %>%
    dplyr::arrange(response, CELLTYPE, HORMONE, EXPOSURE, SEX)
  
  res
}

# ── ART pairwise ──────────────────────────────────────────────────────────────

exposure_art_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                  response_columns = NULL,
                                  ctrl_level = "PBS_Control") {
  filtered_data <- switch(group,
                          "M"   = data %>% dplyr::filter(EXPOSURE != "Untreated_Control", SEX == "M"),
                          "F"   = data %>% dplyr::filter(EXPOSURE != "Untreated_Control", SEX == "F"),
                          "All" = data %>% dplyr::filter(EXPOSURE != "Untreated_Control"),
                          stop("Invalid group. Choose 'All', 'M', or 'F'.")
  ) %>%
    dplyr::mutate(
      SEX         = factor(SEX),
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  if ("TIMEPOINT" %in% names(filtered_data)) {
    filtered_data <- filtered_data %>%
      dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }
  
  if (is.null(response_columns))
    response_columns <- names(filtered_data)[sapply(filtered_data, is.numeric)]
  
  results_list <- list()
  
  for (response in response_columns) {
    if (!response %in% names(filtered_data)) {
      warning("Column ", response, " not found"); next
    }
    
    tp_term <- if ("TIMEPOINT" %in% names(filtered_data) &&
                   dplyr::n_distinct(filtered_data$TIMEPOINT[!is.na(filtered_data$TIMEPOINT)]) > 1) {
      " + TIMEPOINT"
    } else {
      ""
    }
    
    formula <- if (group == "All") {
      as.formula(paste0(response, " ~ EXPOSURE * SEX", tp_term, " + (1|PATIENTCODE)"))
    } else {
      as.formula(paste0(response, " ~ EXPOSURE", tp_term, " + (1|PATIENTCODE)"))
    }
    
    m.art <- try(ARTool::art(formula, data = filtered_data), silent = TRUE)
    if (inherits(m.art, "try-error")) {
      warning("ART failed for ", response); next
    }
    
    anova_res   <- anova(m.art)
    exposure_p  <- {
      r <- anova_res[grepl("^EXPOSURE$", anova_res[[1]], ignore.case = TRUE), ]
      if (nrow(r) > 0) r[["Pr(>F)"]][1] else NA
    }
    interaction_p <- if (group == "All") {
      r <- anova_res[grepl("EXPOSURE:SEX", anova_res[[1]], ignore.case = TRUE), ]
      if (nrow(r) > 0) r[["Pr(>F)"]][1] else NA
    } else NA
    
    # art.con() is the correct way to get pairwise contrasts from an ART model
    # (artlm.con() + emmeans() produces unreliable d.f. for interaction contrasts).
    pairwise_con <- try(
      ARTool::art.con(m.art, "EXPOSURE", adjust = adjust_method),
      silent = TRUE
    )
    if (inherits(pairwise_con, "try-error")) {
      warning("art.con() failed for ", response); next
    }
    ctrl_idx <- which(levels(filtered_data$EXPOSURE) == ctrl_level)
    pairwise <- emmeans::contrast(pairwise_con, "trt.vs.ctrl",
                                  ref = ctrl_idx, adjust = adjust_method)
    
    res               <- as.data.frame(pairwise)
    res$response      <- response
    res$exposure_p    <- exposure_p
    res$interaction_p <- interaction_p
    results_list[[response]] <- res
  }
  
  dplyr::bind_rows(results_list)
}

# ── Convenience wrapper: run all 4 CELLTYPE × HORMONE strata ─────────────────

run_lmer_chunk <- function(label, celltype_filter, hormone_filter,
                           sala_full     = NULL,
                           llod_table    = cytokine_llod,
                           zero_co       = ZERO_CUTOFF,
                           alpha_q       = ALPHA_Q,
                           trend_a       = TREND_ALPHA) {
  cat("\n── LMER:", label, "──\n")
  
  if (is.null(sala_full))
    stop("sala_full must be provided as a data frame containing the full WSTC dataset")
  
  d_sub <- sala_full %>%
    dplyr::filter(CELLTYPE == celltype_filter, HORMONE == hormone_filter)
  
  imp_res    <- impute_lod_sqrt2(d_sub, cytokine_llod = llod_table, zero_cutoff = zero_co)
  d_imp      <- imp_res$data
  valid_cyts <- imp_res$valid_cytokines
  
  if (length(valid_cyts) == 0) {
    warning("run_lmer_chunk(", label, "): no valid cytokines after ZERO_CUTOFF filtering.")
    empty_sig <- tibble::tibble(Measurement = character(), PBS_Control = character())
    return(invisible(list(
      sig_table = empty_sig,
      lmer_interaction = tibble::tibble(),
      lmer_plot_data = tibble::tibble(
        SEX = character(), CELLTYPE = character(), HORMONE = character(),
        EXPOSURE = character(), CYTOKINE = character(),
        estimate = numeric(), SE = numeric(), p.value = numeric(),
        q = numeric(), sig = logical()
      )
    )))
  }
  
  d_filt <- d_imp %>%
    dplyr::filter(EXPOSURE != "Untreated_Control") %>%
    dplyr::mutate(EXPOSURE = factor(EXPOSURE), SEX = factor(SEX), PATIENTCODE = factor(PATIENTCODE))
  
  if ("TIMEPOINT" %in% names(d_filt)) {
    d_filt <- d_filt %>%
      dplyr::mutate(TIMEPOINT = factor(TIMEPOINT, levels = TIMEPOINT_LEVELS))
  }
  
  msd_sum <- summarize_to_wide(d_filt, measure_vars = valid_cyts)
  
  pbs_ctl <- msd_sum %>%
    dplyr::select(Measurement, `PBS_Control`) %>%
    dplyr::arrange(Measurement)
  
  lmer_All <- exposure_lmer_pairwise(d_filt, "All", "fdr", valid_cyts)
  lmer_F   <- exposure_lmer_pairwise(d_filt, "F",   "fdr", valid_cyts)
  lmer_M   <- exposure_lmer_pairwise(d_filt, "M",   "fdr", valid_cyts)
  lmer_plot_All <- exposure_lmer_pairwise(d_filt, "All", "fdr", valid_cyts,
                                          by_timepoint = TRUE)
  lmer_plot_F   <- exposure_lmer_pairwise(d_filt, "F",   "fdr", valid_cyts,
                                          by_timepoint = TRUE)
  lmer_plot_M   <- exposure_lmer_pairwise(d_filt, "M",   "fdr", valid_cyts,
                                          by_timepoint = TRUE)
  lmer_All <- exposure_lmer_pairwise_by_timepoint(d_filt, "All", "fdr", valid_cyts, ctrl_level = PBS_LEVEL)
  lmer_F   <- exposure_lmer_pairwise_by_timepoint(d_filt, "F",   "fdr", valid_cyts, ctrl_level = PBS_LEVEL)
  lmer_M   <- exposure_lmer_pairwise_by_timepoint(d_filt, "M",   "fdr", valid_cyts, ctrl_level = PBS_LEVEL)
  lmer_int <- interaction_lmer_pairwise(d_filt, "All", "fdr", valid_cyts)
  
  # Convert LMER results to plot-compatible format (used by cytokine dotplots
  # and bar plots instead of the retired screen_one_exposure_lmer_log2()).
  lmer_plot_data <- lmer_results_to_plot_format(
    lmer_plot_All, lmer_plot_F, lmer_plot_M,
    celltype = celltype_filter,
    hormone  = hormone_filter,
    alpha    = alpha_q
  )
  
  out_csv <- paste0("MSD_WSTC_", label, "_lmer.csv")
  sig_tbl <- build_lmer_sig_table(
    lmer_All      = lmer_All,
    lmer_F        = lmer_F,
    lmer_M        = lmer_M,
    msd_summary   = msd_sum,
    pbs_control   = pbs_ctl,
    cytokine_llod = llod_table,
    out_filename  = out_csv,
    alpha         = alpha_q,
    trend_alpha   = trend_a
  )
  
  cat("  Saved:", out_csv, "\n")
  # Return significance table, interaction results, and plot-ready LMER data
  invisible(list(sig_table = sig_tbl, lmer_interaction = lmer_int,
                 lmer_plot_data = lmer_plot_data))
}

build_lmer_sig_table <- function(lmer_All, lmer_F, lmer_M,
                                 msd_summary, pbs_control,
                                 cytokine_llod,
                                 out_filename,
                                 alpha       = 0.1,
                                 trend_alpha = 0.2,
                                 epsilon     = 1e-5) {
  
  # Step 6: format with group suffixes
  fmt <- function(df, grp) {
    col <- paste0("response_", grp)
    df %>%
      dplyr::mutate(
        !!col    := paste0(response, "_", grp),
        EXPOSURE  = stringr::str_replace(contrast, " - PBS_Control", "")
      )
  }
  combined <- dplyr::bind_rows(fmt(lmer_All, "All"),
                               fmt(lmer_F,   "F"),
                               fmt(lmer_M,   "M"))
  
  # Step 7: apply FDR correction per group, then assign direction-aware stars
  # Stars are based on the FDR-adjusted q-value, not the raw p.value.
  stars_df <- combined %>%
    dplyr::group_by(response) %>%
    dplyr::mutate(q = p.adjust(p.value, method = "fdr")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(stars = dplyr::case_when(
      q < 0.001 ~ "***",
      q < 0.01  ~ "**",
      q < 0.05  ~ "*",
      TRUE      ~ ""
    ))
  
  # Steps 8-9: pivot stars long
  stars_long <- stars_df %>%
    tidyr::pivot_longer(
      cols           = dplyr::any_of(c("response_All","response_F","response_M")),
      names_to       = "Group",
      values_to      = "Measurement",
      values_drop_na = TRUE
    ) %>%
    dplyr::select(Measurement, EXPOSURE, stars) %>%
    tidyr::pivot_wider(names_from = EXPOSURE, values_from = stars) %>%
    tidyr::pivot_longer(cols = -Measurement,
                        names_to = "Exposure", values_to = "Stars") %>%
    dplyr::mutate(Exposure = stringr::str_trim(
      gsub(" - PBS_Control$", "", Exposure)))
  
  # Step 10: summary without PBS
  wide_no_pbs <- msd_summary %>%
    dplyr::select(-`PBS_Control`) %>%
    dplyr::arrange(Measurement) %>%
    dplyr::mutate(Analyte_Base = gsub("_.*", "", Measurement)) %>%
    dplyr::left_join(cytokine_llod, by = c("Analyte_Base" = "Analyte"))
  
  if (!"LLOD" %in% colnames(wide_no_pbs))
    stop("LLOD column not found. Check cytokine_llod data.")
  
  # Step 11: pivot long + numeric values
  long_vals <- wide_no_pbs %>%
    tidyr::pivot_longer(
      cols      = -c(Measurement, Analyte_Base, LLOD),
      names_to  = "Exposure",
      values_to = "Value_full"
    ) %>%
    dplyr::mutate(
      Value_num = as.numeric(sub("^(\\d*\\.?\\d+).*", "\\1", Value_full))
    )
  
  # Step 12: merge stars + LLOD check
  long_vals <- long_vals %>%
    dplyr::left_join(stars_long, by = c("Measurement","Exposure")) %>%
    dplyr::mutate(Value_final = ifelse(
      !is.na(Stars) & !is.na(Value_num) & Value_num >= LLOD,
      paste0(Value_full, " ", Stars),
      Value_full
    ))
  
  # Step 13: pivot wide
  wide_final <- long_vals %>%
    dplyr::select(Measurement, Exposure, Value_final, LLOD, Analyte_Base) %>%
    tidyr::pivot_wider(names_from = Exposure, values_from = Value_final)
  
  # Step 14: add PBS, apply epsilon, export
  pbs_numeric <- as.numeric(sub(" \u00b1.*", "", pbs_control$`PBS_Control`))
  
  sig_table <- wide_final %>%
    dplyr::left_join(pbs_control, by = "Measurement") %>%
    dplyr::mutate(dplyr::across(
      dplyr::where(is.character) &
        !dplyr::all_of(c("Measurement","PBS_Control")),
      ~ ifelse(pbs_numeric < epsilon, gsub("[*#~]+", "", .), .)
    )) %>%
    dplyr::select(Measurement, `PBS_Control`, dplyr::everything(),
                  -Analyte_Base, -LLOD)
  
  save_table(sig_table, out_filename)
  invisible(sig_table)
}

# ── Convert exposure_lmer_pairwise() output to plot-compatible format ─────────
# Used internally by run_lmer_chunk() to produce lmer_plot_data.
# FDR correction is applied per cytokine and, when present, per TIMEPOINT.
# 
exposure_lmer_pairwise_by_timepoint <- function(
    data,
    group = "All",
    adjust_method = "fdr",
    response_columns = NULL,
    ctrl_level = PBS_LEVEL,
    time_levels = c("4", "144")
) {
  # ---------------------------------------------------------------------------
  # Pair-matched, within-timepoint exposure vs PBS contrasts
  # - For each TIMEPOINT and each treated exposure:
  #   keep only donor-strata (PATIENTCODE,CELLTYPE,HORMONE,SEX,TIMEPOINT)
  #   where BOTH ctrl_level and treated exposure are present.
  # - Fit: response ~ EXPOSURE + (1|PATIENTCODE)
  # - Return trt-vs-ctrl contrast only for that treated exposure.
  # ---------------------------------------------------------------------------
  
  required_cols <- c("TIMEPOINT", "EXPOSURE", "PATIENTCODE")
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("exposure_lmer_pairwise_by_timepoint(): missing columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Optional SEX filter for group-specific runs
  d0 <- data
  if (group != "All") {
    if (!"SEX" %in% names(d0)) {
      stop("exposure_lmer_pairwise_by_timepoint(): SEX column required when group != 'All'")
    }
    d0 <- d0 %>% dplyr::filter(as.character(SEX) == group)
  }
  
  # Standardize core factors
  d0 <- d0 %>%
    dplyr::mutate(
      TIMEPOINT   = factor(as.character(TIMEPOINT), levels = time_levels),
      EXPOSURE    = as.character(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    ) %>%
    dplyr::filter(!is.na(TIMEPOINT), !is.na(EXPOSURE), !is.na(PATIENTCODE))
  
  if (is.null(response_columns)) {
    response_columns <- names(d0)[vapply(d0, is.numeric, logical(1))]
  } else {
    response_columns <- intersect(response_columns, names(d0))
  }
  
  # Helper: keep only matched donor-strata for one treated exposure vs PBS
  match_pairs_for_exposure <- function(df, exposure, pbs = ctrl_level) {
    keys <- intersect(c("PATIENTCODE", "CELLTYPE", "HORMONE", "SEX", "TIMEPOINT"), names(df))
    
    keep <- df %>%
      dplyr::filter(EXPOSURE %in% c(pbs, exposure)) %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(c(keys, "EXPOSURE")))) %>%
      dplyr::mutate(present = 1L) %>%
      tidyr::pivot_wider(
        names_from  = EXPOSURE,
        values_from = present,
        values_fill = 0
      )
    
    if (!(pbs %in% names(keep)) || !(exposure %in% names(keep))) {
      return(df[0, , drop = FALSE])
    }
    
    keep <- keep %>%
      dplyr::filter(.data[[pbs]] == 1L, .data[[exposure]] == 1L) %>%
      dplyr::select(dplyr::all_of(keys))
    
    df %>%
      dplyr::filter(EXPOSURE %in% c(pbs, exposure)) %>%
      dplyr::inner_join(keep, by = keys)
  }
  
  out <- list()
  idx <- 1L
  
  # loop by timepoint
  for (tp in levels(d0$TIMEPOINT)) {
    d_tp <- d0 %>% dplyr::filter(as.character(TIMEPOINT) == tp)
    if (nrow(d_tp) == 0) next
    
    # treated exposures present at this timepoint (excluding control)
    trt_levels <- setdiff(sort(unique(as.character(d_tp$EXPOSURE))), ctrl_level)
    if (length(trt_levels) == 0) next
    
    for (trt in trt_levels) {
      # strict pair-match for THIS exposure vs PBS at THIS timepoint
      d_pair <- match_pairs_for_exposure(d_tp, exposure = trt, pbs = ctrl_level)
      if (nrow(d_pair) == 0) next
      
      # Ensure both levels exist post-match
      ex_levels <- unique(as.character(d_pair$EXPOSURE))
      if (!all(c(ctrl_level, trt) %in% ex_levels)) next
      
      for (resp in response_columns) {
        df_resp <- d_pair %>% dplyr::filter(!is.na(.data[[resp]]))
        if (nrow(df_resp) < 3) next
        if (dplyr::n_distinct(df_resp$EXPOSURE) < 2) next
        
        # Re-factor and set control reference within this pair
        df_resp <- df_resp %>%
          dplyr::mutate(EXPOSURE = factor(EXPOSURE, levels = c(ctrl_level, trt)))
        
        # model
        fit <- try(
          lme4::lmer(
            stats::as.formula(paste0("`", resp, "` ~ EXPOSURE + (1|PATIENTCODE)")),
            data = df_resp,
            REML = FALSE
          ),
          silent = TRUE
        )
        if (inherits(fit, "try-error")) next
        
        # contrast trt vs ctrl
        emm <- try(emmeans::emmeans(fit, ~ EXPOSURE, weights = "equal"), silent = TRUE)
        if (inherits(emm, "try-error")) next
        
        # since levels are explicitly c(ctrl, trt), reference is 1
        pw <- try(
          emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1, adjust = adjust_method),
          silent = TRUE
        )
        if (inherits(pw, "try-error")) next
        
        s <- as.data.frame(summary(pw))
        if (nrow(s) == 0) next
        
        # keep only the target treated row if multiple show up
        s <- s %>%
          dplyr::filter(grepl(paste0("^", trt, " - "), contrast) |
                          grepl(paste0("^", trt, "/"), contrast) |
                          grepl(trt, contrast))
        
        if (nrow(s) == 0) next
        
        # donor count in matched pair set (for diagnostics)
        n_donors <- dplyr::n_distinct(df_resp$PATIENTCODE)
        
        out[[idx]] <- tibble::tibble(
          contrast  = s$contrast[1],
          estimate  = s$estimate[1],
          SE        = s$SE[1],
          df        = if ("df" %in% names(s)) s$df[1] else NA_real_,
          t.ratio   = if ("t.ratio" %in% names(s)) s$`t.ratio`[1] else NA_real_,
          p.value   = s$p.value[1],
          response  = resp,
          TIMEPOINT = as.character(tp),
          SEX       = ifelse(group == "All", "All", group),
          n_donors  = n_donors
        )
        idx <- idx + 1L
      }
    }
  }
  
  dplyr::bind_rows(out)
}

# FDR correction matches build_lmer_sig_table() exactly:
#   group_by(CYTOKINE) %>% p.adjust(p.value, method = "fdr")
# so q-values are IDENTICAL between significance tables and plots.

lmer_results_to_plot_format <- function(lmer_All, lmer_F, lmer_M,
                                        celltype, hormone,
                                        alpha = ALPHA_Q) {
  parse_rows <- function(df, sex_label) {
    df %>%
      dplyr::mutate(
        SEX       = sex_label,
        CELLTYPE  = celltype,
        HORMONE   = hormone,
        TIMEPOINT = as.character(TIMEPOINT),
        EXPOSURE  = stringr::str_replace(contrast, paste0(" - ", stringr::fixed(PBS_LEVEL), "$"), ""),
        CYTOKINE  = response
      ) %>%
      dplyr::select(SEX, CELLTYPE, HORMONE, dplyr::any_of("TIMEPOINT"),
                    EXPOSURE, CYTOKINE,
      dplyr::select(SEX, CELLTYPE, HORMONE, TIMEPOINT, EXPOSURE, CYTOKINE,
                    estimate, SE, p.value)
  }
  
  combined <- dplyr::bind_rows(
    parse_rows(lmer_All, "All"),
    parse_rows(lmer_F,   "F"),
    parse_rows(lmer_M,   "M")
  )
  
  group_vars <- c("CYTOKINE", intersect("TIMEPOINT", names(combined)))

  combined %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
  combined %>%
    dplyr::group_by(SEX, CELLTYPE, HORMONE, TIMEPOINT, CYTOKINE) %>%
    dplyr::mutate(q = p.adjust(p.value, method = "fdr")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sig = q < alpha)
}

# ── LLOD imputation + summary helpers ─────────────────────────────────────────

impute_lod_sqrt2 <- function(input_data, cols = NULL, cytokine_llod, zero_cutoff = ZERO_CUTOFF) {
  if (is.null(cols)) {
    cols <- intersect(names(input_data), cytokine_llod$Analyte)
  }
  
  llod_map <- stats::setNames(cytokine_llod$LLOD, cytokine_llod$Analyte)
  valid_cytokines <- cols[vapply(cols, function(col) {
    vals <- suppressWarnings(as.numeric(input_data[[col]]))
    if (all(is.na(vals))) return(FALSE)
    
    if ("SEX" %in% names(input_data)) {
      max_zero <- input_data %>%
        dplyr::mutate(.val = vals) %>%
        dplyr::group_by(SEX) %>%
        dplyr::summarise(p_zero = mean(is.na(.val) | .val <= 0), .groups = "drop") %>%
        dplyr::summarise(max_p = max(p_zero, na.rm = TRUE)) %>%
        dplyr::pull(max_p)
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
  if (length(measure_vars) == 0) {
    warning("summarize_to_wide(): no measure_vars found in data; returning empty summary.")
    return(tibble::tibble(Measurement = character()))
  }
  
  long <- data %>%
    dplyr::select(dplyr::any_of(c("EXPOSURE", "SEX", measure_vars))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(measure_vars),
      names_to = "CYTOKINE",
      values_to = "VALUE"
    )
  
  by_sex <- long %>%
    dplyr::group_by(CYTOKINE, SEX, EXPOSURE) %>%
    dplyr::summarise(
      mu = mean(VALUE, na.rm = TRUE),
      sd = stats::sd(VALUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Measurement = paste0(CYTOKINE, "_", SEX))
  
  all_rows <- long %>%
    dplyr::group_by(CYTOKINE, EXPOSURE) %>%
    dplyr::summarise(
      mu = mean(VALUE, na.rm = TRUE),
      sd = stats::sd(VALUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Measurement = paste0(CYTOKINE, "_All"))
  
  dplyr::bind_rows(by_sex, all_rows) %>%
    dplyr::mutate(
      mu = ifelse(is.nan(mu), NA_real_, mu),
      sd = ifelse(is.na(sd) | is.nan(sd), 0, sd),
      Value = sprintf("%.3f \u00b1 %.3f", mu, sd)
    ) %>%
    dplyr::select(Measurement, EXPOSURE, Value) %>%
    tidyr::pivot_wider(names_from = EXPOSURE, values_from = Value) %>%
    dplyr::arrange(Measurement)
}
