# =============================================================================
# functions_plots.R
# All ggplot2 helpers: dotplots, bar plots, histograms, boxplots,
# gene expression bars, GSEA barplots
# =============================================================================

UP_COLOR_DEFAULT   <- "#D7191C"
DOWN_COLOR_DEFAULT <- "#2C7BB6"
ESTRADIOL_FILL_DEFAULT <- "#E8F5E9"

add_rescaled_size <- function(df, value_col,
                              clamp_quantile = 0.95,
                              transform = c("sqrt","log1p","none"),
                              out_col = "size_val") {
  transform <- match.arg(transform)
  mag  <- abs(df[[value_col]])
  cap  <- stats::quantile(mag, probs = clamp_quantile, na.rm = TRUE, names = FALSE)
  if (!is.finite(cap) || cap <= 0) cap <- max(mag, na.rm = TRUE)
  mag2 <- pmin(mag, cap)
  df[[out_col]] <- switch(transform,
                          sqrt  = sqrt(mag2),
                          log1p = log1p(mag2),
                          none  = mag2)
  df
}

make_dotplot <- function(df, x_var, title_str,
                         cyt_levels, border_color,
                         up_color   = UP_COLOR_DEFAULT,
                         down_color = DOWN_COLOR_DEFAULT,
                         show_y = TRUE) {
  facet_rows <- if ("TIMEPOINT" %in% names(df)) ggplot2::vars(CELLTYPE, TIMEPOINT) else ggplot2::vars(CELLTYPE)
  fill_scale  <- ggplot2::scale_fill_manual(name = "Direction\n(vs matched PBS)", values = c("Up" = up_color, "Down" = down_color, "Zero" = "grey70"), na.value = "grey70")
  color_scale <- ggplot2::scale_color_manual(name = "Direction\n(vs matched PBS)", values = c("Up" = up_color, "Down" = down_color, "Zero" = "grey70"), na.value = "grey70")
  size_scale  <- ggplot2::scale_size_continuous(name = "Effect size\n(sqrt|mean log2FC|)", range = c(1, 8), limits = c(0, NA))

  if (!x_var %in% names(df)) stop("make_dotplot(): x_var not found in data: ", x_var)
  if (nrow(df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = title_str, subtitle = "No significant cytokines to plot") +
        ggplot2::theme_void()
    )
  }

  ref_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA, colour = border_color, linewidth = 1.2),
      panel.ontop      = TRUE,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = border_color, colour = border_color, linewidth = 1.2),
      strip.text       = ggplot2::element_text(color = "white", face = "bold", size = 11),
      axis.title.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(face = "bold", size = 11),
      axis.text.y      = ggplot2::element_text(face = "bold", size = 10),
      axis.title.y     = ggplot2::element_text(size = 11),
      legend.position  = "right",
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
    )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = factor(CYTOKINE, levels = rev(cyt_levels)))) +
    ggplot2::geom_tile(ggplot2::aes(fill = panel_fill), width = Inf, height = 1,
                       alpha = 1, show.legend = FALSE) +
    ggplot2::scale_fill_identity() +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(ggplot2::aes(fill = direction, color = direction, size = size_val),
                        shape = 21, stroke = 1.5, alpha = 0.9) +
    fill_scale + color_scale + size_scale +
    ggplot2::guides(
      color = "none",
      fill  = ggplot2::guide_legend(title = "Direction\n(vs matched PBS)"),
      size  = ggplot2::guide_legend(title = "Effect size\n(sqrt|mean log2FC|)")
    ) +
    ggplot2::facet_grid(rows = facet_rows, cols = ggplot2::vars(HORMONE), drop = FALSE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.8, 0.8)), drop = TRUE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ref_theme +
    ggplot2::labs(title = title_str, y = "Cytokine")

  if (!show_y) {
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }
  p
}

make_cytokine_barplot <- function(summary_raw, d_raw, cyt,
                                  target_exposure, exp_short,
                                  e2_label, hormone_levels,
                                  border_color,
                                  panel_fill_none,
                                  alpha_q     = ALPHA_Q,
                                  trend_alpha = TREND_ALPHA,
                                  up_color    = UP_COLOR_DEFAULT,
                                  down_color  = DOWN_COLOR_DEFAULT,
                                  estradiol_fill = ESTRADIOL_FILL_DEFAULT) {
  facet_rows <- if ("TIMEPOINT" %in% names(summary_raw)) ggplot2::vars(CELLTYPE, TIMEPOINT) else ggplot2::vars(CELLTYPE)

  y_max     <- max(summary_raw$mu + summary_raw$sd, na.rm = TRUE)
  y_bracket <- y_max * 1.10
  bg_df     <- summary_raw %>% dplyr::distinct(CELLTYPE, HORMONE, panel_fill)

  bracket_df <- summary_raw %>%
    dplyr::filter(EXPOSURE == target_exposure, sig_label != "") %>%
    dplyr::mutate(
      x_left   = as.numeric(x_label) - 1,
      x_right  = as.numeric(x_label),
      x_center = as.numeric(x_label) - 0.5,
      y_top    = y_bracket
    )

  bar_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA, colour = border_color, linewidth = 1),
      panel.ontop      = TRUE,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = border_color, colour = border_color, linewidth = 1),
      strip.text       = ggplot2::element_text(color = "white", face = "bold", size = 11),
      axis.title.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(face = "bold", size = 9, angle = 30, hjust = 1),
      axis.text.y      = ggplot2::element_text(face = "bold", size = 10),
      axis.title.y     = ggplot2::element_text(size = 11),
      legend.position  = "right",
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
    )

  p <- ggplot2::ggplot(summary_raw, ggplot2::aes(x = x_label, y = mu)) +
    ggplot2::geom_rect(
      data = bg_df, ggplot2::aes(fill = panel_fill),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE, show.legend = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_hline(yintercept = pretty(c(0, y_max)), color = "white", linewidth = 0.5) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_col(
      data = summary_raw %>% dplyr::filter(SEX == "M"),
      ggplot2::aes(x = x_label, y = mu, fill = bar_fill, color = bar_color),
      width = 0.9, linewidth = 0.7, alpha = 1) +
    ggplot2::geom_col(
      data = summary_raw %>% dplyr::filter(SEX == "F"),
      ggplot2::aes(x = x_label, y = mu, fill = bar_fill, color = bar_color),
      width = 0.9, linewidth = 0.7, alpha = 0) +
    ggplot2::scale_fill_identity(guide = "none") +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::geom_errorbar(
      ggplot2::aes(x = x_label, ymin = mu - sd, ymax = mu + sd),
      width = 0.2, linewidth = 0.6, color = "grey20") +
    { if (nrow(bracket_df) > 0) list(
      ggplot2::geom_segment(data = bracket_df, ggplot2::aes(x = x_left, xend = x_right, y = y_top, yend = y_top, color = sig_color), linewidth = 0.7, inherit.aes = FALSE, show.legend = FALSE),
      ggplot2::geom_segment(data = bracket_df, ggplot2::aes(x = x_left, xend = x_left, y = y_top, yend = y_top - y_max * 0.04, color = sig_color), linewidth = 0.7, inherit.aes = FALSE, show.legend = FALSE),
      ggplot2::geom_segment(data = bracket_df, ggplot2::aes(x = x_right, xend = x_right, y = y_top, yend = y_top - y_max * 0.04, color = sig_color), linewidth = 0.7, inherit.aes = FALSE, show.legend = FALSE),
      ggplot2::geom_text(data = bracket_df, ggplot2::aes(x = x_center, y = y_top + y_max * 0.03, label = sig_label, color = sig_color), size = 6, fontface = "bold", hjust = 0.5, vjust = 0, inherit.aes = FALSE, show.legend = FALSE)
    ) else list() } +
    ggplot2::geom_jitter(
      data = d_raw,
      ggplot2::aes(x = x_label, y = VALUE, fill = point_fill, color = point_color),
      shape = 21, width = 0.12, size = 2.5,
      stroke = 1.0, alpha = 0.95, inherit.aes = FALSE) +
    ggplot2::scale_fill_identity(guide  = "none") +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
    ggplot2::facet_grid(rows = facet_rows, cols = ggplot2::vars(HORMONE)) +
    bar_theme +
    ggplot2::labs(title = paste0(cyt, " — ", exp_short, " vs Control | M = filled, F = transparent"),
                  y = paste0(cyt, " (pg/mL)"))

  p
}

make_cytokine_histogram <- function(d_raw, cyt, target_exposure) {
  facet_rows <- if ("TIMEPOINT" %in% names(d_raw)) ggplot2::vars(CELLTYPE, TIMEPOINT) else ggplot2::vars(CELLTYPE)
  ggplot2::ggplot(
    d_raw %>% dplyr::filter(EXPOSURE == target_exposure),
    ggplot2::aes(x = VALUE, fill = SEX)
  ) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.65,
                            bins = 20, color = "grey40") +
    ggplot2::facet_grid(rows = facet_rows, cols = ggplot2::vars(HORMONE)) +
    ggplot2::scale_fill_manual(values = c("M" = "lemonchiffon", "F" = "white"), name = "Sex") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = paste("Raw Values —", cyt, "|", target_exposure),
                  x = paste(cyt, "(pg/mL)"), y = "Count")
}

save_plot <- function(filename, plot = ggplot2::last_plot(),
                      width = 10, height = 8, dpi = 300,
                      path = PATH_OUTPUT_FIGURES, bg = "white") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename = file.path(path, filename), plot = plot,
                  width = width, height = height, dpi = dpi,
                  units = "in", bg = bg)
  cat("✓ Saved plot:", file.path(path, filename), "\n")
  invisible(TRUE)
}

save_table <- function(data, filename, path = PATH_OUTPUT_TABLES) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(data, file = file.path(path, filename))
  cat("✓ Saved table:", file.path(path, filename), "\n")
  invisible(TRUE)
}
