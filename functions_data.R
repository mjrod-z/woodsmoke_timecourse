# R/functions_data.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

extract_timepoint_suffix <- function(sample_id) {
  sample_id <- trimws(as.character(sample_id))
  has_suffix <- grepl("(_144|_4)$", sample_id)
  out <- rep(NA_character_, length(sample_id))
  out[has_suffix] <- sub("^.*_(144|4)$", "\\1", sample_id[has_suffix])
  out
}

# ============================================================================
# SALA-SPECIFIC DATA FUNCTIONS
# ============================================================================

average_nonzero_by_sample <- function(data) {
  # Average replicate measurements within each sample.
  # All numeric values (including zeros) are included in the mean so that the
  # downstream LOD-imputation step operates on unbiased sample averages.
  # NA values are excluded via na.rm = TRUE.
  
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  id_col <- names(data)[!names(data) %in% numeric_cols][1]
  
  if (is.na(id_col)) {
    stop("average_nonzero_by_sample(): no sample ID column found")
  }
  
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(numeric_cols),
                    ~ mean(., na.rm = TRUE)),
      .groups = "drop"
    )
}

build_sala_full <- function(averaged_data, metadata) {
  if (!"SAMPLEID" %in% names(averaged_data)) {
    stop("build_sala_full(): averaged_data must contain SAMPLEID.")
  }
  if (!"SAMPLEID" %in% names(metadata)) {
    stop("build_sala_full(): metadata must contain SAMPLEID.")
  }
  
  norm_id <- function(x) {
    x <- trimws(as.character(x))
    x <- toupper(x)
    x <- sub("(_144|_4)$", "", x)
    x <- gsub("[^A-Z0-9]", "", x)
    x
  }
  
  avg <- averaged_data %>%
    dplyr::mutate(
      SAMPLEID      = as.character(SAMPLEID),
      TIMEPOINT     = extract_timepoint_suffix(SAMPLEID),
      SAMPLEID_BASE = sub("(_144|_4)$", "", SAMPLEID),
      JOIN_ID       = norm_id(SAMPLEID_BASE)
    )
  
  meta <- metadata %>%
    dplyr::mutate(
      META_SAMPLEID = as.character(SAMPLEID),
      JOIN_ID       = norm_id(META_SAMPLEID)
    ) %>%
    dplyr::select(-dplyr::any_of("SAMPLEID"))
  
  merged <- avg %>%
    dplyr::left_join(meta, by = "JOIN_ID")
  
  if (!"PATIENTCODE" %in% names(merged)) {
    merged <- merged %>% dplyr::mutate(PATIENTCODE = SAMPLEID_BASE)
  } else {
    merged <- merged %>%
      dplyr::mutate(PATIENTCODE = dplyr::coalesce(as.character(PATIENTCODE), SAMPLEID_BASE))
  }
  
  merged <- merged %>%
    dplyr::mutate(
      CELLTYPE    = factor(CELLTYPE,    levels = CELLTYPE_LEVELS),
      HORMONE     = factor(HORMONE,     levels = HORMONE_LEVELS),
      TIMEPOINT   = factor(TIMEPOINT,   levels = TIMEPOINT_LEVELS),
      SEX         = factor(SEX,         levels = SEX_LEVELS),
      EXPOSURE    = standardize_exposure(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  
  merged <- merged %>%
    dplyr::select(-dplyr::any_of(c("SAMPLEID_BASE", "JOIN_ID", "META_SAMPLEID")))
  
  merged
}

apply_factor_spec <- function(data, 
                              celltype_levels = CELLTYPE_LEVELS,
                              hormone_levels = HORMONE_LEVELS,
                              sex_levels = SEX_LEVELS) {
  # NOTE: This is an alias for the factor-coercion already done inside
  # build_sala_full(). Prefer calling build_sala_full() directly.
  data %>%
    dplyr::mutate(
      CELLTYPE = factor(CELLTYPE, levels = celltype_levels),
      HORMONE  = factor(HORMONE,  levels = hormone_levels),
      SEX      = factor(SEX,      levels = sex_levels)
    )
}

# ============================================================================
# STANDARDIZATION FUNCTIONS
# ============================================================================

#' Standardize sample IDs (remove hyphens, fix prefixes)
standardize_sample_id <- function(x) {
  project_label <- if (exists("WSTC_LABEL", inherits = TRUE)) WSTC_LABEL else "WSTC"
  is_missing <- is.na(x)
  x <- trimws(as.character(x))
  x <- gsub("\\.fastq.*$", "", x, ignore.case = TRUE)
  x <- gsub("\\.bam$", "", x, ignore.case = TRUE)
  x <- gsub("^X", "", x)
  x <- gsub("^(SALA|WSTC)[-_]?", project_label, x, ignore.case = TRUE)
  x <- gsub("[[:space:]]+", "", x)
  x[is_missing] <- NA_character_
  x
}

#' Standardize exposure names (remove E2 suffix, clean formatting)
standardize_exposure <- function(x) {
  is_missing <- is.na(x)
  x <- trimws(as.character(x))
  
  # Remove hormone suffixes first
  x <- gsub(",?\\s*1\\s*nM\\s*B-Estradiol", "", x, ignore.case = TRUE)
  x <- gsub(",?\\s*1\\s*nM\\s*Estradiol", "", x, ignore.case = TRUE)
  x <- gsub(",?\\s*\\+\\s*E2$", "", x, ignore.case = TRUE)
  x <- gsub("\\+\\s*Paraff?in", "", x, ignore.case = TRUE)
  x <- trimws(x)
  
  # Canonical key: lowercase + strip non-alphanumeric
  key <- tolower(gsub("[^a-z0-9]", "", x))
  
  out <- dplyr::case_when(
    key %in% c("pbs", "pbscontrol") ~ "PBS_Control",
    key %in% c("untreated", "untreatedcontrol") ~ "Untreated_Control",
    
    key %in% c("pine5") ~ "Pine_5",
    key %in% c("pine25", "pine") ~ "Pine_25",
    
    key %in% c("peat5") ~ "Peat_5",
    key %in% c("peat25", "peat") ~ "Peat_25",
    
    key %in% c("eucalyptus5") ~ "Eucalyptus_5",
    key %in% c("eucalyptus25", "eucalyptus") ~ "Eucalyptus_25",
    
    key %in% c("redoak5") ~ "RedOak_5",
    key %in% c("redoak25", "redoak", "redoakwood", "redoaksmoke") ~ "RedOak_25",
    
    TRUE ~ x
  )
  
  out <- gsub("[[:space:]-]+", "_", out)
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  out[is_missing] <- NA_character_
  out
}

exposure_display_name <- function(x, include_dose = TRUE) {
  canonical <- standardize_exposure(x)
  display <- dplyr::case_when(
    canonical == "PBS_Control" ~ "PBS Control",
    canonical == "Untreated_Control" ~ "Untreated Control",
    canonical == "Pine_5" ~ "Pine 5",
    canonical == "Pine_25" ~ "Pine 25",
    canonical == "Peat_5" ~ "Peat 5",
    canonical == "Peat_25" ~ "Peat 25",
    canonical == "Eucalyptus_5" ~ "Eucalyptus 5",
    canonical == "Eucalyptus_25" ~ "Eucalyptus 25",
    canonical == "RedOak_5" ~ "Red Oak 5",
    canonical == "RedOak_25" ~ "Red Oak 25",
    TRUE ~ gsub("_", " ", canonical)
  )
  
  if (!include_dose) {
    display <- sub("\\s+[0-9]+$", "", display)
  }
  
  display[is.na(canonical)] <- NA_character_
  display
}

normalize_join_id <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  x <- sub("(_144|_4)$", "", x)   # remove time suffix
  x <- gsub("[^A-Z0-9]", "", x)   # remove underscores, hyphens, spaces, etc.
  x
}

# ============================================================================
# 1. NORMALIZE COUNTS
# ============================================================================

normalize_counts <- function(count_matrix) {
  count_matrix <- as.matrix(count_matrix)
  
  if (!is.numeric(count_matrix)) {
    stop("normalize_counts(): count_matrix must be numeric.")
  }
  
  size_factors <- colSums(count_matrix) / mean(colSums(count_matrix))
  normalized <- sweep(count_matrix, 2, size_factors, "/")
  
  return(normalized)
}

# ============================================================================
# 2. FILTER GENES
# ============================================================================

filter_genes <- function(data, background_threshold = GENE_BACKGROUND_THRESHOLD) {
  if (!"Geneid" %in% names(data)) {
    stop("filter_genes(): data must contain a 'Geneid' column.")
  }
  
  counts <- data %>%
    dplyr::select(-Geneid) %>%
    as.matrix()
  
  if (!is.numeric(counts)) {
    stop("filter_genes(): all columns except 'Geneid' must be numeric counts.")
  }
  
  median_counts  <- apply(counts, 1, median, na.rm = TRUE)
  overall_median <- median(median_counts, na.rm = TRUE)
  
  # Guard: if overall_median is 0 (sparse data), fall back to a count > 0 filter
  # to avoid the threshold collapsing to "> 0 in ≥20% of samples" silently.
  if (!is.finite(overall_median) || overall_median == 0) {
    warning("filter_genes(): overall median is 0 or non-finite; ",
            "applying a simple presence filter (count > 0 in ≥20% of samples).")
    keep_genes <- apply(counts > 0, 1, function(x) sum(x, na.rm = TRUE) > ncol(counts) * 0.20)
  } else {
    keep_genes <- apply(
      counts > (overall_median * background_threshold),
      1,
      function(x) sum(x, na.rm = TRUE) > ncol(counts) * 0.20
    )
  }
  
  filtered_data <- data[keep_genes, , drop = FALSE]
  return(filtered_data)
}

# ============================================================================
# 3. PREPARE METADATA
# ============================================================================

prepare_metadata_seq <- function(metadata_full, sample_ids) {
  id_cols <- c("SAMPLEID", "Sample_ID", "SampleID")
  id_col <- id_cols[id_cols %in% names(metadata_full)][1]
  
  if (is.na(id_col)) {
    stop("prepare_metadata_seq(): no sample ID column found.")
  }
  
  metadata_full[[id_col]] <- standardize_sample_id(metadata_full[[id_col]])
  sample_ids <- standardize_sample_id(sample_ids)
  
  metadata_seq <- metadata_full %>%
    dplyr::filter(.data[[id_col]] %in% sample_ids) %>%
    dplyr::mutate(
      dplyr::across(where(is.character), as.character)
    )
  
  if ("EXPOSURE" %in% names(metadata_seq)) {
    metadata_seq <- metadata_seq %>%
      dplyr::mutate(EXPOSURE = standardize_exposure(EXPOSURE))
  }
  
  metadata_seq
}

# ============================================================================
# 4. CLASSIFY DEGS
# ============================================================================

classify_degs <- function(data, fc_cutoff = LOG2FC_CUTOFF, p_cutoff = ADJ_P_CUTOFF) {
  data %>%
    dplyr::mutate(DEG = dplyr::case_when(
      adj.P.Val <= p_cutoff & log2FC >=  fc_cutoff ~ "UP",
      adj.P.Val <= p_cutoff & log2FC <= -fc_cutoff ~ "DOWN",
      TRUE ~ "NO"
    ))
}

# ============================================================================
# 5. READ GSEA HTML
# ============================================================================

read_gsea_html <- function(html_file) {
  if (!requireNamespace("rvest", quietly = TRUE)) {
    stop("read_gsea_html(): package 'rvest' is required.")
  }
  
  page <- rvest::read_html(html_file)
  tables <- rvest::html_table(page)
  
  if (length(tables) > 0) {
    results <- tables[[1]]
  } else {
    results <- data.frame()
  }
  
  return(results)
}

# ============================================================================
# 6. FILTER GSEA RESULTS
# ============================================================================

filter_gsea_results <- function(gsea_results, fdr_cutoff = FDR_THRESHOLD) {
  gsea_results %>%
    dplyr::mutate(Significant = ifelse(FDR < fdr_cutoff, "Yes", "No")) %>%
    dplyr::arrange(FDR)
}

# ============================================================================
# 7. PREPARE DEG EXPORT
# ============================================================================

prepare_deg_export <- function(deg_data) {
  deg_data %>%
    dplyr::select(Gene = gene_name, log2FC, adj.P.Val) %>%
    dplyr::mutate(
      FC = 2^log2FC,
      log10_padj = -log10(adj.P.Val)
    ) %>%
    dplyr::filter(adj.P.Val < ADJ_P_CUTOFF) %>%
    dplyr::arrange(dplyr::desc(abs(log2FC)))
}

cat("✓ Data / RNA-seq functions loaded\n")
