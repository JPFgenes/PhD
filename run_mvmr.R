#!/usr/bin/env Rscript

# Multivariable Mendelian randomisation example
# Frailty phenotypes adjusted for educational attainment, income, and longevity
# Alzheimer's disease GWAS outcomes
#
# This repository is an example analysis pipeline. Summary-statistics files are
# not included. Update config/analysis_config.R before running.

required_packages <- c("data.table", "dplyr", "TwoSampleMR", "MVMR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(MVMR)
})

source("config/analysis_config.R")
source("R/mvmr_functions.R")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

validate_config(
  exposure_files = exposure_files,
  covariate_files = covariate_files,
  outcome_files = outcome_files,
  plink_bin = plink_bin,
  ld_reference_prefix = ld_reference_prefix
)

message("Reading exposure GWAS files...")
primary_exposures <- lapply(names(exposure_files), function(label) {
  read_standard_exposure(exposure_files[[label]], label)
})
names(primary_exposures) <- names(exposure_files)

covariates <- lapply(names(covariate_files), function(label) {
  read_standard_exposure(covariate_files[[label]], label)
})
names(covariates) <- names(covariate_files)

all_results <- list()
all_strength <- list()
all_qc <- list()

for (primary_label in names(primary_exposures)) {
  message("\nPrimary exposure: ", primary_label)

  exposure_list <- c(
    setNames(list(primary_exposures[[primary_label]]), primary_label),
    covariates
  )

  instruments <- select_union_instruments(
    exposure_list = exposure_list,
    p_threshold = p_threshold,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    plink_bin = plink_bin,
    ld_reference_prefix = ld_reference_prefix
  )

  exposure_wide <- build_exposure_matrices(
    exposure_list = exposure_list,
    instrument_snps = instruments$SNP,
    reference_exposure = primary_label
  )

  for (outcome_label in names(outcome_files)) {
    message("  Outcome: ", outcome_label)

    outcome <- read_outcome(
      path = outcome_files[[outcome_label]]$path,
      format = outcome_files[[outcome_label]]$format,
      column_map = outcome_files[[outcome_label]]$columns,
      label = outcome_label
    )

    aligned <- align_outcome_to_reference(
      exposure_wide = exposure_wide,
      outcome = outcome,
      palindromic_maf_threshold = palindromic_maf_threshold
    )

    analysis_data <- prepare_mvmr_input(aligned$data, names(exposure_list))

    if (nrow(analysis_data$data) <= length(exposure_list)) {
      warning(
        "Skipping ", primary_label, " -> ", outcome_label,
        ": too few complete instruments after harmonisation."
      )
      next
    }

    fit <- MVMR::ivw_mvmr(analysis_data$formatted)
    strength <- MVMR::strength_mvmr(analysis_data$formatted)

    fit_df <- as.data.frame(fit)
    fit_df$exposure <- rownames(fit_df)
    rownames(fit_df) <- NULL
    fit_df$primary_exposure <- primary_label
    fit_df$outcome <- outcome_label
    fit_df$n_instruments <- nrow(analysis_data$data)

    strength_df <- as.data.frame(strength)
    strength_df$exposure <- rownames(strength_df)
    rownames(strength_df) <- NULL
    strength_df$primary_exposure <- primary_label
    strength_df$outcome <- outcome_label

    qc_df <- data.frame(
      primary_exposure = primary_label,
      outcome = outcome_label,
      n_union_instruments = length(instruments$SNP),
      n_after_outcome_alignment = nrow(aligned$data),
      n_complete_for_mvmr = nrow(analysis_data$data),
      n_flipped_outcome = aligned$n_flipped,
      n_removed_allele_mismatch = aligned$n_mismatch,
      n_removed_palindromic = aligned$n_palindromic,
      stringsAsFactors = FALSE
    )

    all_results[[paste(primary_label, outcome_label, sep = "__")]] <- fit_df
    all_strength[[paste(primary_label, outcome_label, sep = "__")]] <- strength_df
    all_qc[[paste(primary_label, outcome_label, sep = "__")]] <- qc_df

    fwrite(
      fit_df,
      file.path(output_dir, paste0(primary_label, "_", outcome_label, "_mvmr.tsv")),
      sep = "\t"
    )
    fwrite(
      strength_df,
      file.path(output_dir, paste0(primary_label, "_", outcome_label, "_strength.tsv")),
      sep = "\t"
    )
  }
}

if (length(all_results) == 0L) {
  stop("No MVMR models were completed. Review the warnings and input files.", call. = FALSE)
}

results_master <- bind_rows(all_results)
strength_master <- bind_rows(all_strength)
qc_master <- bind_rows(all_qc)

# Obtain the p-value column returned by the installed MVMR version.
p_column <- intersect(c("P_value", "p.value", "Pr(>|t|)"), names(results_master))
if (length(p_column) != 1L) {
  warning("Could not identify a unique p-value column; FDR values were not calculated.")
  results_master$q_value <- NA_real_
} else {
  results_master <- results_master %>%
    group_by(exposure) %>%
    mutate(q_value = p.adjust(.data[[p_column]], method = "BH")) %>%
    ungroup()
}

fwrite(results_master, file.path(output_dir, "MVMR_RESULTS_MASTER.tsv"), sep = "\t")
fwrite(strength_master, file.path(output_dir, "MVMR_STRENGTH_MASTER.tsv"), sep = "\t")
fwrite(qc_master, file.path(output_dir, "MVMR_QC_MASTER.tsv"), sep = "\t")

message("\nAnalysis complete. Results written to: ", normalizePath(output_dir))
