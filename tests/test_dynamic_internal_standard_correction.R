suppressPackageStartupMessages(library(tidyverse))

script_path <- normalizePath(file.path(
  dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])),
  "..", "workflow", "scripts", "correct-w-isds-snakemake-v5.R"
))

setClass(
  "MockSnakemake",
  slots = c(input = "list", output = "list", params = "list", config = "list")
)

run_case <- function(number_of_standards) {
  root <- tempfile(paste0("isd-correction-", number_of_standards, "-"))
  dir.create(root, recursive = TRUE)
  old_working_directory <- getwd()
  on.exit({
    setwd(old_working_directory)
    unlink(root, recursive = TRUE)
  }, add = TRUE)
  setwd(root)

  standard_ids <- paste0("Standard-", seq_len(number_of_standards))
  standard_slots <- paste0("isd_", seq_len(number_of_standards))
  standard_hashes <- paste0("internal-hash-", seq_len(number_of_standards))

  asv_table <- bind_rows(
    tibble(
      SampleID = "Sample-1",
      ASV_hash = standard_hashes,
      Corrected_Sequence_Counts = seq_len(number_of_standards) * 100,
      Domain = "Bacteria"
    ),
    tibble(
      SampleID = "Sample-1",
      ASV_hash = "biological-hash",
      Corrected_Sequence_Counts = 500,
      Domain = "Bacteria"
    )
  )
  write_tsv(asv_table, "long_data.tsv")

  internal_stds <- tibble(
    internal_std_ID = standard_ids,
    rRNA_copy_number = 2,
    genome_len_bp = 1000,
    full_16S_sequence = "ACGT"
  )
  write_tsv(internal_stds, "internal_stds.tsv")

  sample_row <- tibble(
    sample = "Sample_1",
    internal_std_normalization_factor = 1,
    units = "L"
  )
  for (standard_id in standard_ids) {
    sample_row[[paste0(standard_id, "_ng")]] <- 1
  }
  write_tsv(sample_row, "samples.tsv")

  asv_paths <- paste0(standard_ids, ".asvs.txt")
  walk2(asv_paths, standard_hashes, ~writeLines(.y, .x))

  pair_count <- choose(number_of_standards, 2)
  method_count <- number_of_standards + 2 + pair_count
  method_outputs <- file.path(root, paste0("method-", seq_len(method_count), ".tsv"))
  id_pairs <- if (number_of_standards >= 2) {
    combn(standard_ids, 2, simplify = FALSE)
  } else {
    list()
  }
  all_ids_stem <- paste(standard_ids, collapse = "_")
  method_stems <- c(
    paste0(standard_ids, "_recovery_ratio"),
    paste0("mean_", all_ids_stem, "_recovery_ratio"),
    paste0("median_", all_ids_stem, "_recovery_ratio"),
    vapply(
      id_pairs,
      function(pair) paste0("mean_", pair[[1]], "_and_", pair[[2]], "_recovery_ratio"),
      character(1)
    )
  )

  snakemake <<- new(
    "MockSnakemake",
    input = list(
      asv_table = file.path(root, "long_data.tsv"),
      standard_asvs = file.path(root, asv_paths),
      isd = file.path(root, "internal_stds.tsv"),
      isd_added = file.path(root, "samples.tsv")
    ),
    output = list(
      corrected = file.path(root, "corrected.tsv"),
      filtered = file.path(root, "isd-removed.tsv"),
      method_tables = method_outputs,
      recovery_plot = file.path(root, "recovery.pdf"),
      domain_plot = file.path(root, "domain.pdf"),
      recovery_plot_png = file.path(root, "recovery.png"),
      domain_plot_png = file.path(root, "domain.png")
    ),
    params = list(
      standard_ids = standard_ids,
      standard_slots = standard_slots,
      method_stems = method_stems
    ),
    config = list(projectName = "TEST")
  )

  suppressMessages(source(script_path, local = globalenv()))

  stopifnot(file.exists(snakemake@output[["corrected"]]))
  stopifnot(file.exists(snakemake@output[["filtered"]]))
  stopifnot(all(file.exists(method_outputs)))
  stopifnot(file.exists(snakemake@output[["recovery_plot_png"]]))
  stopifnot(file.exists(snakemake@output[["domain_plot_png"]]))

  corrected <- read_tsv(snakemake@output[["corrected"]], show_col_types = FALSE)
  filtered <- read_tsv(snakemake@output[["filtered"]], show_col_types = FALSE)
  first_wide <- read_tsv(method_outputs[[1]], show_col_types = FALSE)
  stopifnot(!"internal-hash-1" %in% corrected$ASV_hash)
  stopifnot("biological-hash" %in% corrected$ASV_hash)
  stopifnot(!any(standard_hashes %in% filtered$ASV_hash))
  stopifnot("biological-hash" %in% filtered$ASV_hash)
  stopifnot(identical(names(filtered), names(asv_table)))
  stopifnot(paste0("Copies_Standard-", number_of_standards, "_recovery_ratio") %in% names(corrected))
  stopifnot("Domain" %in% names(first_wide))
  stopifnot("Sample-1" %in% names(first_wide))
  if (number_of_standards >= 2) {
    stopifnot(any(grepl("_and_", names(corrected), fixed = TRUE)))
  }
}

run_case(1)
run_case(3)
run_case(4)
cat("Dynamic internal-standard correction tests passed\n")
