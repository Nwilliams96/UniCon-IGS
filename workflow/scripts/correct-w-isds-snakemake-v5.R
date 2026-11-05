# Convert ASV sequence counts into estimated absolute copies using any positive
# number of genomic internal standards added at DNA extraction.


suppressPackageStartupMessages({
  library(tidyverse)
})

standard_ids <- as.character(unlist(
  snakemake@params[["standard_ids"]], use.names = FALSE
))
standard_slots <- as.character(unlist(
  snakemake@params[["standard_slots"]], use.names = FALSE
))
standard_asv_paths <- as.character(unlist(
  snakemake@input[["standard_asvs"]], use.names = FALSE
))
method_table_outputs <- as.character(unlist(
  snakemake@output[["method_tables"]], use.names = FALSE
))
method_stems <- as.character(unlist(
  snakemake@params[["method_stems"]], use.names = FALSE
))

if (length(standard_ids) < 1) {
  stop("At least one internal standard is required for correction")
}
if (length(unique(standard_ids)) != length(standard_ids)) {
  stop("Configured internal-standard IDs must be unique")
}
if (length(standard_slots) != length(standard_ids) ||
    length(standard_asv_paths) != length(standard_ids)) {
  stop("Internal-standard IDs, slots, and ASV-list inputs are not aligned")
}

configured_standard_names <- paste(standard_ids, collapse = ", ")
cruise_tag <- snakemake@config[["projectName"]]
outdir <- "results/05-internal-std-corrected"
figure_dir <- "results/06-figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_asv_ids <- function(path) {
  values <- trimws(readLines(path, warn = FALSE))
  unique(values[nzchar(values)])
}

row_summary <- function(values, summary_function) {
  apply(as.matrix(values), 1, function(row) {
    if (all(is.na(row))) NA_real_ else summary_function(row, na.rm = TRUE)
  })
}

make_wide_table <- function(data, copy_column, annotation_columns) {
  annotations <- data %>%
    select(ASV_hash, all_of(annotation_columns)) %>%
    distinct(ASV_hash, .keep_all = TRUE)
  wide <- data %>%
    select(SampleID, ASV_hash, all_of(copy_column)) %>%
    group_by(ASV_hash, SampleID) %>%
    summarise(Abundance = sum(.data[[copy_column]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)
  sample_columns <- sort(setdiff(names(wide), "ASV_hash"))
  annotations %>%
    left_join(wide, by = "ASV_hash") %>%
    select(ASV_hash, all_of(annotation_columns), all_of(sample_columns)) %>%
    arrange(ASV_hash)
}

# Input data and validation -------------------------------------------------
asv_table <- read_tsv(snakemake@input[["asv_table"]], show_col_types = FALSE)
formatted_table_columns <- names(asv_table)
# Preserve every ASV-level annotation supplied by the formatted table. This is
# intentionally data-driven so taxonomy columns from future classifiers are
# carried into every wide correction table without another script change.
measurement_columns <- c(
  "SampleID", "Raw_Sequence_Counts", "Corrected_dada2_Sequence_Counts",
  "Corrected_Sequence_Counts", "Relative_Abundance"
)
wide_annotation_columns <- setdiff(names(asv_table), c(measurement_columns, "ASV_hash"))
isd <- read_tsv(snakemake@input[["isd"]], show_col_types = FALSE)
isd_added <- read_tsv(snakemake@input[["isd_added"]], show_col_types = FALSE) %>%
  rename(SampleID = sample) %>%
  mutate(SampleID = str_replace_all(SampleID, "_", "-"))

required_isd_columns <- c(
  "internal_std_ID", "rRNA_copy_number", "genome_len_bp"
)
missing_isd_columns <- setdiff(required_isd_columns, names(isd))
if (length(missing_isd_columns) > 0) {
  stop(
    "config/internal_stds.tsv is missing column(s): ",
    paste(missing_isd_columns, collapse = ", ")
  )
}

amount_columns <- paste0(standard_ids, "_ng")
required_sample_columns <- c(
  "SampleID", amount_columns, "internal_std_normalization_factor", "units"
)
missing_sample_columns <- setdiff(required_sample_columns, names(isd_added))
if (length(missing_sample_columns) > 0) {
  stop(
    "config/samples.tsv is missing column(s): ",
    paste(missing_sample_columns, collapse = ", ")
  )
}

isd_lookup <- isd %>%
  transmute(
    StandardID = as.character(internal_std_ID),
    rRNA_copy_number = as.numeric(rRNA_copy_number),
    genome_len_bp = as.numeric(genome_len_bp)
  )
if (anyDuplicated(isd_lookup$StandardID)) {
  stop("config/internal_stds.tsv contains duplicate internal_std_ID values")
}
missing_standard_rows <- setdiff(standard_ids, isd_lookup$StandardID)
if (length(missing_standard_rows) > 0) {
  stop(
    "Configured internal standard(s) missing from config/internal_stds.tsv: ",
    paste(missing_standard_rows, collapse = ", ")
  )
}
if (any(
  is.na(isd_lookup$rRNA_copy_number) | isd_lookup$rRNA_copy_number <= 0 |
    is.na(isd_lookup$genome_len_bp) | isd_lookup$genome_len_bp <= 0
)) {
  stop("Internal-standard copy numbers and genome lengths must be positive numbers")
}

standard_lookup <- tibble(
  StandardID = standard_ids,
  StandardSlot = standard_slots,
  amount_column = amount_columns
)

# Copies added and recovered -----------------------------------------------
bp_weight <- 617.9
avogadro <- 6.022e23

copies_added <- isd_added %>%
  select(SampleID, all_of(amount_columns)) %>%
  pivot_longer(
    cols = all_of(amount_columns),
    names_to = "amount_column",
    values_to = "amount_ng"
  ) %>%
  left_join(standard_lookup, by = "amount_column") %>%
  left_join(isd_lookup, by = "StandardID") %>%
  mutate(
    amount_ng = as.numeric(amount_ng),
    copies_added = (((amount_ng / 1e9) /
      (bp_weight * genome_len_bp)) * avogadro) * rRNA_copy_number
  )

recovered_by_standard <- map2_dfr(
  standard_slots,
  standard_asv_paths,
  function(slot, path) {
    standard_asvs <- read_asv_ids(path)
    recovered <- asv_table %>%
      filter(ASV_hash %in% standard_asvs) %>%
      group_by(SampleID) %>%
      summarise(copies_recovered = sum(Corrected_Sequence_Counts, na.rm = TRUE),
                .groups = "drop")

    isd_added %>%
      distinct(SampleID) %>%
      left_join(recovered, by = "SampleID") %>%
      mutate(
        StandardSlot = slot,
        copies_recovered = replace_na(copies_recovered, 0)
      )
  }
)

recovery_long <- copies_added %>%
  left_join(recovered_by_standard, by = c("SampleID", "StandardSlot")) %>%
  mutate(
    copies_recovered = replace_na(copies_recovered, 0),
    recovery_ratio = if_else(
      !is.na(copies_added) & copies_added > 0,
      copies_recovered / copies_added,
      NA_real_
    )
  )

recovery_wide <- recovery_long %>%
  select(SampleID, StandardSlot, recovery_ratio) %>%
  pivot_wider(
    names_from = StandardSlot,
    values_from = recovery_ratio,
    names_glue = "{StandardSlot}_recovery_ratio"
  )

standard_ratio_columns <- paste0(standard_slots, "_recovery_ratio")
recovery_wide$recovery_mean <- row_summary(
  recovery_wide[standard_ratio_columns], mean
)
recovery_wide$recovery_median <- row_summary(
  recovery_wide[standard_ratio_columns], median
)

pair_indexes <- if (length(standard_slots) >= 2) {
  combn(seq_along(standard_slots), 2, simplify = FALSE)
} else {
  list()
}
pair_ratio_columns <- character()
pair_labels <- character()
for (pair in pair_indexes) {
  pair_slots <- standard_slots[pair]
  pair_column <- paste0(
    paste(pair_slots, collapse = "_"), "_mean_recovery_ratio"
  )
  recovery_wide[[pair_column]] <- row_summary(
    recovery_wide[paste0(pair_slots, "_recovery_ratio")], mean
  )
  pair_ratio_columns <- c(pair_ratio_columns, pair_column)
  pair_labels <- c(pair_labels, paste(standard_ids[pair], collapse = " + "))
}

# Define all correction methods in one ordered table. This order matches the
# output list assembled in the Snakemake rule.
method_specs <- bind_rows(
  tibble(
    stem = paste0(standard_slots, "_recovery_ratio"),
    ratio_column = standard_ratio_columns,
    label = paste0(standard_ids, " only")
  ),
  tibble(
    stem = c("mean_recovery_ratio", "median_recovery_ratio"),
    ratio_column = c("recovery_mean", "recovery_median"),
    label = c(
      paste0("Mean: ", configured_standard_names),
      paste0("Median: ", configured_standard_names)
    )
  ),
  tibble(
    stem = pair_ratio_columns,
    ratio_column = pair_ratio_columns,
    label = pair_labels
  )
) %>%
  mutate(stem = method_stems) %>%
  mutate(copy_column = paste0("Copies_", stem))

if (length(method_table_outputs) != nrow(method_specs)) {
  stop("The number of correction-table outputs does not match the methods")
}

# Recovery-ratio figure -----------------------------------------------------
plot_data <- recovery_wide %>%
  select(SampleID, all_of(method_specs$ratio_column)) %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "ratio_column",
    values_to = "Recovery"
  ) %>%
  left_join(method_specs %>% select(ratio_column, Method = label),
            by = "ratio_column") %>%
  filter(!is.na(Recovery))

means_df <- plot_data %>%
  filter(ratio_column == "recovery_mean") %>%
  distinct(SampleID, Recovery)

colourblind_palette <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
  "#DDCC77", "#CC6677", "#882255", "#AA4499", "#0072B2",
  "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#64748B"
)
method_colours <- setNames(
  rep(colourblind_palette, length.out = nrow(method_specs)),
  method_specs$label
)

plot1 <- ggplot(plot_data, aes(x = SampleID, y = Recovery, colour = Method)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2.5) +
  geom_line(
    data = means_df,
    aes(x = SampleID, y = Recovery, group = 1),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1
  ) +
  theme_minimal() +
  scale_y_log10(expand = expansion(mult = c(0.15, 0.05))) +
  scale_colour_manual(values = method_colours, name = "Correction method") +
  labs(
    title = "Recovery Ratios per Sample",
    subtitle = paste0(
      "Configured standards: ", configured_standard_names,
      ". Black line = per-sample mean."
    ),
    y = "Recovery Ratio (log10 scale)",
    x = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

# Apply each correction method ---------------------------------------------
asv_table <- asv_table %>%
  left_join(recovery_wide, by = "SampleID") %>%
  left_join(
    isd_added %>%
      select(SampleID, internal_std_normalization_factor, units),
    by = "SampleID"
  )

for (index in seq_len(nrow(method_specs))) {
  ratio_column <- method_specs$ratio_column[[index]]
  copy_column <- method_specs$copy_column[[index]]
  asv_table[[copy_column]] <-
    (asv_table$Corrected_Sequence_Counts / asv_table[[ratio_column]]) /
    asv_table$internal_std_normalization_factor
}

all_internal_standard_asvs <- unique(unlist(
  map(standard_asv_paths, read_asv_ids), use.names = FALSE
))
isd_filtered_asv_table <- asv_table %>%
  filter(!ASV_hash %in% all_internal_standard_asvs) %>%
  select(all_of(formatted_table_columns))
asv_table <- asv_table %>%
  filter(!ASV_hash %in% all_internal_standard_asvs)

# Domain totals figure ------------------------------------------------------
domain_totals_long <- asv_table %>%
  filter(!Domain %in% "Unassigned") %>%
  group_by(SampleID, Domain) %>%
  summarise(
    across(all_of(method_specs$copy_column), ~sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = all_of(method_specs$copy_column),
    names_to = "copy_column",
    values_to = "Total_Copies"
  ) %>%
  left_join(method_specs %>% select(copy_column, Method = label),
            by = "copy_column")

plot2 <- ggplot(
  domain_totals_long,
  aes(x = SampleID, y = Total_Copies, colour = Domain, group = Domain)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10(expand = expansion(mult = c(0.15, 0.05))) +
  facet_wrap(~Method, ncol = 1) +
  theme_minimal() +
  scale_colour_manual(values = c(
    "Bacteria" = "#0072B2",
    "Archaea" = "#E69F00",
    "Eukaryota" = "#009E73",
    "Unassigned" = "#64748B"
  )) +
  labs(
    title = "Total Copies per Unit by Domain and Correction Method",
    y = "Total Copies per unit (log10 scale)",
    x = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom"
  )

# Tables --------------------------------------------------------------------
write_tsv(isd_filtered_asv_table, snakemake@output[["filtered"]])
write_tsv(asv_table, snakemake@output[["corrected"]])
walk2(
  method_table_outputs,
  method_specs$copy_column,
  function(output_path, copy_column) {
    write_tsv(
      make_wide_table(asv_table, copy_column, wide_annotation_columns),
      output_path
    )
  }
)

# Figures -------------------------------------------------------------------
domain_plot_height <- max(8, 2.4 * nrow(method_specs))

pdf(snakemake@output[["recovery_plot"]], width = 12, height = 6)
print(plot1)
dev.off()
pdf(snakemake@output[["domain_plot"]], width = 12, height = domain_plot_height)
print(plot2)
dev.off()

sample_plot_width <- max(1800, 42 * n_distinct(asv_table$SampleID))
png(
  snakemake@output[["recovery_plot_png"]],
  width = sample_plot_width, height = 900, res = 150
)
print(plot1)
dev.off()
png(
  snakemake@output[["domain_plot_png"]],
  width = sample_plot_width,
  height = max(1200, 360 * nrow(method_specs)),
  res = 150
)
print(plot2)
dev.off()
