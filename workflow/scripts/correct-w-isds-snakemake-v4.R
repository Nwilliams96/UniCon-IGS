#This script was written to convert ASV copies into ASV absolute copies by correcting data with internal standards.
#It can only be used if genomic internal standards were added at the time of DNA extraction.
#Written by Nathan Williams 18/02/2026.

suppressPackageStartupMessages({
library(tidyverse)
})

#Set paths
isd_path  <- snakemake@input[["isd"]]
isd_added_path <- snakemake@input[["isd_added"]]
cruise_tag <- snakemake@config[["projectName"]]
isd_1_id <- snakemake@params[["intstd1name"]]
isd_2_id <- snakemake@params[["intstd2name"]]
isd_3_id <- snakemake@params[["intstd3name"]]
configured_standard_names <- paste(c(isd_1_id, isd_2_id, isd_3_id), collapse = ", ")
outdir <- "results/05-internal-std-corrected"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- "results/06-figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_isd_table <- function(df, stem) {
  readr::write_tsv(df, file.path(outdir, paste0(cruise_tag, ".", stem, ".tsv")))
}

#Import Data
asv_table <- read_tsv(snakemake@input[["asv_table"]])
isd <- read_tsv(isd_path , show_col_types = FALSE)
isd_added <- read_tsv(isd_added_path, show_col_types = FALSE) %>%
  rename(
    SampleID = sample,
    isd_1_ng = all_of(paste0(isd_1_id, "_ng")),
    isd_2_ng = all_of(paste0(isd_2_id, "_ng")),
    isd_3_ng = all_of(paste0(isd_3_id, "_ng"))
  ) %>%
  mutate(SampleID = str_replace_all(SampleID, "_", "-"))
samples <- isd_added

#Import Data local
isd_1_asvs <- read_delim(snakemake@input[[2]], delim = "\n", col_names = FALSE)
isd_2_asvs <- read_delim(snakemake@input[[3]], delim = "\n", col_names = FALSE)
isd_3_asvs <- read_delim(snakemake@input[[4]], delim = "\n", col_names = FALSE)

#Make the ISD dataframe lookup vectors
genome_len  <- setNames(isd$genome_len_bp, isd$internal_std_ID)
copy_number <- setNames(isd$rRNA_copy_number, isd$internal_std_ID)

#Calculate copies of each ISD added
#Set parameters
bp_weight <- 617.9
avogadro <- 6.022 * 1e23
# 1e9 is to convert to copies added per L.

isd_1len=isd$genome_len_bp[isd$internal_std_ID == isd_1_id]
isd_2len=isd$genome_len_bp[isd$internal_std_ID == isd_2_id]
isd_3len=isd$genome_len_bp[isd$internal_std_ID == isd_3_id]

isd_1copynum=isd$rRNA_copy_number[isd$internal_std_ID == isd_1_id]
isd_2copynum=isd$rRNA_copy_number[isd$internal_std_ID == isd_2_id]
isd_3copynum=isd$rRNA_copy_number[isd$internal_std_ID == isd_3_id]

# Do calculation
isd_copies_added <- isd_added %>%
  select(SampleID, isd_3_ng, isd_1_ng, isd_2_ng) %>%
  mutate(isd_3_copies = ((((isd_3_ng/1e9) / (bp_weight * isd_3len)) * avogadro) * isd_3copynum)) %>%
  mutate(isd_2_copies = ((((isd_2_ng/1e9) / (bp_weight * isd_2len)) * avogadro) * isd_2copynum)) %>%
  mutate(isd_1_copies = ((((isd_1_ng/1e9) / (bp_weight * isd_1len)) * avogadro) * isd_1copynum))

#Pull internal standard copies out of ASV table frame
isd_1_ids <- pull(isd_1_asvs) %>% as.character()
isd_2_ids <- pull(isd_2_asvs) %>% as.character()
isd_3_ids <- pull(isd_3_asvs) %>% as.character()

isd_1_by_sample <- asv_table %>%
  filter(ASV_hash %in% isd_1_ids) %>%
  group_by(SampleID) %>%
  summarize(isd_1_copies_recovered = sum(Corrected_Sequence_Counts, na.rm = TRUE), .groups = "drop")

isd_2_by_sample <- asv_table %>%
  filter(ASV_hash %in% isd_2_ids) %>%
  group_by(SampleID) %>%
  summarize(isd_2_copies_recovered = sum(Corrected_Sequence_Counts, na.rm = TRUE), .groups = "drop")

isd_3_by_sample <- asv_table %>%
  filter(ASV_hash %in% isd_3_ids) %>%
  group_by(SampleID) %>%
  summarize(isd_3_copies_recovered = sum(Corrected_Sequence_Counts, na.rm = TRUE), .groups = "drop")

#Now calculate recovery ratio
merged_isd_data <- isd_copies_added %>%
  left_join(isd_1_by_sample) %>%
  left_join(isd_2_by_sample) %>%
  left_join(isd_3_by_sample)

#Calculate per-ISD recovery ratio columns (recovered / added)
merged_isd_data <- merged_isd_data %>%
  mutate(isd_3_recovery_ratio = isd_3_copies_recovered / isd_3_copies) %>%
  mutate(isd_1_recovery_ratio = isd_1_copies_recovered / isd_1_copies) %>%
  mutate(isd_2_recovery_ratio = isd_2_copies_recovered / isd_2_copies)

#Calculate the mean and median of all three combinations, as well as each combination of two.
merged_isd_data <- merged_isd_data %>%
  rowwise() %>%
  mutate(recovery_mean = mean( c(isd_1_recovery_ratio, isd_2_recovery_ratio, isd_3_recovery_ratio), na.rm = TRUE),
  recovery_median = median(c(isd_1_recovery_ratio, isd_2_recovery_ratio, isd_3_recovery_ratio),na.rm = TRUE),
  isd_1_isd_2_mean_recovery_ratio = mean( c(isd_1_recovery_ratio, isd_2_recovery_ratio), na.rm = TRUE),
  isd_1_isd_3_mean_recovery_ratio = mean( c(isd_1_recovery_ratio, isd_3_recovery_ratio), na.rm = TRUE),
  isd_2_isd_3_mean_recovery_ratio = mean( c(isd_2_recovery_ratio, isd_3_recovery_ratio), na.rm = TRUE)
  ) %>%
  ungroup()

#Make a subset of data for the plots
plot_data <- merged_isd_data %>%
  select(SampleID, isd_3_recovery_ratio, isd_1_recovery_ratio, isd_2_recovery_ratio, recovery_mean,
         isd_1_isd_2_mean_recovery_ratio, isd_1_isd_3_mean_recovery_ratio, isd_2_isd_3_mean_recovery_ratio) %>%
  pivot_longer(
    cols = c(isd_3_recovery_ratio, isd_1_recovery_ratio, isd_2_recovery_ratio, recovery_mean,
             isd_1_isd_2_mean_recovery_ratio, isd_1_isd_3_mean_recovery_ratio, isd_2_isd_3_mean_recovery_ratio),
    names_to = "Method",
    values_to = "Recovery"
  ) %>%
  filter(!is.na(Recovery))

#means dataframe
means_df <- plot_data %>% filter(Method %in% c("recovery_mean")) %>%
  distinct(SampleID, Recovery)

plot_data <- plot_data %>%
  mutate(Method = recode(
    Method,
    "isd_1_recovery_ratio" = paste0(isd_1_id, " only"),
    "isd_2_recovery_ratio" = paste0(isd_2_id, " only"),
    "isd_3_recovery_ratio" = paste0(isd_3_id, " only"),
    "recovery_mean" = "Mean of all three standards",
    "isd_1_isd_2_mean_recovery_ratio" = paste(isd_1_id, "+", isd_2_id),
    "isd_1_isd_3_mean_recovery_ratio" = paste(isd_1_id, "+", isd_3_id),
    "isd_2_isd_3_mean_recovery_ratio" = paste(isd_2_id, "+", isd_3_id)
  ))

#Generate a plot to compare ratios
plot1 <- ggplot(plot_data, aes(x = SampleID, y = Recovery, colour = Method)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2.5) +
  geom_line(data = means_df, aes(x = SampleID, y = Recovery, group = 1),
            inherit.aes = FALSE, colour = "black", linewidth = 1) +
  theme_minimal() +
  scale_colour_brewer(palette = "Dark2", name = "Correction method") +
  labs(title = "Recovery Ratios per Sample",
       subtitle = paste0("Configured standards: ", configured_standard_names, ". Black line = per-sample mean."),
       y = "Recovery Ratio",
       x = "Sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
plot1

#Choose your method
recovery_ratios <- merged_isd_data %>% select(c("SampleID", "isd_3_recovery_ratio", "isd_1_recovery_ratio", "isd_2_recovery_ratio" ,
                                                "recovery_mean", "recovery_median", "isd_1_isd_2_mean_recovery_ratio",
                                                "isd_1_isd_3_mean_recovery_ratio", "isd_2_isd_3_mean_recovery_ratio"))

#Join in recovery ratio
asv_table <- asv_table %>%
  left_join(recovery_ratios)


#add in unit for normalisation
isd_norm_fact <- samples %>%
  select(SampleID,internal_std_normalization_factor,units)

asv_table <- asv_table %>% left_join(isd_norm_fact)

#Calculate recovery ratio
asv_table <- asv_table %>%
  mutate(Copies_isd_1_RR = (Corrected_Sequence_Counts / isd_1_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_2_RR = (Corrected_Sequence_Counts / isd_2_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_3_RR = (Corrected_Sequence_Counts / isd_3_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_mean_RR = (Corrected_Sequence_Counts / recovery_mean)/internal_std_normalization_factor) %>%
  mutate(Copies_median_RR = (Corrected_Sequence_Counts / recovery_median)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_1_isd_2_mean_recovery_ratio_RR = (Corrected_Sequence_Counts /isd_1_isd_2_mean_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_1_isd_3_mean_recovery_ratio_RR = (Corrected_Sequence_Counts /isd_1_isd_3_mean_recovery_ratio)/internal_std_normalization_factor) %>%
  mutate(Copies_isd_2_isd_3_mean_recovery_ratio_RR = (Corrected_Sequence_Counts /isd_2_isd_3_mean_recovery_ratio)/internal_std_normalization_factor)

#Remove ISDs from the data
asv_table <- asv_table %>%
  filter(!ASV_hash %in% isd_1_ids) %>%
  filter(!ASV_hash %in% isd_2_ids) %>%
  filter(!ASV_hash %in% isd_3_ids)

Domain_Totals <- asv_table %>%
  filter(!Domain %in% c("Unassigned")) %>%
  group_by(SampleID, Domain) %>%
  summarise(
    Total_isd_3_RR     = sum(Copies_isd_3_RR, na.rm = TRUE),
    Total_isd_1_RR     = sum(Copies_isd_1_RR, na.rm = TRUE),
    Total_isd_2_RR     = sum(Copies_isd_2_RR, na.rm = TRUE),
    Total_mean_RR   = sum(Copies_mean_RR, na.rm = TRUE),
    Total_median_RR = sum(Copies_median_RR, na.rm = TRUE),
    Total_mean_isd_1_isd_3_RR = sum(Copies_isd_1_isd_2_mean_recovery_ratio_RR, na.rm = TRUE),
    Total_mean_isd_1_isd_2_RR = sum(Copies_isd_1_isd_3_mean_recovery_ratio_RR, na.rm = TRUE),
    Total_mean_isd_2_isd_3_RR = sum(Copies_isd_2_isd_3_mean_recovery_ratio_RR, na.rm = TRUE),
    .groups = "drop"
  )

Domain_Totals_long <- Domain_Totals %>%
  pivot_longer(
    cols = starts_with("Total_"),
    names_to = "Method",
    values_to = "Total_Copies"
  ) %>%
  mutate(Method = recode(Method, "Total_isd_3_RR" = paste0(isd_3_id, " only"), "Total_isd_1_RR" = paste0(isd_1_id, " only"), "Total_isd_2_RR" = paste0(isd_2_id, " only"),
                         "Total_mean_RR" = paste0("Mean: ", configured_standard_names), "Total_median_RR" = paste0("Median: ", configured_standard_names), "Total_mean_isd_1_isd_3_RR" = paste(isd_1_id, "+", isd_3_id),
                         "Total_mean_isd_1_isd_2_RR" = paste(isd_1_id, "+", isd_2_id), "Total_mean_isd_2_isd_3_RR" = paste(isd_2_id, "+", isd_3_id)))

plot2 <- ggplot(
  Domain_Totals_long,
  aes(x = SampleID, y = Total_Copies, colour = Domain, group = Domain)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~Method, ncol = 2) +
  theme_minimal() +
  labs(
    title = "Total Copies per unit by Domain and Correction Method",
    y = "Total Copies per unit (log10 scale)",
    x = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
plot2

# isd_3_RR wide table
asv_table_isd_3_RR <- asv_table %>% select(SampleID, ASV_hash, Copies_isd_3_RR) %>%
  group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_3_RR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# isd_1_RR wide table
asv_table_isd_1_RR <- asv_table %>% select(SampleID, ASV_hash, Copies_isd_1_RR) %>%
  group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_1_RR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# isd_2_RR wide table
asv_table_isd_2_RR <- asv_table %>% select(SampleID, ASV_hash, Copies_isd_2_RR) %>%
  group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_2_RR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# mean_RR wide table
asv_table_mean_RR <- asv_table %>%
  select(SampleID, ASV_hash, Copies_mean_RR) %>% group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_mean_RR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# median_RR wide table
asv_table_median_RR <- asv_table %>%
  select(SampleID, ASV_hash, Copies_median_RR) %>%
  group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_median_RR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# mean_isd_1_isd_2_RR wide table
asv_table_mean_isd_1_isd_2_RR <- asv_table %>%
  select(SampleID, ASV_hash, Copies_isd_1_isd_2_mean_recovery_ratio_RR) %>% group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_1_isd_2_mean_recovery_ratio_RR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# mean_isd_1_isd_3_RR wide table
asv_table_mean_isd_1_isd_3_RR <- asv_table %>%
  select(SampleID, ASV_hash, Copies_isd_1_isd_3_mean_recovery_ratio_RR) %>% group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_1_isd_3_mean_recovery_ratio_RR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

# mean_isd_2_isd_3_RR wide table
asv_table_mean_isd_2_isd_3_RR <- asv_table %>%
  select(SampleID, ASV_hash, Copies_isd_2_isd_3_mean_recovery_ratio_RR) %>% group_by(ASV_hash, SampleID) %>%
  summarise(Abundance = sum(Copies_isd_2_isd_3_mean_recovery_ratio_RR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance, values_fill = 0)

#Write tsv files - one for each correction method.
write_isd_table(asv_table, "ISD_corrected_asv_table")
write_isd_table(asv_table_isd_3_RR, "asv_table_isd_3_RR")
write_isd_table(asv_table_isd_1_RR, "asv_table_isd_1_RR")
write_isd_table(asv_table_isd_2_RR, "asv_table_isd_2_RR")
write_isd_table(asv_table_mean_RR, "asv_table_mean_RR")
write_isd_table(asv_table_median_RR, "asv_table_median_RR")
write_isd_table(asv_table_mean_isd_1_isd_2_RR, "asv_table_mean_isd_1_isd_2_RR")
write_isd_table(asv_table_mean_isd_1_isd_3_RR, "asv_table_mean_isd_1_isd_3_RR")
write_isd_table(asv_table_mean_isd_2_isd_3_RR, "asv_table_mean_isd_2_isd_3_RR")

#Write plots
pdf(snakemake@output[["recovery_plot"]], width = 12, height = 6)
print(plot1)
dev.off()
pdf(snakemake@output[["domain_plot"]], width = 12, height = 8)
print(plot2)
dev.off()
