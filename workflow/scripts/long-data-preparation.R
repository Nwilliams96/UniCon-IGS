#!/usr/bin/env Rscript

#Load dependencies
library(tidyverse)
library(data.table)
library("Biostrings")

#Import data
counts <- readr::read_tsv(snakemake@input[["mergedtabledada218Scorrected"]], show_col_types = FALSE) %>% as.data.table()

counts_dada2_corrected <- readr::read_tsv(snakemake@input[["mergedtabledada2"]], show_col_types = FALSE) %>% as.data.table()

counts_uncorrected <- readr::read_tsv(snakemake@input[["mergedtableuncorrected"]], show_col_types = FALSE) %>% as.data.table()

#Import all the ASV sequences from the 16S and 18S data
prokaryote_asv_sequences <- readDNAStringSet(snakemake@input[["fasta16S"]])
ASV_hash <- names(prokaryote_asv_sequences)
ASV <- paste(prokaryote_asv_sequences)
prokaryote_asv_sequences <- data.frame(ASV_hash, ASV)
eukaryote_asv_sequences  <- readDNAStringSet(snakemake@input[["fasta18S"]])
ASV_hash <- names(eukaryote_asv_sequences)
ASV <- paste(eukaryote_asv_sequences)
eukaryote_asv_sequences <- data.frame(ASV_hash, ASV)

#Join together asv sequences
asv_sequences <- bind_rows(eukaryote_asv_sequences,prokaryote_asv_sequences)
asv_sequences <- write_tsv(asv_sequences, snakemake@output[["asvsequences"]])
asv_sequences <- as.data.frame(asv_sequences)

# Extract a named QIIME 2 taxonomy rank without relying on its position in the
# lineage. SILVA 144 adds k__ (Kingdom), while older SILVA releases did not, so
# positional splitting would shift every rank after Domain.
extract_prefixed_rank <- function(taxonomy, prefix) {
  rank_value <- stringr::str_match(
    taxonomy,
    paste0("(?:^|;\\s*)", prefix, "([^;]*)")
  )[, 2]
  dplyr::na_if(trimws(rank_value), "")
}

parse_prefixed_taxonomy <- function(taxonomy_table) {
  taxonomy_table %>%
    mutate(
      Domain = extract_prefixed_rank(Taxonomy, "d__"),
      Kingdom = extract_prefixed_rank(Taxonomy, "k__"),
      Phylum = extract_prefixed_rank(Taxonomy, "p__"),
      Class = extract_prefixed_rank(Taxonomy, "c__"),
      Order = extract_prefixed_rank(Taxonomy, "o__"),
      Family = extract_prefixed_rank(Taxonomy, "f__"),
      Genus = extract_prefixed_rank(Taxonomy, "g__"),
      Species = extract_prefixed_rank(Taxonomy, "s__")
    ) %>%
    select(-Taxonomy)
}

# Parse out plastid labels to record whether a sequence came from a plastid.
# Confident PR2 calls contain :plas.  A SILVA o__Chloroplast label is retained
# when PR2 does not make a confident plastid call, so recognize both forms.
Taxonomy <- counts %>%
  select(Taxonomy, ProPortal_ASV_Ecotype, ASV_hash)
Taxonomy <- Taxonomy %>%
  mutate(
    plastid_16S_rRNA = case_when(
      str_detect(Taxonomy, regex(":plas|(?:^|;\\s*)o__Chloroplast(?:\\s*;|\\s*$)", ignore_case = TRUE)) ~ "yes",
      TRUE ~ "no"
    )
  )
Taxonomy <- Taxonomy %>% 
  mutate(Source_database = case_when(str_detect(Taxonomy, "d__") ~ "SILVA", TRUE ~ "PR2"))

# This is to create a dataframe for proportal assigned taxa - this is requried for the "source_database" column
ProPortal <- Taxonomy %>% 
  filter(!ProPortal_ASV_Ecotype %in% (NA)) %>% mutate(Source_database = c("ProPortal"))
ProPortal$ProPortal_ASV_Ecotype <- as.character(ProPortal$ProPortal_ASV_Ecotype)
Taxonomy <- Taxonomy %>% 
  filter(ProPortal_ASV_Ecotype %in% (NA))

ProPortal <- parse_prefixed_taxonomy(ProPortal)

# This is to create a dataframe for SILVA assigned taxa - this is requried for the "source_database" column
SILVA <- Taxonomy %>% 
  filter(Source_database %in% c("SILVA"))
SILVA <- parse_prefixed_taxonomy(SILVA)

# This is to create a dataframe for PR2 assigned taxa - this is requried for the "source_database" column
PR2 <- Taxonomy %>% 
  filter(Source_database %in% c("PR2"))
PR2 <- PR2 %>% 
  separate(
    Taxonomy,
    c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species"),
    ";",
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(Kingdom = NA_character_, Phylum = NA_character_)
PR2 <- lapply(PR2, gsub, pattern = c(":plas"), replacement = '')
PR2 <- as.data.frame(PR2)

#Join all the Taxonomy together
Taxonomy <- bind_rows(SILVA,PR2,ProPortal)

#Left join taxonomy to asv_sequences
Taxonomy <- Taxonomy %>% 
  left_join(asv_sequences)

#Join asv table to Taxonomy}
counts <- Taxonomy %>% 
  left_join(counts)

# ASV data long. Select sample columns by name rather than by position so the
# new Kingdom annotation (and future taxonomy columns) cannot be mistaken for
# sample counts.
count_annotation_columns <- c(
  "Taxonomy", "Domain", "Kingdom", "Supergroup", "Division", "Subdivision", "Phylum",
  "Class", "Order", "Family", "Genus", "Species", "ProPortal_ASV_Ecotype",
  "plastid_16S_rRNA", "ASV_hash", "ASV", "Source_database"
)
counts_long <- counts %>%
  pivot_longer(
    cols = -all_of(intersect(count_annotation_columns, names(counts))),
    names_to = "SampleID",
    values_to = "Corrected_Sequence_Counts"
  ) %>%
  filter(Corrected_Sequence_Counts != 0)

merged_annotation_columns <- c("ASV_hash", "Taxonomy", "ProPortal_ASV_Ecotype")
counts_dada2 <- counts_dada2_corrected %>%
  pivot_longer(
    cols = -all_of(intersect(merged_annotation_columns, names(counts_dada2_corrected))),
    names_to = "SampleID",
    values_to = "Corrected_dada2_Sequence_Counts"
  ) %>%
  filter(Corrected_dada2_Sequence_Counts != 0)

counts_raw <- counts_uncorrected %>%
  pivot_longer(
    cols = -all_of(intersect(merged_annotation_columns, names(counts_uncorrected))),
    names_to = "SampleID",
    values_to = "Raw_Sequence_Counts"
  ) %>%
  filter(Raw_Sequence_Counts !=0)

asv_long <- counts_long %>% 
  left_join(counts_dada2)

asv_long <- asv_long %>% 
  left_join(counts_raw)

asv_long <- asv_long %>%
  group_by(SampleID, ASV_hash) %>%
  distinct(.keep_all = TRUE)

#Calculate Relative Abundance
asv_long <- asv_long %>% 
  group_by(SampleID) %>% 
  mutate(TC = sum(Corrected_Sequence_Counts)) %>% 
  group_by(SampleID,ASV) %>% 
  mutate(Relative_Abundance = (Corrected_Sequence_Counts/TC))

Check <- asv_long %>% 
  ungroup() %>% 
  group_by(SampleID) %>% 
  summarise(Check = sum(Relative_Abundance))

#Make Sequence Type Column.  Plastid status takes priority over Domain because
# SILVA chloroplast lineages can retain d__Bacteria while PR2 plastid lineages
# normally resolve to Eukaryota.
asv_long <- asv_long %>%
  mutate(
    Sequence_Type = case_when(
      plastid_16S_rRNA == "yes" ~ "Chloroplast_16S",
      Domain %in% c("Bacteria", "Archaea") ~ "Prokaryotic_16S",
      Domain == "Eukaryota" ~ "Eukaryote_18S",
      is.na(Domain) | Domain == "" | Domain == "Unassigned" ~ "Unassigned",
      TRUE ~ "Unassigned"
    )
  )

#Tidy order of columns and what's included in final sheet
asv_long <- asv_long %>% 
  select(SampleID, Domain, Kingdom, Supergroup, Division, Subdivision, Phylum, Class, Order, Family, Genus,
         Species, ProPortal_ASV_Ecotype, Sequence_Type, plastid_16S_rRNA, ASV_hash, ASV, Raw_Sequence_Counts, Corrected_dada2_Sequence_Counts, Corrected_Sequence_Counts, Relative_Abundance, Source_database)

#Last thing is to remove 0's for our use to prevent having such a large file
asv_long <- asv_long %>% filter(!Corrected_Sequence_Counts %in% (0))

#Write asv long
write_tsv(asv_long,snakemake@output[["longdata"]])
