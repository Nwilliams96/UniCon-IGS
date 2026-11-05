#!/usr/bin/env bash

#This script searches the representative ASVs against three internal-standard BLAST databases and extracts the identifiers of matching ASVs.
#Requires complete query coverage and at least 99% identity; for an approximately 400-base amplicon, this permits roughly four mismatches.

db=config/intstd_fastas/`basename ${snakemake_input[intstd1]} .nhr` #Removes the `.nhr` extension from the first database index and reconstructs its BLAST database prefix.
echo blastn -query ${snakemake_input[latestseqs]} -db $db -outfmt 6 -perc_identity 99 -qcov_hsp_perc 100 > ${snakemake_output[intstd1tsv]} #Writes the first BLAST command to the result path before the next command replaces it.
blastn -query ${snakemake_input[latestseqs]} -db $db -outfmt 6 -perc_identity 99 -qcov_hsp_perc 100 > ${snakemake_output[intstd1tsv]} #Searches the ASVs against the first standard and writes full-coverage matches with at least 99% identity in tabular format.
cut -f1 ${snakemake_output[intstd1tsv]} > ${snakemake_output[intstd1asvs]} #Extracts the query ASV identifier from each first-standard match.

db=config/intstd_fastas/`basename ${snakemake_input[intstd2]} .nhr` #Constructs the BLAST database prefix for the second internal standard.
blastn -query ${snakemake_input[latestseqs]} -db $db -outfmt 6 -perc_identity 99 -qcov_hsp_perc 100 > ${snakemake_output[intstd2tsv]} #Searches the ASVs against the second standard using the same thresholds.
cut -f1 ${snakemake_output[intstd2tsv]} > ${snakemake_output[intstd2asvs]} #Extracts the ASV identifiers matching the second standard.

db=config/intstd_fastas/`basename ${snakemake_input[intstd3]} .nhr` #Constructs the BLAST database prefix for the third internal standard.
blastn -query ${snakemake_input[latestseqs]} -db $db -outfmt 6 -perc_identity 99 -qcov_hsp_perc 100 > ${snakemake_output[intstd3tsv]} #Searches the ASVs against the third standard using the same thresholds.
cut -f1 ${snakemake_output[intstd3tsv]} > ${snakemake_output[intstd3asvs]} #Extracts the ASV identifiers matching the third standard.