#This rule combines taxonomy-annotated 16S and 18S ASV tables, adds ProPortal annotations, and exports raw, DADA2-adjusted and cross-marker-corrected abundance tables with the applied correction factors.
rule merge_prok_euk: #Defines the rule that merges and corrects the prokaryotic and eukaryotic ASV tables.
    input:
        raw16S="results/02-proks/10-exports/" + config["projectName"] + ".all-16S-seqs.with-tax.tsv", #Uses the complete taxonomy-annotated 16S abundance table.
        raw18S="results/02-euks/15-exports/" + config["projectName"] + ".all-18S-seqs.with-PR2-tax.tsv", #Uses the complete PR2-annotated 18S abundance table.
        stats16S="results/02-proks/04-DADA2d-plaintext-exports/" + config["projectName"] + ".16S.latest_stats.tsv", #Provides per-sample 16S input and non-chimeric read counts for DADA2-retention correction.
        stats18S="results/02-euks/09-DADA2d-plaintext-exports/" + config["projectName"] + ".18S.latest_stats.tsv", #Provides per-sample 18S input and non-chimeric read counts for DADA2-retention correction.
        read_summary="results/" + config["projectName"] + ".eukfrac-all.tsv", #Provides the overall read fractions assigned to the prokaryotic and eukaryotic branches.
        amplicon_concentrations=AMPLICON_CONCENTRATIONS_FILE, #Provides the starting 16S and 18S amplicon concentrations used to calculate cross-marker correction factors.
        proportalclassification="results/02-proks/10-exports/" + config["projectName"] + ".Synechococcales.proportal-classified.tsv" #Provides optional ProPortal ecotype annotations for matching 16S ASVs.
    output:
        mergedtableuncorrected="results/03-merged/" + config["projectName"] + ".merged_uncorrected.tsv", #Stores the combined raw 16S and 18S ASV counts without abundance correction.
        mergedtabledada2="results/03-merged/" + config["projectName"] + ".merged_dada2_corrected.tsv", #Stores counts adjusted for the fraction of reads retained by DADA2 in each branch.
        mergedtabledada218Scorrected="results/03-merged/" + config["projectName"] + ".merged_dada2_18S_corrected.tsv", #Stores DADA2-adjusted counts after also correcting the observed 16S/18S balance to the starting molar balance.
        correction_factors="results/03-merged/" + config["projectName"] + ".16S_18S_correction_factors.tsv" #Records the concentrations, read fractions and correction factors applied to both marker branches.
    conda:
        "../envs/r-tidyverse-2.0.0.yml" #Provides R, tidyverse and data.table for the merge and correction calculations.
    log:
        "logs/03-merging/merge_prok_euk.log" #Declares the log path available to the R helper.
    script:
        "../scripts/correct_16S_18S_ASV-snakemake-v1.R" #Calculates the correction factors, merges the tables, adds ProPortal annotations and writes all four outputs.
