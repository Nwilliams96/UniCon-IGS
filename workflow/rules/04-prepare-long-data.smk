#This rule reshapes the uncorrected and corrected 16S/18S abundance tables into analysis-ready long format and creates a lookup table linking each ASV identifier to its sequence.
rule prepare_long_data: #Defines the rule that prepares merged ASV data for downstream analysis.
    input:
        mergedtableuncorrected="results/03-merged/" + config["projectName"] + ".merged_uncorrected.tsv", #Uses the merged 16S/18S table containing the original read counts.
        mergedtabledada2="results/03-merged/" + config["projectName"] + ".merged_dada2_corrected.tsv", #Uses the merged table corrected for DADA2 read retention.
        mergedtabledada218Scorrected="results/03-merged/" + config["projectName"] + ".merged_dada2_18S_corrected.tsv", #Uses the table additionally corrected for the observed 16S/18S amplification balance.
        fasta16S="results/02-proks/04-DADA2d-plaintext-exports/" + config["projectName"] + ".16S.latest_seqs.fasta", #Provides the representative 16S ASV sequences.
        fasta18S="results/02-euks/09-DADA2d-plaintext-exports/" + config["projectName"] + ".18S.latest_seqs.fasta" #Provides the representative 18S ASV sequences.
    output:
        longdata="results/04-formatted/" + config["projectName"] + ".long_data.tsv", #Stores the abundance tables in a combined long format suitable for analysis and plotting.
        asvsequences="results/04-formatted/" + config["projectName"] + ".asv_sequences.tsv" #Stores the sequence associated with each 16S and 18S ASV identifier.
    conda:
        "../envs/r-tidyverse-2.0.0.yml" #Provides R and tidyverse packages for reshaping and combining the data.
    log:
        "logs/04-formatting/prepare_long_data.log" #Specifies the log file for this formatting step.
    script:
        "../scripts/long-data-preparation.R" #Reshapes the abundance tables and extracts the ASV sequences from the FASTA files.
