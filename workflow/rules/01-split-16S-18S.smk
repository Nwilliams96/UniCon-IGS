#Separates primer-trimmed reads into prokaryotic and eukaryotic groups, restores separate R1/R2 files, and summarizes the split.

#This rule maps each sample’s primer-trimmed read pairs against the shared BBsplit database and separates them into prokaryotic and eukaryotic FASTQ files.
rule bbsplit_prok_euk: #Defines the Snakemake rule that performs reference-based prokaryotic/eukaryotic read splitting.
    input:
        database=BBSPLIT_DB_DIR, # Uses the prepared BBsplit index containing the prokaryotic and eukaryotic references.
        r1="results/00-trimmed/{sample}.1.fastq", #Uses the sample’s primer-trimmed forward reads.
        r2="results/00-trimmed/{sample}.2.fastq" #Uses the corresponding primer-trimmed reverse reads.
    output:
        prok="results/01-split/{sample}.prok.fastq", #Stores read pairs assigned to the SILVA 132 prokaryotic reference.
        euk="results/01-split/{sample}.euk.fastq", #Stores read pairs assigned to the combined SILVA 132 and PR2 eukaryotic reference.
    conda:
        "../envs/bbmap.yaml" #Uses the project’s Conda environment containing BBTools and BBsplit.
    threads: 8 # Allows BBsplit to use eight processing threads.
    resources:
        #The shared BBsplit reference used approximately 55 GB of Java heap during the SR03 run.
        #Reserving 64 GB prevents multiple memory-heavy split jobs from running simultaneously within the 128 GB CARC allocation.
        mem_mb=64000, # Requests 64,000 MB of memory for this rule.
    log:
        "logs/01-splitting/{sample}_bbsplit.log" #Records mapping and splitting messages separately for each sample.
    shell:
        "bbsplit.sh threads={threads} usequality=f qtrim=f minratio=0.30 minid=0.30 pairedonly=f path={input.database} in={input.r1} in2={input.r2} out_SILVA_132_PROK.cdhit95pc={output.prok} out_SILVA_132_and_PR2_EUK.cdhit95pc={output.euk} 2>&1 | tee -a {log}" # Maps without quality-based seed selection or trimming, uses permissive 0.30 mapping thresholds, allows mapping without requiring both mates to map and appends all messages to the log.
        
#This rule separates each sample’s interleaved eukaryotic BBsplit output into compressed forward and reverse FASTQ files.
rule deinterleave_split_reads_euk: #Defines the Snakemake rule that deinterleaves the eukaryotic read pairs.
    input:
        "results/01-split/{sample}.euk.fastq" #Uses the eukaryotic FASTQ produced by the BBsplit rule.
    output:
        out="results/01-split/{sample}.euk.R1.fastq.gz", #Writes the eukaryotic forward reads as a gzip-compressed FASTQ.
        out2="results/01-split/{sample}.euk.R2.fastq.gz" #Writes the corresponding eukaryotic reverse reads as a gzip-compressed FASTQ.
    conda:
        "../envs/bbmap.yaml" #Uses the project’s Conda environment containing BBTools Reformat.
    threads: 4 #Allows Reformat to use four processing threads.
    resources:
        mem_mb=4000, #Requests 4,000 MB of memory for deinterleaving and compression.
    shell:
        "reformat.sh threads={threads} in={input:q} out1={output.out:q} out2={output.out2:q} overwrite=t" # Splits the shell-quoted interleaved input into R1 and R2 outputs and permits replacement of existing output files.
        
#This rule separates each sample’s interleaved prokaryotic BBsplit output into compressed forward and reverse FASTQ files.
rule deinterleave_split_reads_prok: #Defines the Snakemake rule that deinterleaves the prokaryotic read pairs.
    input:
        "results/01-split/{sample}.prok.fastq" #Uses the prokaryotic FASTQ produced by the BBsplit rule.
    output:
        out="results/01-split/{sample}.prok.R1.fastq.gz", #Writes the prokaryotic forward reads as a gzip-compressed FASTQ.
        out2="results/01-split/{sample}.prok.R2.fastq.gz" #Writes the corresponding prokaryotic reverse reads as a gzip-compressed FASTQ.
    conda:
        "../envs/bbmap.yaml" #Uses the project’s Conda environment containing BBTools Reformat.
    threads: 4 #Allows Reformat to use four processing threads.
    resources:
        mem_mb=4000, #Requests 4,000 MB of memory for deinterleaving and compression.
    shell:
        "reformat.sh threads={threads} in={input:q} out1={output.out:q} out2={output.out2:q} overwrite=t" # Splits the shell-quoted interleaved input into R1 and R2 outputs and permits replacement of existing output files.

#This rule counts the prokaryotic and eukaryotic read pairs assigned to each sample and calculates the sample’s eukaryotic fraction.
rule count_seqs: #Defines the Snakemake rule that creates a per-sample split-read count record.
    input:
        prok="results/01-split/{sample}.prok.R1.fastq.gz", #Uses the prokaryotic R1 FASTQ, where each R1 record represents one retained read pair.
        euk="results/01-split/{sample}.euk.R1.fastq.gz" #Uses the eukaryotic R1 FASTQ, where each R1 record represents one retained read pair.
    output:
        eukfrac="results/01-split/counts/{sample}.eukfrac", #Stores the prokaryotic count, eukaryotic count, combined count and eukaryotic fraction.
    params:
        "{sample}" #Passes the sample identifier as the script’s first unnamed parameter.
    conda:
        config["qiime2version"] #Uses the configured environment containing the command-line utilities required by the counting script.
    script:
        "../scripts/count-seqs.sh" #Counts FASTQ records and calculates eukaryotic reads divided by total assigned reads.

#This rule combines every sample’s split-read counts and eukaryotic fraction into one project-level per-sample table.
rule concatenate_eukfrac_data: #Defines the Snakemake rule that assembles the per-sample eukaryotic-fraction summary.
    input:
        eukfrac=expand("results/01-split/counts/{sample}.eukfrac", sample=samples["sample"]) #Expands the filename pattern across the sample table’s complete ordered list of sample identifiers.
    output:
        eukfracpersample="results/" + config["projectName"] + ".eukfrac-per-sample.tsv" # Names the tab-separated output using the configured project name.
    conda:
        config["qiime2version"] #Uses the configured environment containing the command-line utilities required by the concatenation script.
    script:
        "../scripts/concat-eukfrac.sh" #Adds the column header and combines the individual sample records into one table.

#This rule aggregates the per-sample split-read counts into a project-level table containing total prokaryotic and eukaryotic reads.
rule calc_eukfrac_overall: #Defines the Snakemake rule that calculates project-wide split-read totals.
    input:
        eukfracpersample="results/" + config["projectName"] + ".eukfrac-per-sample.tsv" # Uses the project-named table containing counts and eukaryotic fractions for every sample.
    output:
        eukfracall="results/" + config["projectName"] + ".eukfrac-all.tsv" # Names the tab-separated summary using the configured project name.
    script:
        "../scripts/calc-EUK-fraction.sh" #Sums the prokaryotic and eukaryotic counts across samples and writes the two totals.
