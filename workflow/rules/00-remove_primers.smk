#Removes primers from paired-end raw reads and compresses the trimmed FASTQ files.

#This rule removes the configured primers from paired-end raw reads, discards pairs lacking either expected primer and records trimming statistics.
rule cutadapt: #Defines the Snakemake rule that performs paired-end primer trimming with Cutadapt.
    input:
        r1=lambda wildcards: os.path.join( # Constructs the R1 input path separately for each sample.
            config["rawdatadir"], wildcards.sample + config["R1file_ending"] # Joins the raw-data directory with the sample name and configured R1 filename ending.
        ), # Completes the R1 path and assigns it to the named input.
        r2=lambda wildcards: os.path.join( # Constructs the R2 input path separately for each sample.
            config["rawdatadir"], wildcards.sample + config["R2file_ending"] # Joins the raw-data directory with the sample name and configured R2 filename ending.
        ), # Completes the R2 path and assigns it to the named input.
        #Requires database capability and setup checks to finish before FASTQ processing begins.
        #The shell command explicitly passes only the two read files because the pinned Cutadapt wrapper expects exactly two sequence inputs.
        database_checks=DATABASE_SETUP_CHECKS, # Adds the database readiness files as scheduling dependencies without passing them to Cutadapt.
    output:
        fastq1=temp("results/00-trimmed/{sample}.1.fastq"), #Writes the uncompressed trimmed R1 reads and marks the file for removal after compression.
        fastq2=temp("results/00-trimmed/{sample}.2.fastq"), #Writes the uncompressed trimmed R2 reads and marks the file for removal after compression.
        qc="results/00-trimmed/{sample}.qc.txt", #Saves Cutadapt’s trimming summary as a persistent quality-control report.
    params:
        #The anchored adapter syntax follows Cutadapt’s documented adapter-type conventions.
        adapters="-g ^" + config["fwdPrimer"] + " -G ^" + config["revPrimer"], #Searches for the forward primer at the start of R1 and the reverse primer at the start of R2.
        #These settings control primer matching and removal of untrimmed read pairs.
        extra="--no-indels --pair-filter=any --error-rate=0.2 --discard-untrimmed", #Disallows indels, permits a 20% mismatch rate and discards a pair when either mate lacks its primer.
    log:
        "logs/00-trimming/{sample}.log", #Records Cutadapt’s standard-error messages separately for each sample.
    threads: 4  #Allows Cutadapt to use four processing threads.
    conda:
        "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v7.0.0/bio/cutadapt/pe/environment.yaml" # Uses the pinned paired-end Cutadapt wrapper environment definition.
    shell:
        "cutadapt --cores {threads} {params.adapters} {params.extra} " # Runs Cutadapt with the declared threads, anchored primers and filtering settings.
        "-o {output.fastq1:q} -p {output.fastq2:q} " # Writes the shell-quoted trimmed R1 and R2 output files.
        "{input.r1:q} {input.r2:q} > {output.qc:q} 2> {log:q}" # Supplies the shell-quoted raw reads, sends the summary to the QC file and errors to the log.


#This rule compresses each sample’s primer-trimmed R1 FASTQ with BGZIP for storage and downstream processing.
rule compress_trimmed_output_r1: #Defines the Snakemake rule that compresses the trimmed forward reads.
    input:
        "results/00-trimmed/{sample}.1.fastq", #Uses the temporary trimmed R1 FASTQ produced by the Cutadapt rule.
    output:
        "results/00-trimmed/{sample}.1.fastq.gz", #Writes the BGZIP-compressed R1 FASTQ for the sample.
    threads: 1 #Allocates one processing thread to compression.
    log:
        "logs/bgzip/{sample}_compress-trimmed-r1.log", #Records BGZIP messages and errors separately for each sample’s R1 file.
    wrapper:
        "v7.2.0/bio/bgzip" #Runs the pinned BGZIP wrapper from Snakemake Wrappers version 7.2.0.

#This rule compresses each sample’s primer-trimmed R2 FASTQ with BGZIP for storage and downstream processing.
rule compress_trimmed_output_r2: #Defines the Snakemake rule that compresses the trimmed reverse reads.
    input:
        "results/00-trimmed/{sample}.2.fastq" #Uses the temporary trimmed R2 FASTQ produced by the Cutadapt rule.
    output:
        "results/00-trimmed/{sample}.2.fastq.gz" #Writes the BGZIP-compressed R2 FASTQ for the sample.
    threads: 1 #Allocates one processing thread to compression.
    log:
        "logs/bgzip/{sample}_compress-trimmed-r2.log", #Records BGZIP messages and errors separately for each sample’s R2 file.
    wrapper:
        "v7.2.0/bio/bgzip" #Runs the pinned BGZIP wrapper from Snakemake Wrappers version 7.2.0.
