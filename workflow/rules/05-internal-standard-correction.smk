#Activates the optional internal-standard workflow. Prepares reference FASTAs, identifies matching 16S/18S ASVs with BLAST, applies the configured correction methods, removes standard-derived ASVs from biological data, and generates diagnostic plots.

#This block prepares, detects and applies the configured internal standards only when internal-standard correction is enabled.
if USE_INTERNAL_STANDARDS: #Includes the following definitions and rules only when internal standards are enabled.

    INTERNAL_STANDARD_FASTAS = expand( #Builds the expected reference-FASTA path for every configured standard.
        "config/intstd_fastas/{standard}.fasta", #Defines the filename pattern for each standard FASTA.
        standard=INTERNAL_STANDARD_IDS, #Expands the pattern using the configured standard identifiers.
    )
    INTERNAL_STANDARD_ASV_LISTS = expand( #Builds the expected matched-ASV list path for every standard.
        "config/intstd_fastas/{standard}.asvs.txt", #Defines the filename pattern for each standard’s ASV list.
        standard=INTERNAL_STANDARD_IDS, #Creates one ASV-list path per configured standard.
    )
    INTERNAL_STANDARD_METHOD_TABLES = expand( #Builds the output path for every configured correction method.
        "results/05-internal-std-corrected/" + config["projectName"] + ".asv_table_{method}.tsv", #Includes the project and correction-method names in each table filename.
        method=INTERNAL_STANDARD_METHOD_STEMS, #Creates one table path per correction method.
    )


    #This rule converts the internal-standard definitions into marker-specific reference sequences suitable for matching against pipeline ASVs.
    rule prepare_internal_standard_fastas: #Defines the reference-FASTA preparation rule.
        input:
            table="config/internal_stds.tsv" #Uses the table containing the internal-standard identifiers and reference sequences.
        output:
            fastas=INTERNAL_STANDARD_FASTAS #Creates one prepared FASTA for every configured standard.
        params:
            standard_ids=INTERNAL_STANDARD_IDS, #Passes the ordered standard identifiers to the preparation script.
            forward_primer=config["fwdPrimer"], #Passes the forward primer used by the workflow.
            reverse_primer=config["revPrimer"], #Passes the reverse primer used by the workflow.
            trunc_r1=config["trunclens"]["truncR1"], #Passes the configured forward-read truncation length.
            trunc_r2=config["trunclens"]["truncR2"] #Passes the configured reverse-read truncation length.
        script:
            "../scripts/prepare_internal_standard_fastas.py" #Extracts the expected 16S or 18S amplicon representation for each standard.


    #This rule creates a nucleotide BLAST database from each prepared internal-standard reference FASTA.
    rule prepare_16S_BLASTdb: #Defines the BLAST-database preparation rule.
        input:
            fasta="config/intstd_fastas/{standard}.fasta" #Uses the prepared reference sequence for the selected standard.
        output:
            nhr="config/intstd_fastas/{standard}.fasta.nhr" #Tracks the BLAST nucleotide-header file as evidence that the database was created.
        wildcard_constraints:
            standard="|".join(re.escape(value) for value in INTERNAL_STANDARD_IDS) #Restricts the wildcard to configured standard identifiers.
        conda:
            "../envs/blast-env.yaml" #Provides BLAST+ and its required dependencies.
        shell:
            "makeblastdb -dbtype nucl -in {input.fasta:q} && test -s {output.nhr:q}" #Builds the nucleotide database and verifies that its tracked index file is nonempty.


    #This rule searches both 16S and 18S representative sequences for near-exact, full-coverage matches to each internal standard.
    rule identify_intsd_ASVS: #Defines the internal-standard ASV identification rule.
        input:
            latestseqs_16s="results/02-proks/04-DADA2d-plaintext-exports/" + config["projectName"] + ".16S.latest_seqs.fasta", #Provides the representative 16S ASV sequences.
            latestseqs_18s="results/02-euks/09-DADA2d-plaintext-exports/" + config["projectName"] + ".18S.latest_seqs.fasta", #Provides the representative 18S ASV sequences.
            database="config/intstd_fastas/{standard}.fasta.nhr" #Requires the selected standard’s BLAST database to be ready.
        output:
            matches_16s="config/intstd_fastas/{standard}.16S.asvs.outfmt6.tsv", #Stores detailed 16S BLAST matches in tabular format.
            matches_18s="config/intstd_fastas/{standard}.18S.asvs.outfmt6.tsv", #Stores detailed 18S BLAST matches in tabular format.
            matches="config/intstd_fastas/{standard}.asvs.outfmt6.tsv", #Stores the combined 16S and 18S BLAST results.
            asvs="config/intstd_fastas/{standard}.asvs.txt" #Stores the unique ASV identifiers assigned to the standard.
        wildcard_constraints:
            standard="|".join(re.escape(value) for value in INTERNAL_STANDARD_IDS) #Limits processing to configured internal standards.
        params:
            database=lambda wildcards: f"config/intstd_fastas/{wildcards.standard}.fasta" #Builds the BLAST database prefix for the selected standard.
        conda:
            "../envs/blast-env.yaml" #Provides the BLAST and standard Unix command-line tools.
        shell:
            "blastn -query {input.latestseqs_16s:q} -db {params.database:q} " #Searches the 16S ASVs against the selected standard.
            "-outfmt 6 -perc_identity 99 -qcov_hsp_perc 100 > {output.matches_16s:q} " #Retains matches with at least 99% identity and complete query coverage.
            "&& blastn -query {input.latestseqs_18s:q} -db {params.database:q} " #Searches the 18S ASVs against the same standard.
            "-outfmt 6 -perc_identity 99 -qcov_hsp_perc 100 > {output.matches_18s:q} " #Applies the same identity and query-coverage thresholds.
            "&& cat {output.matches_16s:q} {output.matches_18s:q} > {output.matches:q} " #Combines the marker-specific BLAST results.
            "&& cut -f1 {output.matches:q} | sort -u > {output.asvs:q}" #Extracts a unique, sorted list of matching ASV identifiers.


    #This rule applies the configured internal-standard corrections, removes standard-derived ASVs from biological data and creates recovery diagnostics.
    rule intstd_correct_data: #Defines the internal-standard correction and reporting rule.
        input:
            asv_table="results/04-formatted/" + config["projectName"] + ".long_data.tsv", #Uses the combined long-format 16S and 18S abundance data.
            standard_asvs=INTERNAL_STANDARD_ASV_LISTS, #Uses the matched-ASV list for every configured standard.
            isd="config/internal_stds.tsv", #Provides the standard definitions and reference information.
            isd_added="config/samples.tsv" #Provides the amount of each standard added to each sample.
        params:
            standard_ids=INTERNAL_STANDARD_IDS, #Passes the ordered internal-standard identifiers.
            standard_slots=INTERNAL_STANDARD_SLOTS, #Passes the internal column names used for the configured standards.
            method_stems=INTERNAL_STANDARD_METHOD_STEMS #Passes the names of the correction methods to generate.
        output:
            corrected=ISD_CORRECTED_LONG_TABLE, #Stores the corrected long-format data including internal-standard records.
            filtered=ISD_FILTERED_LONG_TABLE, #Stores corrected data with internal-standard ASVs removed.
            method_tables=INTERNAL_STANDARD_METHOD_TABLES, #Stores a separate wide ASV table for each correction method.
            recovery_plot="results/06-figures/" + config["projectName"] + ".recovery_ratios.pdf", #Stores the internal-standard recovery diagnostic as a PDF.
            domain_plot="results/06-figures/" + config["projectName"] + ".Domain_by_sampleID.pdf", #Stores domain-level corrected composition by sample as a PDF.
            recovery_plot_png="results/06-figures/" + config["projectName"] + ".recovery_ratios.png", #Stores a PNG copy of the recovery diagnostic.
            domain_plot_png="results/06-figures/" + config["projectName"] + ".Domain_by_sampleID.png" #Stores a PNG copy of the domain-composition plot.
        conda:
            "../envs/r-tidyverse-2.0.0.yml" #Provides R, tidyverse and the plotting dependencies.
        log:
            "logs/05-internal-std-correction/prepare_long_data_corrected.log" #Specifies the log file for correction and plotting messages.
        script:
            "../scripts/correct-w-isds-snakemake-v5.R" #Calculates recovery-based corrections, filters standard ASVs and creates the tables and figures.
