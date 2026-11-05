#Processes the eukaryotic reads from the BBsplit stage into denoised 18S ASVs, assigns PR2 taxonomy, generates diagnostics, and exports analysis-ready tables.

#This rule creates a QIIME 2 paired-end manifest linking every sample identifier to its split eukaryotic R1 and R2 FASTQ files for initial quality visualization.
rule create_manifest_euk_viz: #Defines the Snakemake rule that builds the paired-end eukaryotic manifest.
    input: 
        r1=expand("results/01-split/{sample}.euk.R1.fastq.gz", sample=samples["sample"]), #Expands the R1 filename pattern across every sample in the sample table.
        r2=expand("results/01-split/{sample}.euk.R2.fastq.gz", sample=samples["sample"]) #Expands the corresponding R2 filename pattern across the same samples.
    output: 
        "results/02-euks/manifest-viz.tsv" #Writes the sample IDs and paired FASTQ paths to a tab-separated QIIME 2 manifest.
    conda: 
        config["qiime2version"] #Uses the configured QIIME 2 environment when running the manifest helper.
    script: 
        "../scripts/E00-create-manifest.sh" #Constructs the paired-end manifest and replaces underscores with hyphens in imported sample IDs.
        
#This rule imports the manifest-listed eukaryotic read pairs into a QIIME 2 artifact for paired-end quality visualization.
rule import_euk_viz: #Defines the Snakemake rule that imports the split eukaryotic FASTQs.
    input:
        "results/02-euks/manifest-viz.tsv" #Uses the paired-end manifest containing sample IDs and R1/R2 file paths.
    output:
        "results/02-euks/18S-viz.qza" #Stores the reads as a `SampleData[PairedEndSequencesWithQuality]` artifact.
    conda: 
        config["qiime2version"] #Uses the configured QIIME 2 environment for importing the reads.
    script: 
        "../scripts/E01-import.sh" #Imports the manifest as Phred33-encoded paired-end FASTQ data.

#This rule generates paired-end quality summaries for the split eukaryotic R1 and R2 reads before truncation and concatenation.
rule visualize_euk_seq_quality:
    input:
        "results/02-euks/18S-viz.qza" #Uses the paired-end eukaryotic reads imported through the visualization manifest.
    output:
        directory("results/02-euks/02-quality-plots-R1-R2/") #Stores QIIME 2’s paired-end demultiplexing summary and quality visualization.
    conda:
        config["qiime2version"]
    script:
        "../scripts/E02-visualize-quality_R1-R2.sh" #Runs `qiime demux summarize` on the paired-end read artifact.


#This rule truncates eukaryotic R1 and R2 reads to their configured lengths, removes reads that are too short and repairs the surviving pairs.
rule bbduk_cut_reads:
    input:
        r1="results/01-split/{sample}.euk.R1.fastq.gz", #Uses each sample’s split eukaryotic forward reads.
        r2="results/01-split/{sample}.euk.R2.fastq.gz" #Uses the corresponding eukaryotic reverse reads.
    output:
        r1=temp("results/02-euks/03-size-selected/{sample}.euk.R1.trimmed.fastq"), #Stores the truncated R1 reads temporarily.
        r2=temp("results/02-euks/03-size-selected/{sample}.euk.R2.trimmed.fastq") #Stores the paired truncated R2 reads temporarily.
    params:
        truncR1=config["trunclens"]["truncR1"], #Passes the configured retained length for R1.
        truncR2=config["trunclens"]["truncR2"] #Passes the configured retained length for R2.
    resources:
        mem_mb=2500, # Requests 2,500 MB of memory.
    conda:
        "../envs/bbmap.yaml"
    log:
        "logs/euk-trimming.{sample}.log"
    script:
        "../scripts/E03-bbduk-cut-reads.sh" #Runs BBDuk separately on both mates and uses BBTools Repair to restore synchronized pairing.


#This rule concatenates each truncated R1 with the reverse-complemented R2, without inserting spacer bases, to create one 18S sequence per pair.
rule fuse_trimmed_euk_seqs:
    input:
        r1=rules.bbduk_cut_reads.output.r1, #Uses the truncated and pair-repaired forward reads.
        r2=rules.bbduk_cut_reads.output.r2 #Uses the corresponding truncated reverse reads.
    output:
        "results/02-euks/04-concatenated/{sample}.euk.concatenated.fastq", #Stores the concatenated sequences as a single-end FASTQ.
    resources:
        mem_mb=2500, # Requests 2,500 MB of memory.
    conda:
        "../envs/bbmap.yaml"
    log:
        "logs/euk-fusing.{sample}.log"
    script:
        "../scripts/E04-fuse-EUKs-withoutNs.sh" #Uses BBTools Fuse with `fusepairs=t` and `pad=0`, so no `N` spacer is added.


#This rule creates the single-end QIIME 2 manifest for all concatenated eukaryotic sequences.
rule create_manifest_euk_concat:
    input:
        files=expand("results/02-euks/04-concatenated/{sample}.euk.concatenated.fastq", sample=samples["sample"]) #Collects the concatenated FASTQ for every sample.
    output:
        "results/02-euks/manifest-concat.tsv" #Writes sample IDs and concatenated FASTQ paths in single-end manifest format.
    conda:
        config["qiime2version"]
    script:
        "../scripts/E05-create-manifest-concat.sh" #Constructs the manifest and replaces underscores with hyphens in imported sample IDs.


#This rule imports the concatenated eukaryotic reads into QIIME 2 as single-end sequence data.
rule import_euk_concat:
    input:
        "results/02-euks/manifest-concat.tsv" #Uses the manifest containing each sample’s concatenated FASTQ path.
    output:
        "results/02-euks/18S-concat.qza" #Stores the reads as a `SampleData[SequencesWithQuality]` artifact.
    conda:
        config["qiime2version"]
    script:
        "../scripts/E06-import-concat.sh" #Imports the data using the Phred33 single-end FASTQ manifest format.


#This rule summarizes the quality distribution of the concatenated sequences before DADA2 denoising.
rule visualize_quality_single_seqs:
    input:
        "results/02-euks/18S-concat.qza" #Uses the imported concatenated 18S reads.
    output:
        directory("results/02-euks/07-quality-plots-concat/") #Stores QIIME 2’s single-end sequence-count and quality visualization.
    conda:
        config["qiime2version"]
    script:
        "../scripts/E07-visualize-quality-single-seqs.sh" #Runs `qiime demux summarize` on the concatenated reads.


#This rule denoises the concatenated 18S reads with single-end DADA2 and produces the ASV table, representative sequences and diagnostics.
rule denoise_euk_dada2:
    input:
        "results/02-euks/18S-concat.qza" #Uses the imported concatenated eukaryotic reads.
    params:
        max_ee=DADA2_EUK["max_ee"], # Sets the maximum expected-error threshold allowed per read.
        trunc_q=DADA2_EUK["trunc_q"], # Sets the quality score at which DADA2 truncates a read.
        pooling_method=DADA2_EUK["pooling_method"], # Selects independent, pseudo-pooled or pooled sample inference.
        chimera_method=DADA2_EUK["chimera_method"], # Selects the configured chimera-detection method.
        min_fold_parent_over_abundance=DADA2_EUK["min_fold_parent_over_abundance"], # Sets the abundance ratio required for a sequence to serve as a chimera parent.
        n_reads_learn=DADA2_EUK["n_reads_learn"] # Sets the number of reads used to learn the error model.
    output:
        directory("results/02-euks/08-DADA2d/"), #Declares the complete DADA2 result directory.
        eukrepseqs="results/02-euks/08-DADA2d/representative_sequences.qza", #Stores one representative sequence for each inferred 18S ASV.
        eukstats="results/02-euks/08-DADA2d/denoising_stats.qza", #Stores per-sample filtering, denoising and chimera statistics.
        euktable="results/02-euks/08-DADA2d/table.qza", #Stores the sample-by-ASV abundance table.
        eukbasetransitions="results/02-euks/08-DADA2d/base_transition_stats.qza" # Stores the learned base-transition statistics.
    conda:
        config["qiime2version"]
    threads: 8 # Allows DADA2 to use eight processing threads.
    resources:
        mem_mb=96000, # Requests 96,000 MB of memory for error learning and ASV inference.
    log:
        "logs/02-denoise-and-export-euk/08-DADA2/DADA2.stderrout"
    script:
        "../scripts/E08-DADA2.sh" #Runs `qiime dada2 denoise-single` without additional left trimming or fixed-length truncation.


#This rule exports the 18S DADA2 artifacts into readable FASTA, TSV and BIOM files while maintaining project-named current copies.
rule export_DADA2_results_euk:
    input:
        "results/02-euks/08-DADA2d/" # Uses the directory containing all eukaryotic DADA2 artifacts.
    params:
        projectName=config["projectName"] # Provides the project name used in exported filenames.
    output:
        directory("results/02-euks/09-DADA2d-plaintext-exports/"), #Stores the complete set of plaintext and BIOM exports.
        lateststats="results/02-euks/09-DADA2d-plaintext-exports/" + config["projectName"] + ".18S.latest_stats.tsv", # Provides the current denoising statistics under a stable project-named filename.
        latestseqs="results/02-euks/09-DADA2d-plaintext-exports/" + config["projectName"] + ".18S.latest_seqs.fasta" # Provides the current representative ASV sequences under a stable project-named filename.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-euk/09_export_DADA2_results/DADA2_export.stderrout"
    script:
        "../scripts/E09-export-DADA2-results.sh" #Exports QIIME 2 artifacts, converts the BIOM table to TSV and archives timestamped copies.


#This rule classifies each 18S ASV with PR2 and applies the established whitespace-normalization workaround.
rule classify_18S_with_PR2:
    input:
        sequences="results/02-euks/08-DADA2d/representative_sequences.qza", # Uses the denoised 18S representative ASV sequences.
        database_ready=rules.ensure_pr2_classifier.output.marker, # Requires the resolved PR2 classifier to be validated and ready.
    params:
        classDB=PR2_CLASSIFIER, # Passes the resolved PR2 classifier path to q2-feature-classifier.
    output:
        classified="results/02-euks/10-classified/" + config["projectName"] + "_PR2.classified.qza", # Stores the original PR2 classification artifact.
        exported_taxonomy=temp(directory("results/02-euks/10-classified/PR2-normalization/exported-taxonomy/")), # Temporarily stores the exported taxonomy TSV.
        taxonomy_metadata=temp("results/02-euks/10-classified/PR2-normalization/taxonomy-as-metadata.qzv"), # Temporarily converts the taxonomy to a metadata visualization.
        exported_metadata=temp(directory("results/02-euks/10-classified/PR2-normalization/exported-metadata/")), # Temporarily stores the re-exported metadata table.
        normalized=temp("results/02-euks/10-classified/PR2-normalization/whitespace-normalized.qza"), # Stores the whitespace-normalized taxonomy artifact.
    conda:
        config["qiime2version"]
    threads: 2 # Uses two classifier workers.
    resources:
        mem_mb=110000, # Requests 110,000 MB of memory for PR2 classification.
    log:
        "logs/02-denoise-and-export-euk/10-classify-18S-with-PR2.log"
    script:
        "../scripts/E10-classify-18S-with-PR2.sh" # Runs `classify-sklearn`, normalizes whitespace through metadata export and reimports the taxonomy.


#This rule exports the whitespace-normalized PR2 classification so its rank prefixes can be removed without losing confidence values.
rule export_18S_PR2_taxonomy_for_normalization:
    input:
        rules.classify_18S_with_PR2.output.normalized # Uses the project-local whitespace-normalized PR2 classification.
    output:
        temp(".snakemake/classification-exports/18S/PR2-raw/taxonomy.tsv") # Stores the exported raw taxonomy table temporarily.
    params:
        export_dir=".snakemake/classification-exports/18S/PR2-raw" # Defines the project-local export directory.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-euk/10-export-PR2-taxonomy.log"
    shell:
        "mkdir -p {params.export_dir:q} && " # Creates the project-local export directory if necessary.
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} " # Exports the PR2 taxonomy artifact as TSV.
        "> {log:q} 2>&1" # Redirects command output and errors to the log.


#This rule removes RESCRIPt-style PR2 rank prefixes while preserving sequence IDs, confidence values and every other exported column.
rule normalize_18S_PR2_rank_prefixes:
    input:
        taxonomy=rules.export_18S_PR2_taxonomy_for_normalization.output[0] # Uses the exported PR2 taxonomy TSV.
    output:
        normalized=temp(".snakemake/classification-exports/18S/PR2/taxonomy.tsv") # Stores the prefix-normalized taxonomy table temporarily.
    log:
        "logs/02-denoise-and-export-euk/10-normalize-PR2-rank-prefixes.log"
    #RESCRIPt’s `edit-taxonomy` action returns only a taxonomy Series and would drop classification confidence.
    #The dedicated TSV transformation therefore modifies only the taxonomy text.
    script:
        "../scripts/normalize_pr2_taxonomy.py" # Removes PR2 rank markers and validates the table structure without parsing numeric-looking fields.


#This rule imports the prefix-normalized PR2 taxonomy back into QIIME 2 for downstream plots and exports.
rule import_normalized_18S_PR2_taxonomy:
    input:
        rules.normalize_18S_PR2_rank_prefixes.output.normalized # Uses the normalized taxonomy TSV.
    output:
        normalized="results/02-euks/10-classified/" + config["projectName"] + "_PR2.normalized.classified.qza" # Stores the final normalized `FeatureData[Taxonomy]` artifact.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-euk/10-import-normalized-PR2-taxonomy.log"
    shell:
        "qiime tools import --type 'FeatureData[Taxonomy]' " # Selects the QIIME 2 taxonomy artifact type.
        "--input-path {input:q} --output-path {output.normalized:q} " # Imports the normalized TSV into the declared artifact.
        "> {log:q} 2>&1" # Redirects command output and errors to the log.


#This rule combines sample metadata with DADA2 retention statistics and upstream prokaryotic/eukaryotic split counts.
rule create_sample_metadata_file_euk:
    input:
        manifest="results/02-euks/manifest-concat.tsv", #Provides the imported sample IDs and concatenated-read paths.
        samplesdottsv="config/samples.tsv", #Provides the user-supplied sample metadata.
        eukstats=rules.export_DADA2_results_euk.output.lateststats, #Provides each sample’s non-chimeric DADA2 read count.
        eukfracpersample="results/" + config["projectName"] + ".eukfrac-per-sample.tsv" # Provides the upstream prokaryotic/eukaryotic read counts and fractions.
    output:
        "results/02-euks/sample-metadata.tsv" #Stores the merged metadata used for the 18S bar plot.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-euk/11-make-sample-metadata-file.log"
    script:
        "../scripts/E11-make-sample-metadata-file.py" #Normalizes underscores to hyphens and merges the tables by sample identifier.


#This rule creates an interactive PR2 taxonomic-composition bar plot for the denoised 18S feature table.
rule make_PR2_euk_barplots:
    input:
        euktable="results/02-euks/08-DADA2d/table.qza", #Uses the complete sample-by-18S-ASV abundance table.
        euktax=rules.import_normalized_18S_PR2_taxonomy.output.normalized, # Uses the normalized PR2 assignments.
        eukmetadata="results/02-euks/sample-metadata.tsv" #Uses the combined sample metadata and processing statistics.
    params:
        projectName=config["projectName"] # Provides the project name for the archived visualization filename.
    output:
        directory("results/02-euks/12-PR2-barplots/") # Stores the QIIME 2 bar-plot visualization and a timestamped copy.
    conda:
        config["qiime2version"]
    script:
        "../scripts/E12-make-barplot.sh" #Runs `qiime taxa barplot` using the table, PR2 taxonomy and sample metadata.


#This rule exports the normalized PR2 taxonomy and complete 18S feature table in preparation for taxonomy-enriched BIOM and TSV outputs.
rule export_tax_convert_biom_euk:
    input:
        PR2classified=rules.import_normalized_18S_PR2_taxonomy.output.normalized, # Uses the final normalized PR2 classification.
        all18Stable=rules.denoise_euk_dada2.output.euktable # Uses the complete denoised 18S feature table.
    output:
        PR2taxdir=temp(directory("results/02-euks/15-exports/" + config["projectName"] + ".taxonomy-PR2/")), # Temporarily stores the raw taxonomy export.
        PR2taxfile="results/02-euks/15-exports/" + config["projectName"] + ".taxonomy-PR2.tsv", # Stores the persistent PR2 taxonomy table.
        all18Stablebiom=temp("results/02-euks/15-exports/" + config["projectName"] + ".all-18S-seqs.biom") # Stores the complete BIOM feature table temporarily.
    conda:
        config["qiime2version"]
    script:
        "../scripts/E15a-generate-biom-tables.sh" #Exports both artifacts and formats the taxonomy header for BIOM metadata import.


#This rule attaches the normalized PR2 taxonomy to each observation in the complete 18S BIOM table.
rule add_tax_to_biom_euk:
    input:
        PR2taxfile="results/02-euks/15-exports/" + config["projectName"] + ".taxonomy-PR2.tsv", # Uses the exported ASV-to-PR2 taxonomy table.
        all18Stablebiom="results/02-euks/15-exports/" + config["projectName"] + ".all-18S-seqs.biom" # Uses the complete 18S abundance table.
    output:
        all18StablebiomPR2tax=temp("results/02-euks/15-exports/" + config["projectName"] + ".all-18S-seqs.with-PR2-tax.biom") # Stores the taxonomy-enriched BIOM table temporarily.
    conda:
        config["qiime2version"]
    script:
        "../scripts/E15b-add-tax-to-biom.sh" #Uses `biom add-metadata` to attach semicolon-separated taxonomy to each ASV.


#This rule converts the taxonomy-enriched 18S BIOM table into a readable tab-separated abundance table.
rule export_biom_tsv_euk:
    input:
        all18StablebiomPR2tax="results/02-euks/15-exports/" + config["projectName"] + ".all-18S-seqs.with-PR2-tax.biom" # Uses the complete taxonomy-enriched BIOM table.
    output:
        all18StablebiomPR2taxtsv="results/02-euks/15-exports/" + config["projectName"] + ".all-18S-seqs.with-PR2-tax.tsv" # Writes all ASV abundances with PR2 taxonomy in the final column.
    conda:
        config["qiime2version"]
    script:
        "../scripts/E15c-export-tax-tsvs.sh" #Runs `biom convert --to-tsv` with taxonomy as the observation header metadata.


#This rule divides the complete taxonomy-enriched 18S table into Metazoa-only and Metazoa-excluded tables.
rule subset_euk_taxonomy_tables:
    input:
        all18SPR2="results/02-euks/15-exports/" + config["projectName"] + ".all-18S-seqs.with-PR2-tax.tsv" # Uses the complete PR2-annotated 18S abundance table.
    output:
        excludemetazoaPR2tablebiomtaxtsv="results/02-euks/15-exports/" + config["projectName"] + ".exclude_Metazoa_PR2_filtered_table.with-tax.tsv", # Stores all ASVs whose taxonomy does not contain `Metazoa`.
        includemetazoaPR2tablebiomtaxtsv="results/02-euks/15-exports/" + config["projectName"] + ".include_Metazoa_PR2_filtered_table.with-tax.tsv" # Stores only ASVs whose taxonomy contains `Metazoa`.
    script:
        "../scripts/E15d-subset-tax-tsvs.py" # Preserves valid BIOM TSV headers even when either taxonomic subset is empty.
