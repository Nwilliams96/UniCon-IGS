#processes prokaryotic 16S reads through import, quality assessment, DADA2 denoising, SILVA/PR2 taxonomy classification, and plastid reconciliation. It then generates taxonomically filtered feature tables, BIOM/TSV exports, species assignments, and supplemental ProPortal classifications.

#This rule creates a paired-end QIIME 2 manifest linking every sample to its split prokaryotic R1 and R2 FASTQ files.
rule create_manifest_prok: #Defines the rule that builds the prokaryotic read manifest.
    input:
        r1=expand("results/01-split/{sample}.prok.R1.fastq.gz", sample=samples["sample"]), #Collects the prokaryotic forward-read FASTQs for every sample.
        r2=expand("results/01-split/{sample}.prok.R2.fastq.gz", sample=samples["sample"]) #Collects the corresponding prokaryotic reverse-read FASTQs.
    output:
        "results/02-proks/manifest.tsv" #Writes sample IDs and absolute paired-read paths in QIIME 2 manifest format.
    conda:
        config["qiime2version"]
    script:
        "../scripts/P00-create-manifest.sh" #Builds the paired-end manifest and converts underscores to hyphens in imported sample IDs.

#This rule imports the manifest-listed prokaryotic read pairs into a QIIME 2 artifact for denoising and quality inspection.
rule import_prok: #Defines the rule that imports the split prokaryotic FASTQs.
    input:
        "results/02-proks/manifest.tsv" #Uses the paired-end manifest created for the prokaryotic reads.
    output:
        "results/02-proks/16S.qza" #Stores the reads as Phred33-encoded paired-end QIIME 2 data.
    conda:
        config["qiime2version"]
    script:
        "../scripts/P01-import.sh" #Imports the manifest as `SampleData[PairedEndSequencesWithQuality]`.

#This rule summarizes the sequence counts and positional quality of the imported prokaryotic R1 and R2 reads.
rule visualize_prok_seq_quality: #Defines the rule that generates the paired-end 16S quality visualization.
    input:
        "results/02-proks/16S.qza" #Uses the imported paired-end prokaryotic reads.
    output:
        directory("results/02-proks/02-quality-plots-R1-R2/") #Stores the QIIME 2 demultiplexed-sequence summary and quality plot.
    conda:
        config["qiime2version"]
    script:
        "../scripts/P02-visualize-quality_R1-R2.sh" #Runs `qiime demux summarize` on the imported paired reads.

#This rule performs the canonical paired-end DADA2 workflow, retaining only successfully merged 16S read pairs as ASVs.
rule denoise_prok_dada2: #Defines the standard prokaryotic denoising rule.
    input:
        "results/02-proks/16S.qza" #Uses the imported paired-end prokaryotic reads.
    params:
        truncR1=DADA2_PROK["trunc_len_f"], #Sets the retained forward-read length from the effective 16S DADA2 configuration.
        truncR2=DADA2_PROK["trunc_len_r"], #Sets the retained reverse-read length.
        max_ee_f=DADA2_PROK["max_ee_f"], #Sets the maximum expected errors allowed in a forward read.
        max_ee_r=DADA2_PROK["max_ee_r"], #Sets the maximum expected errors allowed in a reverse read.
        trunc_q=DADA2_PROK["trunc_q"], #Sets the quality threshold at which a read is truncated.
        min_overlap=DADA2_PROK["min_overlap"], #Sets the minimum overlap required to merge a read pair.
        max_merge_mismatch=DADA2_PROK["max_merge_mismatch"], #Sets the maximum mismatches permitted within the overlap.
        trim_overhang=DADA2_PROK["trim_overhang"], #Controls whether overhanging sequence is trimmed during paired-read merging.
        retain_unmerged=False, #Excludes pairs that cannot be merged from the canonical ASV outputs.
        pooling_method=DADA2_PROK["pooling_method"], #Selects how samples are pooled during sequence-variant inference.
        chimera_method=DADA2_PROK["chimera_method"], #Selects the configured chimera-detection method.
        min_fold_parent_over_abundance=DADA2_PROK["min_fold_parent_over_abundance"], #Sets the minimum abundance ratio required of a potential chimera parent.
        n_reads_learn=DADA2_PROK["n_reads_learn"] #Sets the number of reads used to learn the error model.
    output:
        directory("results/02-proks/03-DADA2d/"), #Declares the complete canonical DADA2 result directory.
        prokrepseqs="results/02-proks/03-DADA2d/representative_sequences.qza", #Stores one representative sequence for each merged 16S ASV.
        prokstats="results/02-proks/03-DADA2d/denoising_stats.qza", #Stores per-sample filtering, merging and chimera-removal statistics.
        proktable="results/02-proks/03-DADA2d/table.qza", #Stores the canonical sample-by-ASV abundance table.
        prokbasetransitions="results/02-proks/03-DADA2d/base_transition_stats.qza" #Stores the learned base-transition statistics.
    conda:
        config["qiime2version"]
    threads: 8 #Allows DADA2 to use eight processing threads.
    resources:
        mem_mb=96000, #Requests 96,000 MB of memory for error learning and ASV inference.
    log:
        "logs/02-denoise-and-export-prok/03-DADA2/DADA2.stderrout" #Records the canonical DADA2 command output and errors.
    script:
        "../scripts/P03-DADA2.sh" #Runs paired-end DADA2 with unmerged-pair retention disabled.

#This rule repeats paired-end DADA2 while retaining unmerged linked read pairs for a separate supplemental 16S rescue analysis.
rule denoise_prok_dada2_with_unmerged: #Defines the experimental retained-pair denoising rule.
    input:
        "results/02-proks/16S.qza" #Uses the same imported prokaryotic read pairs as the canonical run.
    params:
        truncR1=DADA2_PROK["trunc_len_f"], #Uses the canonical forward-read truncation length.
        truncR2=DADA2_PROK["trunc_len_r"], #Uses the canonical reverse-read truncation length.
        max_ee_f=DADA2_PROK["max_ee_f"], #Uses the canonical forward-read expected-error limit.
        max_ee_r=DADA2_PROK["max_ee_r"], #Uses the canonical reverse-read expected-error limit.
        trunc_q=DADA2_PROK["trunc_q"], #Uses the canonical quality-truncation threshold.
        min_overlap=DADA2_PROK["min_overlap"], #Uses the canonical minimum merge overlap.
        max_merge_mismatch=DADA2_PROK["max_merge_mismatch"], #Uses the canonical overlap mismatch limit.
        trim_overhang=DADA2_PROK["trim_overhang"], #Uses the canonical overhang-trimming behavior.
        retain_unmerged=True, #Retains pairs that cannot merge as linked sequence features.
        pooling_method=DADA2_PROK["pooling_method"], #Uses the canonical sample-pooling method.
        chimera_method=DADA2_PROK["chimera_method"], #Uses the canonical chimera-detection method.
        min_fold_parent_over_abundance=DADA2_PROK["min_fold_parent_over_abundance"], #Uses the canonical parent-abundance threshold.
        n_reads_learn=DADA2_PROK["n_reads_learn"] #Uses the canonical number of reads for error-model learning.
    output:
        directory("results/02-proks/03-DADA2d-with-unmerged/"), #Declares the retained-pair DADA2 result directory.
        prokrepseqs="results/02-proks/03-DADA2d-with-unmerged/representative_sequences.qza", #Stores merged and retained linked representative features.
        prokstats="results/02-proks/03-DADA2d-with-unmerged/denoising_stats.qza", #Stores statistics for the retained-pair DADA2 run.
        proktable="results/02-proks/03-DADA2d-with-unmerged/table.qza", #Stores the abundance table that includes retained linked features.
        prokbasetransitions="results/02-proks/03-DADA2d-with-unmerged/base_transition_stats.qza" #Stores the learned base-transition statistics for this run.
    conda:
        config["qiime2version"]
    threads: 8 #Allows DADA2 to use eight processing threads.
    resources:
        mem_mb=96000, #Requests 96,000 MB of memory.
    log:
        "logs/02-denoise-and-export-prok/03_DADA2_with_unmerged/DADA2.stderrout" #Records output and errors from the retained-pair run.
    script:
        "../scripts/P03-DADA2.sh" #Runs paired-end DADA2 with unmerged-pair retention enabled.

#This rule classifies all representative features from the retained-pair DADA2 run with the validated SILVA classifier.
rule classify_unmerged_16S_with_SILVA: #Defines SILVA classification for the supplemental retained-pair analysis.
    input:
        sequences=rules.denoise_prok_dada2_with_unmerged.output.prokrepseqs, #Uses the representative features from the DADA2 run that retained unmerged pairs.
        # A local readiness marker schedules download/training/validation while
        # keeping the shared classifier outside Snakemake's output cleanup.
        database_ready=rules.ensure_silva_classifier.output.marker, #Requires the SILVA classifier to be validated before classification begins.
    params:
        classDB=SILVA_CLASSIFIER, #Passes the resolved SILVA classifier path to q2-feature-classifier.
    output:
        classified="results/02-proks/11-unmerged-16S-rescue/SILVA.classified.qza" #Stores SILVA taxonomy and confidence for the retained-run features.
    conda:
        config["qiime2version"]
    threads: 8 #Uses eight sklearn classification workers.
    resources:
        mem_mb=64000, #Requests 64,000 MB of memory for SILVA classification.
    log:
        "logs/02-denoise-and-export-prok/11-unmerged-SILVA.log" #Declares the SILVA classification log path.
    script:
        "../scripts/P05-classify-eASVs.sh" #Runs QIIME 2's sklearn taxonomy classifier.

#This rule classifies the same retained-run features with PR2 so plastid lineages can be evaluated in the supplemental analysis.
rule classify_unmerged_16S_with_PR2: #Defines PR2 classification for the supplemental retained-pair analysis.
    input:
        sequences=rules.denoise_prok_dada2_with_unmerged.output.prokrepseqs, #Uses the representative features from the retained-pair DADA2 run.
        database_ready=rules.ensure_pr2_classifier.output.marker, #Requires the PR2 classifier to be validated before classification begins.
    params:
        classDB=PR2_CLASSIFIER, #Passes the resolved PR2 classifier path to q2-feature-classifier.
    output:
        classified="results/02-proks/11-unmerged-16S-rescue/PR2.classified.qza" #Stores PR2 taxonomy and confidence for the retained-run features.
    conda:
        config["qiime2version"]
    # Each sklearn worker loads its own classifier state. Two workers keep the large unmerged-ASV classification below the CARC job memory limit.
    threads: 2 #Uses two workers to limit the classifier's combined memory footprint.
    resources:
        mem_mb=110000, #Requests 110,000 MB of memory for PR2 classification.
    log:
        "logs/02-denoise-and-export-prok/11-unmerged-PR2.log" #Declares the PR2 classification log path.
    script:
        "../scripts/P05-classify-eASVs.sh" #Runs QIIME 2's sklearn taxonomy classifier.


#This rule exports the retained-pair DADA2 feature table as BIOM for comparison with the canonical merged-only table.
rule export_unmerged_16S_retained_table: #Defines the retained-run feature-table export.
    input:
        rules.denoise_prok_dada2_with_unmerged.output.proktable #Uses the abundance table from the DADA2 run that retained unmerged pairs.
    output:
        temp(".snakemake/classification-exports/unmerged-16S/retained-table/feature-table.biom") #Stores the exported retained-run BIOM table temporarily.
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/retained-table" #Defines the project-local export directory.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-retained-table.log" #Records messages and errors from the QIIME 2 export.
    shell:
        "mkdir -p {params.export_dir:q} && " #Creates the project-local export directory if needed.
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} " #Exports the QIIME 2 feature table as BIOM.
        "> {log:q} 2>&1" #Redirects command output and errors to the log.


#This rule exports the canonical merged-only DADA2 feature table so features unique to the retained-pair run can be identified.
rule export_unmerged_16S_standard_table: #Defines the canonical feature-table export used by the supplemental comparison.
    input:
        rules.denoise_prok_dada2.output.proktable #Uses the canonical merged-only 16S abundance table.
    output:
        temp(".snakemake/classification-exports/unmerged-16S/standard-table/feature-table.biom") #Stores the exported canonical BIOM table temporarily.
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/standard-table" #Defines the project-local export directory.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-standard-table.log" #Records messages and errors from the QIIME 2 export.
    shell:
        "mkdir -p {params.export_dir:q} && " #Creates the project-local export directory if needed.
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} " #Exports the canonical QIIME 2 table as BIOM.
        "> {log:q} 2>&1" #Redirects command output and errors to the log.


#This rule exports the representative features from the retained-pair run in QIIME 2's linked-sequence FASTA format.
rule export_unmerged_16S_sequences: #Defines the retained-run representative-sequence export.
    input:
        rules.denoise_prok_dada2_with_unmerged.output.prokrepseqs #Uses the representative features from the retained-pair run.
    output:
        temp(".snakemake/classification-exports/unmerged-16S/sequences/linked-dna-sequences.fasta") #Stores the linked forward/reverse representations temporarily in FASTA format.
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/sequences" #Defines the project-local sequence-export directory.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-sequences.log" #Records messages and errors from the sequence export.
    shell:
        "mkdir -p {params.export_dir:q} && " #Creates the project-local export directory if needed.
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} " #Exports the linked representative features as FASTA.
        "> {log:q} 2>&1" #Redirects command output and errors to the log.


#This rule exports the SILVA assignments for retained-run features as a temporary TSV used by the supplemental resolver.
rule export_unmerged_16S_SILVA_taxonomy: #Defines the supplemental SILVA taxonomy export.
    input:
        rules.classify_unmerged_16S_with_SILVA.output.classified #Uses the SILVA classification artifact for the retained-run features.
    output:
        temp(".snakemake/classification-exports/unmerged-16S/SILVA/taxonomy.tsv") #Stores feature IDs, taxonomy and confidence temporarily as TSV.
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/SILVA" #Defines the project-local SILVA export directory.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-SILVA-taxonomy.log" #Records messages and errors from the SILVA taxonomy export.
    shell:
        "mkdir -p {params.export_dir:q} && " #Creates the project-local export directory if needed.
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} " #Exports the SILVA classification as taxonomy TSV.
        "> {log:q} 2>&1" #Redirects command output and errors to the log.


#This rule exports the raw PR2 assignments for retained-run features so rank prefixes can be normalized without losing confidence values.
rule export_unmerged_16S_PR2_taxonomy: #Defines the supplemental raw PR2 taxonomy export.
    input:
        rules.classify_unmerged_16S_with_PR2.output.classified #Uses the PR2 classification artifact for the retained-run features.
    output:
        temp(".snakemake/classification-exports/unmerged-16S/PR2-raw/taxonomy.tsv") #Stores the unmodified PR2 taxonomy table temporarily.
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/PR2-raw" #Defines the project-local raw PR2 export directory.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-PR2-taxonomy.log" #Records messages and errors from the PR2 taxonomy export.
    shell:
        "mkdir -p {params.export_dir:q} && " #Creates the project-local export directory if needed.
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} " #Exports the raw PR2 classification as taxonomy TSV.
        "> {log:q} 2>&1" #Redirects command output and errors to the log.


#This rule removes PR2 rank prefixes from the supplemental linked-read taxonomy while preserving feature IDs, confidence values and all other columns.
rule normalize_unmerged_16S_PR2_rank_prefixes:
    input:
        taxonomy=rules.export_unmerged_16S_PR2_taxonomy.output[0] #Uses the raw exported PR2 classification for the linked-read sensitivity run.
    output:
        normalized=temp(".snakemake/classification-exports/unmerged-16S/PR2/taxonomy.tsv") #Stores the prefix-normalized taxonomy temporarily for the supplemental export.
    log:
        "logs/02-denoise-and-export-prok/11-normalize-PR2-rank-prefixes.log"
    #The helper changes only the Taxon column because RESCRIPt's edit-taxonomy action would discard the Confidence column.
    script:
        "../scripts/normalize_pr2_taxonomy.py" #Removes recognized PR2 rank markers without parsing or rewriting the remaining TSV fields.


#This rule converts the experimental linked-read DADA2 feature table from BIOM to TSV for comparison with the canonical table.
rule convert_unmerged_16S_retained_biom:
    input:
        biom=rules.export_unmerged_16S_retained_table.output #Uses the table from the DADA2 run that retains unmerged linked pairs.
    output:
        temp(".snakemake/classification-exports/unmerged-16S/retained-table.tsv") #Stores a temporary feature-by-sample count table.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-convert-retained-biom.log"
    shell:
        "biom convert -i {input.biom:q} -o {output:q} --to-tsv " #Converts the shell-quoted retained BIOM table to tab-separated text.
        "> {log:q} 2>&1" #Records command output and errors in the rule log.


#This rule converts the canonical merged-read DADA2 feature table from BIOM to TSV for comparison with the linked-read sensitivity run.
rule convert_unmerged_16S_standard_biom:
    input:
        biom=rules.export_unmerged_16S_standard_table.output #Uses the feature table from the standard merged-read DADA2 run.
    output:
        temp(".snakemake/classification-exports/unmerged-16S/standard-table.tsv") #Stores a temporary feature-by-sample count table.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-convert-standard-biom.log"
    shell:
        "biom convert -i {input.biom:q} -o {output:q} --to-tsv " #Converts the shell-quoted canonical BIOM table to tab-separated text.
        "> {log:q} 2>&1" #Records command output and errors in the rule log.


#This rule exports linked 16S features found only by the unmerged-read sensitivity run, together with sequences, abundance counts, reconciled taxonomy and recovery statistics.
rule export_unmerged_16S_supplement:
    input:
        retained_table_tsv=rules.convert_unmerged_16S_retained_biom.output[0], #Provides counts from the run that retains unmerged linked pairs.
        standard_table_tsv=rules.convert_unmerged_16S_standard_biom.output[0], #Provides canonical counts used to identify features unique to the sensitivity run.
        sequence_fasta=rules.export_unmerged_16S_sequences.output[0], #Provides the linked sequence representation for each retained feature.
        silva_export=rules.export_unmerged_16S_SILVA_taxonomy.output[0], #Provides the SILVA classifications and confidence values.
        pr2_export=rules.normalize_unmerged_16S_PR2_rank_prefixes.output.normalized, #Provides normalized PR2 classifications for plastid reconciliation.
    output:
        table=UNMERGED_16S_TABLE, #Writes the separate per-feature, per-sample linked-read supplement.
        summary=UNMERGED_16S_SUMMARY, #Writes overall recovery and supplemental chloroplast counts.
    params:
        min_confidence=CHLOROPLAST_PR2_MIN_CONFIDENCE, #Applies the configured minimum confidence before accepting a PR2 plastid assignment.
    script:
        "../scripts/export_unmerged_16S.py" #Finds features absent from the canonical run and builds the supplemental table and summary.

#This rule exports the canonical 16S DADA2 artifacts as FASTA, TSV and BIOM files and maintains stable project-named copies of the latest sequences and statistics.
rule export_DADA2_results:
    input:
        "results/02-proks/03-DADA2d/" #Uses the directory containing the canonical DADA2 artifacts.
    params:
        projectName=config["projectName"] #Provides the project prefix for stable and timestamped export filenames.
    output:
        directory("results/02-proks/04-DADA2d-plaintext-exports/"), #Stores the complete set of plaintext and BIOM exports.
        lateststats="results/02-proks/04-DADA2d-plaintext-exports/" + config["projectName"] + ".16S.latest_stats.tsv", #Provides the current denoising statistics under a stable project-named path.
        latestseqs="results/02-proks/04-DADA2d-plaintext-exports/" + config["projectName"] + ".16S.latest_seqs.fasta" #Provides the current representative ASV sequences under a stable project-named path.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/04_export_DADA2_results/DADA2_export.stderrout"
    script:
        "../scripts/P04-export-DADA2-results.sh" #Exports QIIME 2 artifacts, converts the feature table to TSV and archives timestamped copies.

#This rule classifies the canonical 16S ASVs with the resolved SILVA classifier.
rule classify_ASVs:
    input:
        sequences="results/02-proks/03-DADA2d/representative_sequences.qza", #Uses the canonical DADA2 representative ASV sequences.
        database_ready=rules.ensure_silva_classifier.output.marker, #Requires the selected SILVA classifier to be validated and ready.
    params:
        classDB=SILVA_CLASSIFIER, #Passes the resolved SILVA classifier path to q2-feature-classifier.
    output:
        directory("results/02-proks/05-classified/"), #Declares the classification output directory.
        classified="results/02-proks/05-classified/" + config["projectName"] + "_SILVA.classified.qza" #Stores the SILVA taxonomy assignment for every canonical ASV.
    conda:
        config["qiime2version"]
    threads: 8 #Allows the classifier to use eight parallel workers.
    resources:
        mem_mb=64000, #Requests 64,000 MB of memory for SILVA classification.
    log:
        "logs/02-denoise-and-export-prok/05-classify-ASVs.log"
    script:
        "../scripts/P05-classify-eASVs.sh" #Runs QIIME 2's scikit-learn taxonomy classifier.

#This rule combines sample metadata with canonical DADA2 retention statistics and upstream prokaryotic/eukaryotic split counts for downstream visualizations.
rule create_sample_metadata_file:
    input:
        manifest="results/02-proks/manifest.tsv", #Provides the imported sample IDs and paired-read paths.
        samplesdottsv="config/samples.tsv", #Provides the user-supplied sample metadata.
        prokstats=rules.export_DADA2_results.output.lateststats, #Provides each sample's non-chimeric 16S read count.
        eukfracpersample="results/" + config["projectName"] + ".eukfrac-per-sample.tsv" #Provides the upstream prokaryotic/eukaryotic counts and fractions.
    output:
        "results/02-proks/sample-metadata.tsv" #Stores the merged QIIME 2-compatible metadata table.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/06-make-sample-metadata-file.log"
    script:
        "../scripts/P06-make-sample-metadata-file.py" #Normalizes sample IDs, removes manifest paths and joins the four input tables by sample ID.

#This rule creates an interactive SILVA taxonomic-composition bar plot for the canonical 16S feature table.
rule make_SILVA_only_prok_barplots:
    input:
        proktable="results/02-proks/03-DADA2d/table.qza", #Uses the canonical sample-by-ASV abundance table.
        proktax=rules.classify_ASVs.output.classified, #Uses the SILVA taxonomy assignments.
        prokmetadata="results/02-proks/sample-metadata.tsv" #Uses the combined sample metadata and processing statistics.
    params:
        projectName=config["projectName"] #Provides the project name for the timestamped visualization filename.
    output:
        directory("results/02-proks/07-SILVA-only-barplots/") #Stores the QIIME 2 taxonomic bar-plot visualization.
    conda:
        config["qiime2version"]
    script:
        "../scripts/P07-make-barplot.sh" #Runs `qiime taxa barplot` and renames the visualization with its timestamp and project name.

#This rule classifies every canonical 16S ASV with PR2 at an unrestricted confidence threshold so the resolver can evaluate all possible plastid matches.
rule classify_all_16S_with_PR2:
    input:
        database_ready=rules.ensure_pr2_classifier.output.marker, #Requires the selected PR2 classifier to be validated and ready.
        prokseqs=rules.denoise_prok_dada2.output.prokrepseqs, #Uses the canonical DADA2 representative 16S sequences.
    params:
        PR2classifier=PR2_CLASSIFIER, #Passes the resolved PR2 classifier path to q2-feature-classifier.
    output:
        classified="results/02-proks/09-subsetting/reclassified/all_16S_ASVs_PR2.classified.qza", #Stores PR2 assignments and confidence values for all canonical ASVs.
    #Two workers limit simultaneous classifier copies because this PR2 classification can approach the CARC job memory limit.
    threads: 2 #Uses two classifier workers to control memory consumption.
    resources:
        mem_mb=110000, #Requests 110,000 MB of memory for PR2 classification.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/09-classify-all-16S-with-PR2.log"
    script:
        "../scripts/P09a-classify-all-16S-with-PR2.sh" #Runs `classify-sklearn` at confidence 0.0 so the resolver can apply the configured acceptance threshold itself.


#This rule exports the canonical SILVA taxonomy artifact as TSV for plastid reconciliation.
rule export_SILVA_taxonomy_for_plastid_resolution:
    input:
        rules.classify_ASVs.output.classified #Uses the SILVA classification of every canonical 16S ASV.
    output:
        temp(".snakemake/classification-exports/canonical-16S/SILVA/taxonomy.tsv") #Stores the exported SILVA taxonomy temporarily.
    params:
        export_dir=".snakemake/classification-exports/canonical-16S/SILVA" #Defines the project-local export directory expected to contain taxonomy.tsv.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/09-export-SILVA-taxonomy.log"
    shell:
        "mkdir -p {params.export_dir:q} && " #Creates the project-local export directory if necessary.
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} " #Exports the shell-quoted SILVA taxonomy artifact.
        "> {log:q} 2>&1" #Records command output and errors in the rule log.


#This rule exports the all-ASV PR2 classification as TSV so its rank prefixes can be normalized before plastid reconciliation.
rule export_PR2_taxonomy_for_plastid_resolution:
    input:
        rules.classify_all_16S_with_PR2.output.classified #Uses the PR2 classification of every canonical 16S ASV.
    output:
        temp(".snakemake/classification-exports/canonical-16S/PR2-raw/taxonomy.tsv") #Stores the raw exported PR2 taxonomy temporarily.
    params:
        export_dir=".snakemake/classification-exports/canonical-16S/PR2-raw" #Defines the project-local export directory expected to contain taxonomy.tsv.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/09-export-PR2-taxonomy.log"
    shell:
        "mkdir -p {params.export_dir:q} && " #Creates the project-local export directory if necessary.
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} " #Exports the shell-quoted PR2 taxonomy artifact.
        "> {log:q} 2>&1" #Records command output and errors in the rule log.


#This rule removes PR2 rank prefixes while preserving feature IDs and classifier confidence for canonical plastid reconciliation.
rule normalize_canonical_16S_PR2_rank_prefixes:
    input:
        taxonomy=rules.export_PR2_taxonomy_for_plastid_resolution.output[0] #Uses the raw exported all-ASV PR2 taxonomy table.
    output:
        normalized=temp(".snakemake/classification-exports/canonical-16S/PR2/taxonomy.tsv") #Stores the prefix-normalized table temporarily for reconciliation.
    log:
        "logs/02-denoise-and-export-prok/09-normalize-PR2-rank-prefixes.log"
    #Only prefixed PR2 rank notation changes; feature IDs, confidence values and all other columns remain intact.
    script:
        "../scripts/normalize_pr2_taxonomy.py" #Normalizes the Taxon column without applying the later plastid-confidence threshold.


#This rule conservatively replaces SILVA assignments with PR2 plastid lineages that meet the configured confidence threshold and records every routing decision.
rule resolve_PR2_plastids:
    input:
        silva_taxonomy=rules.export_SILVA_taxonomy_for_plastid_resolution.output[0], #Provides the baseline SILVA taxonomy and confidence for each canonical ASV.
        pr2_taxonomy=rules.normalize_canonical_16S_PR2_rank_prefixes.output.normalized, #Provides normalized PR2 candidates and confidence values for the same ASVs.
    output:
        taxonomy_tsv=temp("results/02-proks/09-subsetting/tax-merged/" + config["projectName"] + ".PR2-plastid-resolved-taxonomy.tsv"), #Stores the final import-ready taxonomy table temporarily.
        audit="results/02-proks/09-subsetting/tax-merged/" + config["projectName"] + ".PR2-plastid-routing-audit.tsv", #Records both classifications, the decision and the final assignment for every ASV.
        summary="results/02-proks/09-subsetting/tax-merged/" + config["projectName"] + ".PR2-plastid-routing-summary.tsv", #Summarizes confirmed, rescued, rejected and retained classifications.
    params:
        min_confidence=CHLOROPLAST_PR2_MIN_CONFIDENCE, #Accepts only PR2 lineages containing `:plas` at or above this threshold.
    script:
        "../scripts/resolve_pr2_plastids.py" #Keeps SILVA by default and substitutes only accepted PR2 plastid assignments.


#This rule imports the reconciled SILVA/PR2 taxonomy table into QIIME 2 for downstream table splitting, plots and exports.
rule import_resolved_PR2_plastid_taxonomy:
    input:
        rules.resolve_PR2_plastids.output.taxonomy_tsv #Uses the final reconciled taxonomy TSV.
    output:
        taxonomy="results/02-proks/09-subsetting/tax-merged/chloroplasts-PR2-reclassified-merged-classification.qza" #Stores the resolved taxonomy as a QIIME 2 `FeatureData[Taxonomy]` artifact.
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/09-import-resolved-taxonomy.log"
    shell:
        "qiime tools import --type 'FeatureData[Taxonomy]' " #Selects the QIIME 2 taxonomy artifact type.
        "--input-format TSVTaxonomyFormat --input-path {input:q} " #Reads the shell-quoted reconciled TSV using QIIME 2's taxonomy format.
        "--output-path {output.taxonomy:q} > {log:q} 2>&1" #Writes the artifact and records both command streams in the log.


# This rule divides the complete 16S feature table into organelle, cyanobacterial and archaeal categories using the resolved SILVA/PR2 taxonomy.
rule split_resolved_prok_categories:
    input:
        proktable=rules.denoise_prok_dada2.output.proktable, # Uses the complete denoised 16S ASV table.
        resolved_taxonomy=rules.import_resolved_PR2_plastid_taxonomy.output.taxonomy, # Uses taxonomy in which SILVA plastid assignments have been resolved with PR2 where possible.
    output:
        includechlorotable="results/02-proks/09-subsetting/split-tables/include_o__Chloroplast_filtered_table.qza", # Retains PR2 plastids and legacy SILVA chloroplast assignments.
        excludechlorotable="results/02-proks/09-subsetting/split-tables/exclude_o__Chloroplast_filtered_table.qza", # Removes PR2 plastids and legacy SILVA chloroplast assignments.
        onlymitotable="results/02-proks/09-subsetting/split-tables/include_f__Mitochondria_filtered_table.qza", # Retains only mitochondrial assignments.
        onlyalgaetable="results/02-proks/09-subsetting/split-tables/include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.qza", # Retains plastids plus Cyanobacteria and Cyanobacteriota.
        onlycyanotable="results/02-proks/09-subsetting/split-tables/include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.qza", # Retains Cyanobacteria and Cyanobacteriota while excluding plastids.
        nomitotable="results/02-proks/09-subsetting/split-tables/exclude_f__Mitochondria_filtered_table.qza", # Removes mitochondrial assignments.
        nomitonochlorotable="results/02-proks/09-subsetting/split-tables/exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.qza", # Removes mitochondrial and plastid assignments.
        nomitonochloronocyanotable="results/02-proks/09-subsetting/split-tables/exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.qza", # Removes mitochondria, plastids, Cyanobacteria and Cyanobacteriota.
        onlyarchaeatable="results/02-proks/09-subsetting/split-tables/include_d__Archaea_filtered_table.qza", # Retains only archaeal assignments.
        noarchaeatable="results/02-proks/09-subsetting/split-tables/exclude_d__Archaea_filtered_table.qza" # Removes archaeal assignments.
    conda:
        config["qiime2version"]
    script:
        "../scripts/P09b-split-resolved-prok-categories.sh" # Applies the category-specific QIIME 2 taxonomy filters.

# This rule exports the resolved taxonomy and converts the complete and category-specific QIIME 2 feature tables to BIOM format.
rule export_tax_convert_biom:
    input:
        mergedtax=rules.import_resolved_PR2_plastid_taxonomy.output.taxonomy, # Uses the resolved taxonomy for the plaintext taxonomy export.
        all16Stable="results/02-proks/03-DADA2d/table.qza", # Uses the complete 16S feature table.
        noarch="results/02-proks/09-subsetting/split-tables/exclude_d__Archaea_filtered_table.qza", # Uses the Archaea-excluded table.
        nomito="results/02-proks/09-subsetting/split-tables/exclude_f__Mitochondria_filtered_table.qza", # Uses the mitochondria-excluded table.
        nochloronomito="results/02-proks/09-subsetting/split-tables/exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.qza", # Uses the table excluding plastids and mitochondria.
        nochloro="results/02-proks/09-subsetting/split-tables/exclude_o__Chloroplast_filtered_table.qza", # Uses the plastid-excluded table.
        nochloronocyanonomito="results/02-proks/09-subsetting/split-tables/exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.qza", # Uses the table excluding plastids, mitochondria and cyanobacteria.
        onlyarch="results/02-proks/09-subsetting/split-tables/include_d__Archaea_filtered_table.qza", # Uses the Archaea-only table.
        onlymito="results/02-proks/09-subsetting/split-tables/include_f__Mitochondria_filtered_table.qza", # Uses the mitochondria-only table.
        onlychloro="results/02-proks/09-subsetting/split-tables/include_o__Chloroplast_filtered_table.qza", # Uses the plastid-only table.
        onlycyano="results/02-proks/09-subsetting/split-tables/include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.qza", # Uses the cyanobacteria-only table with plastids excluded.
        onlyalgae="results/02-proks/09-subsetting/split-tables/include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.qza", # Uses the combined cyanobacteria and plastid table.
    output:
        genus_mergedtaxtsv=temp("results/02-proks/10-exports/" + config["projectName"] + ".taxonomy-through-genus.tsv"), # Stores taxonomy with BIOM-compatible column names for later species assignment.
        all16Stable_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".all-16S-seqs.biom"), # Stores the complete feature table in BIOM format.
        noarch_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_d__Archaea_filtered_table.biom"), # Stores the Archaea-excluded BIOM table.
        nomito_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_f__Mitochondria_filtered_table.biom"), # Stores the mitochondria-excluded BIOM table.
        nochloronomito_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.biom"), # Stores the plastid- and mitochondria-excluded BIOM table.
        nochloro_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_filtered_table.biom"), # Stores the plastid-excluded BIOM table.
        nochloronocyanonomito_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.biom"), # Stores the plastid-, mitochondria- and cyanobacteria-excluded BIOM table.
        onlyarch_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_d__Archaea_filtered_table.biom"), # Stores the Archaea-only BIOM table.
        onlymito_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_f__Mitochondria_filtered_table.biom"), # Stores the mitochondria-only BIOM table.
        onlychloro_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_o__Chloroplast_filtered_table.biom"), # Stores the plastid-only BIOM table.
        onlycyano_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.biom"), # Stores the cyanobacteria-only BIOM table with plastids excluded.
        onlyalgae_biom=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.qza.biom"), # Stores the combined cyanobacteria and plastid BIOM table.
    conda:
        config["qiime2version"]
    script:
        "../scripts/P10a-generate-biom-tables.sh" # Exports the taxonomy and converts every QIIME 2 feature table to BIOM.

# This rule adds unambiguous exact-match SILVA 144 species labels when they agree with the classified genus, or passes older-release taxonomy through unchanged.
rule assign_SILVA_species:
    input:
        taxonomy=rules.export_tax_convert_biom.output.genus_mergedtaxtsv, # Uses the resolved taxonomy table through genus level.
        sequences=rules.export_DADA2_results.output.latestseqs, # Uses the representative 16S ASV sequences for exact matching.
        # The exact-match reference belongs to release 144 only. Older
        # classifiers pass their taxonomy through without loading this file.
        species_reference=([SILVA_SPECIES_REFERENCE] if RUN_SILVA_SPECIES_MATCHING else []), # Requires the SILVA 144 species reference only when exact matching is enabled.
    output:
        taxonomy="results/02-proks/10-exports/" + config["projectName"] + ".taxonomy.tsv", # Stores the final taxonomy used in downstream BIOM exports.
        summary="results/02-proks/10-exports/" + config["projectName"] + ".SILVA-" + SILVA_VERSION + ".species-assignment-summary.tsv", # Reports matching status and species-assignment counts.
    params:
        chunk_size=500, # Processes at most 500 ASVs per assignSpecies batch to limit peak memory.
        enabled=RUN_SILVA_SPECIES_MATCHING, # Enables exact species matching only for the supported SILVA release.
        silva_version=SILVA_VERSION, # Records the configured SILVA release in the summary.
    conda:
        config["qiime2version"]
    threads: 1 # Runs species matching in a single process.
    resources:
        # Reserve the node for this step. assignSpecies is also batched below,
        # so its actual peak memory no longer grows with the full ASV count.
        mem_mb=110000, # Reserves 110 GB while the reference is queried in batches.
        runtime=360, # Requests up to 360 minutes.
    log:
        "logs/02-denoise-and-export-prok/10-assign-SILVA-species.log"
    script:
        "../scripts/assign_silva_species.R" # Stages the reference locally, performs exact matching and appends genus-consistent species labels.

# This rule attaches the final taxonomy as semicolon-separated observation metadata to the complete and category-specific BIOM tables.
rule add_tax_to_biom:
    input:
        mergedtaxtsv=rules.assign_SILVA_species.output.taxonomy, # Uses taxonomy after optional SILVA 144 species assignment.
        all16Stable_biom="results/02-proks/10-exports/" + config["projectName"] + ".all-16S-seqs.biom",
        noarch_biom="results/02-proks/10-exports/" + config["projectName"] + ".exclude_d__Archaea_filtered_table.biom",
        nomito_biom="results/02-proks/10-exports/" + config["projectName"] + ".exclude_f__Mitochondria_filtered_table.biom",
        nochloronomito_biom="results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.biom",
        nochloro_biom="results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_filtered_table.biom",
        nochloronocyanonomito_biom="results/02-proks/10-exports/" + config["projectName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.biom",
        onlyarch_biom="results/02-proks/10-exports/" + config["projectName"] + ".include_d__Archaea_filtered_table.biom",
        onlymito_biom="results/02-proks/10-exports/" + config["projectName"] + ".include_f__Mitochondria_filtered_table.biom",
        onlychloro_biom="results/02-proks/10-exports/" + config["projectName"] + ".include_o__Chloroplast_filtered_table.biom",
        onlycyano_biom="results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.biom",
        onlyalgae_biom="results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.qza.biom",
    output:
        all16Stable_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".all-16S-seqs.with-tax.biom"),
        noarch_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_d__Archaea_filtered_table.with-tax.biom"),
        nomito_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_f__Mitochondria_filtered_table.with-tax.biom"),
        nochloronomito_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.with-tax.biom"),
        nochloro_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_filtered_table.with-tax.biom"),
        nochloronocyanonomito_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.with-tax.biom"),
        onlyarch_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_d__Archaea_filtered_table.with-tax.biom"),
        onlymito_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_f__Mitochondria_filtered_table.with-tax.biom"),
        onlychloro_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_o__Chloroplast_filtered_table.with-tax.biom"),
        onlycyano_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.with-tax.biom"),
        onlyalgae_biomtax=temp("results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.with-tax.biom")
    conda:
        config["qiime2version"]
    script:
        "../scripts/P10b-add-tax-to-biom.sh" # Runs `biom add-metadata` on the complete table and every taxonomic subset.

# This rule converts the taxonomy-enriched complete and category-specific BIOM tables into persistent tab-separated abundance tables.
rule export_biom_tsv:
    input:
        all16Stable_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".all-16S-seqs.with-tax.biom",
        noarch_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".exclude_d__Archaea_filtered_table.with-tax.biom",
        nomito_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".exclude_f__Mitochondria_filtered_table.with-tax.biom",
        nochloronomito_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.with-tax.biom",
        nochloro_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_filtered_table.with-tax.biom",
        nochloronocyanonomito_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.with-tax.biom",
        onlyarch_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".include_d__Archaea_filtered_table.with-tax.biom",
        onlymito_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".include_f__Mitochondria_filtered_table.with-tax.biom",
        onlychloro_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".include_o__Chloroplast_filtered_table.with-tax.biom",
        onlycyano_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.with-tax.biom",
        onlyalgae_biomtax="results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.with-tax.biom",
    output:
        all16Stable_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".all-16S-seqs.with-tax.tsv",
        noarch_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".exclude_d__Archaea_filtered_table.with-tax.tsv",
        nomito_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".exclude_f__Mitochondria_filtered_table.with-tax.tsv",
        nochloronomito_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.with-tax.tsv",
        nochloro_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".exclude_o__Chloroplast_filtered_table.with-tax.tsv",
        nochloronocyanonomito_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.with-tax.tsv",
        onlyarch_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".include_d__Archaea_filtered_table.with-tax.tsv",
        onlymito_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".include_f__Mitochondria_filtered_table.with-tax.tsv",
        onlychloro_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".include_o__Chloroplast_filtered_table.with-tax.tsv",
        onlycyano_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.with-tax.tsv",
        onlyalgae_biomtaxtsv="results/02-proks/10-exports/" + config["projectName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.with-tax.tsv",
    conda:
        config["qiime2version"]
    script:
        "../scripts/P10c-export-tax-tsvs.sh" # Converts every BIOM table to TSV with taxonomy included as observation metadata.

# This rule extracts Synechococcales ASVs and classifies their exact full-length matches against the external ProPortal 16S reference.
rule proportal_classify:
    input:
        proktable="results/02-proks/10-exports/" + config["projectName"] + ".all-16S-seqs.with-tax.tsv", # Uses the full taxonomy-enriched table to identify Synechococcales ASVs.
        dnaseqs="results/02-proks/04-DADA2d-plaintext-exports/" + config["projectName"] + ".16S.latest_seqs.fasta" # Uses the representative sequences for the selected ASV IDs.
    output:
        taxonomy="results/02-proks/10-exports/" + config["projectName"] + ".Synechococcales.proportal-classified.tsv", # Stores the separate ProPortal classification results.
    conda:
        "../envs/proportal.yml"
    script:
        "../scripts/add-proportal-taxonomy.sh" # Clones ProPortal, performs 100% identity and coverage BLAST matching and formats the classifications.
