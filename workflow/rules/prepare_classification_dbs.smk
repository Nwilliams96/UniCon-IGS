#This smk prepares the SILVA and PR2 taxonomy classifiers that downstream rules use to classify ASVs. It:
	# Downloads the separate SILVA 144 exact-species reference when that analysis is enabled—even if the classifier itself already exists.
    # Keeps construction intermediates in project-local working directories.
    # Uses save_classifier.py to validate completed classifiers and safely save shared models using locking and atomic replacement.
    # Creates readiness markers that downstream classification rules require.
#In short, common.smk resolves the choices and paths; this file contains the rules that obtain, prepare, validate and save the classifiers.

#This function chooses which classifier file the readiness rule should use. It returns the final classifier path when the user selected an existing classifier or when a reusable completed file is available. Otherwise, it returns the project-local staging path, which the download or training rules must produce.
def classifier_source(record): # Defines a function named classifier_source. record contains the resolved settings for one database—SILVA or PR2.
    if record["mode"] == "existing" or database_config.reusable_classifier( # Starts checking two alternative conditions: Is the selected mode "existing" or does the helper consider the completed classifier reusable.
        record["classifier"], record["checksum"] # Passes these values to helper script.
    ):
        return record["classifier"] # If either condition is satisfied, returns the final classifier path and ends the function.
    return record["staged"] # Selects the project-local staging file that downloading or training must create when no reusable completed classifier exists.

#This defines a rule only when SILVA is configured for local training. The rule checks whether the RESCRIPt installed in the chosen QIIME environment supports downloading the requested SILVA release.
if SILVA_DB["mode"] == "train": # Includes the following rule only when SILVA’s selected mode is "train". 
    rule check_silva_reference_capability: # Define the snakemake rule.
        output:
            os.path.join(SILVA_BUILD_ROOT, "capability.json") # Constructs the output path inside the SILVA build directory.
        params:
            database="SILVA", # Tells the helper to inspect SILVA acquisition rather than PR2 acquisition.
            version=SILVA_VERSION # Passes the resolved SILVA release requested in the configuration.
        conda:
            config["qiime2version"] #Selects the configured QIIME environment
        script:
            "../scripts/check_reference_capability.py" # Runs the helper relative to this rule file. check_reference_capability.py is a preflight check for local database training. It checks whether the installed RESCRIPt supports acquiring the requested SILVA or PR2 release.

#This defines a capability-check rule only when PR2 is configured for local training. It checks whether RESCRIPt in the selected QIIME environment supports acquiring the requested PR2 release.
if PR2_DB["mode"] == "train": # Includes the following rule only when PR2’s selected mode is "train".
    rule check_pr2_reference_capability: # Defines snakemake rule.
        output:
            os.path.join(PR2_BUILD_ROOT, "capability.json") # Constructs the output path inside the PR2 build directory.
        params:
            database="PR2", # Selects RESCRIPt's PR2 acquisition action for the capability check.
            version=PR2_VERSION # Supplies the configured PR2 release whose availability must be checked.
        conda:
            config["qiime2version"] #Begins specifying the environment used to execute the check.
        script:
            "../scripts/check_reference_capability.py" # checks the permitted release values of RESCRIPt’s get_pr2_data action.

#This defines the SILVA classifier-download rule when “Download pre-trained database” is selected. It downloads an already-trained classifier into a project-local staging location, verifies its checksum when supplied, and checks that it is a valid QIIME taxonomy-classifier artifact.
if SILVA_DB["mode"] == "download": # Includes the following rule only when SILVA’s resolved mode is "download".
    rule download_SILVA_classifier: # Defines rule.
        output:
            temp(SILVA_DB["staged"]) # Uses SILVA’s staging path as the output.
        params:
            url=SILVA_DB["url"], # Passes the resolved classifier-download URL.
            md5=SILVA_DB["checksum"], # Passes the expected MD5 checksum, if available.
            validate_classifier=True # Requires the helper to check that the file is a TaxonomicClassifier artifact and run QIIME validation.
        conda:
            config["qiime2version"] #Selects the configured QIIME environment, providing the QIIME command needed for artifact validation.
        log:
            "logs/SILVA_classifier_download.log" # Specifies where download and validation messages are recorded.
        script:
            "../scripts/download_verified_reference.py" # Runs the verified downloader relative to this rule file.

#This defines the PR2 download rule when “Download pre-trained database” is selected. It downloads an already-trained classifier into a project-local staging location, checks its checksum when supplied, and validates its QIIME artifact type and integrity.
if PR2_DB["mode"] == "download": # Includes the following rule only when PR2’s resolved mode is "download".
    rule download_PR2_classifier: # Defines rule.
        output:
            temp(PR2_DB["staged"]) # Provides staging path.
        params:
            url=PR2_DB["url"], # Passes the resolved PR2 classifier-download URL.
            md5=PR2_DB["checksum"], # Passes the expected MD5 checksum, if supplied.
            validate_classifier=True # Requires the helper to confirm that the download is a TaxonomicClassifier artifact and run QIIME validation.
        conda:
            config["qiime2version"] #Selects the configured QIIME environment.
        log:
            "logs/PR2_classifier_download.log" # Specifies where download and validation messages are recorded.
        script:
            "../scripts/download_verified_reference.py" # Runs the verified downloader relative to this rule file.

#This defines a rule to download the separate SILVA 144 reference used for exact species matching. In the current workflow, that step is enabled when the selected SILVA release is "144".
if RUN_SILVA_SPECIES_MATCHING: # Includes the following rule only when SILVA species matching is enabled.
    rule download_SILVA_144_species_reference: # Defines rule.
        input:
            rules.initialize_database_directories.output.marker # Requires the database-directory initialization marker, ensuring initialization completes before downloading.
        output:
            SILVA_SPECIES_REFERENCE # Uses the full destination path defined in common.smk for silva_v144_assignSpecies.fa.gz. This is a gzip-compressed FASTA reference. It is not marked temp(), so Snakemake does not automatically delete it after use.
        params:
            url=SILVA_SPECIES_REFERENCE_URL, # Passes the SILVA 144 species-reference download URL defined in common.smk.
            md5=SILVA_SPECIES_REFERENCE_MD5 # Passes the expected MD5 checksum. The helper uses it to verify that the file’s contents match the intended reference.
        log:
            "logs/SILVA_144_species_reference_download.log" # Specifies where download messages are recorded.
        script:
            "../scripts/download_verified_reference.py" # Runs the verified downloader relative to this rule file. It handles locking, reuse of a checksum-matching file, downloading and safe saving.

# Defines the SILVA training outputs and taxonomy ranks when the reference database must be trained.
if SILVA_DB["mode"] == "train": # Runs this configuration only in training mode.
    SILVA_FILES = { # Begins mapping each intermediate dataset name to its QIIME 2 artifact path.
        name: os.path.join(SILVA_BUILD_ROOT, name + ".qza") # Builds the artifact path from the output directory and dataset name.
        for name in ("rna", "taxonomy", "dna", "culled", "filtered", "discarded", # Begins iterating over the required SILVA artifacts.
                     "derep", "derep-taxonomy", "sliced", "extraction-stats", # Includes dereplication, slicing and extraction-statistics artifacts.
                     "sliced-derep", "sliced-taxonomy") # Includes the final sliced sequence and taxonomy artifacts. "Sliced" is the trimmed in silico to region (simulated amplicon)
    } # Completes the artifact-path dictionary.
    SILVA_RANKS = ["domain"] + (["kingdom"] if SILVA_VERSION == "144" else []) + [ # Includes kingdom only for SILVA 144.
        "phylum", "class", "order", "family", "genus" # Adds the remaining supported taxonomic ranks.
    ] # Completes the ordered rank list.

 #This rule downloads the selected SILVA reference sequences and taxonomy through QIIME 2 RESCRIPt for subsequent classifier training.
    rule download_SILVA: #Defines the Snakemake rule that downloads the SILVA reference data.
        input:
            rules.check_silva_reference_capability.output # Requires successful confirmation that the selected environment can retrieve the requested SILVA release.
        output:
            seqs=temp(SILVA_FILES["rna"]), # Creates the SILVA RNA sequence artifact and marks it for automatic removal after its downstream use.
            taxonomy=temp(SILVA_FILES["taxonomy"]) # Creates the corresponding taxonomy artifact and marks it as temporary.
        params:
            version=SILVA_VERSION, # Passes the configured SILVA release number.
            ranks=SILVA_RANKS # Passes the ordered taxonomic ranks appropriate for that SILVA release.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        log:
            "logs/SILVA_classification_db_prep/01-download-reference.log" # Records the command’s standard output and error messages.
        shell:
            "qiime rescript get-silva-data --p-version {params.version:q} " # Starts RESCRIPt’s SILVA downloader and supplies the shell-quoted release number.
            "--p-target SSURef_NR99 --p-ranks {params.ranks:q} " # Requests the NR99 small-subunit reference and the configured taxonomy ranks.
            "--p-rank-propagation --p-no-include-species-labels " # Propagates higher-rank labels into missing ranks and excludes SILVA species labels.
            "--o-silva-sequences {output.seqs:q} " # Writes the downloaded sequences to the declared sequence artifact.
            "--o-silva-taxonomy {output.taxonomy:q} > {log:q} 2>&1" # Writes the taxonomy artifact and redirects both output streams to the log.
    
#This rule converts the downloaded SILVA RNA reference sequences into DNA sequences for downstream filtering and classifier preparation.
    rule reverse_transcribe: #Defines the Snakemake rule that converts the reference sequence alphabet from RNA to DNA.
        input:
            rules.download_SILVA.output.seqs # Uses the SILVA RNA sequence artifact produced by the download rule.
        output:
            temp(SILVA_FILES["dna"]) # Writes the DNA sequence artifact and marks it for automatic removal after downstream use.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        log:
            "logs/SILVA_classification_db_prep/02-reverse-transcribe.log" # Records messages and errors from the conversion command.
        shell:
            "qiime rescript reverse-transcribe --i-rna-sequences {input:q} " # Runs RESCRIPt’s reverse-transcription action on the shell-quoted RNA artifact.
            "--o-dna-sequences {output:q} > {log:q} 2>&1" # Writes the DNA artifact and redirects both output streams to the log.
       
            
#This rule removes SILVA reference sequences containing excessive ambiguous bases or long homopolymer runs before further database preparation.
    rule qc_seqs_cull: #Defines the Snakemake rule that screens the DNA reference sequences for quality problems.
        input:
            rules.reverse_transcribe.output # Uses the DNA sequences produced by the reverse-transcription rule.
        output:
            temp(SILVA_FILES["culled"]) # Writes the sequences that pass screening and marks the artifact as temporary.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        log:
            "logs/SILVA_classification_db_prep/03-cull.log" # Records messages and errors from sequence screening.
        shell:
            "qiime rescript cull-seqs --i-sequences {input:q} " # Runs RESCRIPt’s sequence-culling action on the shell-quoted input artifact.
            "--p-num-degenerates 5 --p-homopolymer-length 8 " # Removes sequences with at least five degenerate bases or a homopolymer of at least eight bases.
            "--o-clean-sequences {output:q} > {log:q} 2>&1" # Writes the retained sequences and redirects both output streams to the log.
            
            
    #This rule removes SILVA sequences below taxon-specific minimum lengths. Rejected sequences are temporary; use --notemp if they must remain available for inspection.
    rule qc_seqs_filter: #Defines the Snakemake rule that filters reference sequences by length and taxonomy.
        input:
            sequences=rules.qc_seqs_cull.output, # Uses the sequences that passed degenerate-base and homopolymer screening.
            taxonomy=rules.download_SILVA.output.taxonomy # Uses the downloaded SILVA taxonomy to apply group-specific thresholds.
        output:
            filteredDNA=temp(SILVA_FILES["filtered"]), # Stores sequences that meet their applicable minimum length and marks the artifact as temporary.
            discardedDNA=temp(SILVA_FILES["discarded"]) # Stores rejected sequences temporarily; Snakemake may remove this artifact unless temporary outputs are retained.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        log:
            "logs/SILVA_classification_db_prep/04-filter-length.log" # Records messages and errors from taxon-specific length filtering.
        shell:
            "qiime rescript filter-seqs-length-by-taxon " # Starts RESCRIPt’s taxon-aware sequence-length filtering action.
            "--i-sequences {input.sequences:q} --i-taxonomy {input.taxonomy:q} " # Supplies the shell-quoted sequence and taxonomy artifacts.
            "--p-labels Bacteria Archaea Eukaryota --p-min-lens 1200 900 1400 " # Requires minimum lengths of 1,200, 900 and 1,400 bases for Bacteria, Archaea and Eukaryota, respectively.
            "--o-filtered-seqs {output.filteredDNA:q} " # Writes sequences that meet the applicable minimum length.
            "--o-discarded-seqs {output.discardedDNA:q} > {log:q} 2>&1" # Writes rejected sequences and redirects both output streams to the log.
            
            
    #This rule collapses duplicate full-length SILVA references while retaining distinct taxonomic annotations associated with identical sequences.
    rule qc_seqs_dereplicate: #Defines the Snakemake rule that dereplicates the quality-filtered SILVA reference.
        input:
            sequences=rules.qc_seqs_filter.output.filteredDNA, # Uses the sequences that passed taxon-specific length filtering.
            taxonomy=rules.download_SILVA.output.taxonomy # Uses the downloaded SILVA taxonomic assignments associated with the sequences.
        output:
            sequences=temp(SILVA_FILES["derep"]), # Stores the dereplicated full-length sequences and marks the artifact as temporary.
            taxonomy=temp(SILVA_FILES["derep-taxonomy"]) # Stores the corresponding taxonomy and marks the artifact as temporary.
        params:
            ranks=SILVA_RANKS # Provides the ordered taxonomic ranks used to fill missing ranks in the output taxonomy.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        log:
            "logs/SILVA_classification_db_prep/05-dereplicate-full-reference.log" # Records messages and errors from full-reference dereplication.
        shell:
            "qiime rescript dereplicate --i-sequences {input.sequences:q} " # Starts RESCRIPt dereplication using the shell-quoted sequence artifact.
            "--i-taxa {input.taxonomy:q} --p-mode uniq --p-rank-handles {params.ranks:q} " # Supplies taxonomy, retains each unique taxonomic affiliation and backfills the configured ranks.
            "--o-dereplicated-sequences {output.sequences:q} " # Writes the dereplicated sequence artifact.
            "--o-dereplicated-taxa {output.taxonomy:q} > {log:q} 2>&1" # Writes the matched taxonomy artifact and redirects both output streams to the log.

#This section completes local SILVA classifier training, builds the selected PR2 classifier and validates both classifiers for downstream taxonomy assignment.
#This rule extracts the configured primer-bounded region from the dereplicated SILVA reference sequences.
    rule extract_primers: #Defines the Snakemake rule that creates simulated SILVA amplicons.
        input:
            rules.qc_seqs_dereplicate.output.sequences # Uses the dereplicated full-length SILVA sequences.
        output:
            slicedDNA=temp(SILVA_FILES["sliced"]), # Stores the extracted primer-bounded sequences and marks the artifact as temporary.
            extractionStats=temp(SILVA_FILES["extraction-stats"]) # Stores statistics describing primer-region extraction and marks them as temporary.
        params:
            fwdPrimer=config["fwdPrimer"], #Uses the configured forward-primer sequence.
            revPrimer=config["revPrimer"] #Uses the configured reverse-primer sequence.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing q2-feature-classifier.
        threads: 8 # Requests up to eight workers; Snakemake can reduce this to match available cores.
        log:
            "logs/SILVA_classification_db_prep/06-extract-primer-region.log" # Records messages and errors from SILVA primer-region extraction.
        shell:
            "qiime feature-classifier extract-reads --i-sequences {input:q} " # Runs q2-feature-classifier’s read-extraction action on the SILVA reference.
            "--p-f-primer {params.fwdPrimer:q} --p-r-primer {params.revPrimer:q} " # Supplies the shell-quoted forward and reverse primers.
            "--p-n-jobs {threads} --p-read-orientation forward " # Uses the declared workers and searches for the primers in the forward orientation.
            "--o-reads {output.slicedDNA:q} " # Writes the extracted simulated amplicons.
            "--o-read-extraction-stats {output.extractionStats:q} > {log:q} 2>&1" # Writes extraction statistics and redirects both output streams to the log.

    #RESCRIPt: This rule dereplicates the extracted SILVA amplicons while retaining distinct taxonomic affiliations.
    rule dereplicated_sliced_data: #Defines the Snakemake rule that dereplicates the primer-region reference.
        input:
            sequences=rules.extract_primers.output.slicedDNA, # Uses the simulated amplicons produced by primer extraction.
            taxonomy=rules.qc_seqs_dereplicate.output.taxonomy # Uses the taxonomy associated with the full-length dereplicated reference.
        output:
            sequences=temp(SILVA_FILES["sliced-derep"]), # Stores the dereplicated simulated amplicons and marks the artifact as temporary.
            taxonomy=temp(SILVA_FILES["sliced-taxonomy"]) # Stores their corresponding taxonomy and marks the artifact as temporary.
        params:
            ranks=SILVA_RANKS # Provides the SILVA ranks used to backfill missing taxonomic levels.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        log:
            "logs/SILVA_classification_db_prep/07-dereplicate-primer-region.log" # Records messages and errors from amplicon dereplication.
        shell:
            "qiime rescript dereplicate --i-sequences {input.sequences:q} " # Starts RESCRIPt dereplication using the extracted SILVA sequences.
            "--i-taxa {input.taxonomy:q} --p-mode uniq --p-rank-handles {params.ranks:q} " # Retains unique taxonomic affiliations and backfills the configured ranks.
            "--o-dereplicated-sequences {output.sequences:q} " # Writes the dereplicated amplicon sequences.
            "--o-dereplicated-taxa {output.taxonomy:q} > {log:q} 2>&1" # Writes their taxonomy and redirects both output streams to the log.

    #This rule trains the SILVA naïve Bayes classifier using either the amplicon-specific or full-length reference.
    rule train_classifier: #Defines the Snakemake rule that trains the SILVA classifier.
        input:
            sequences=(SILVA_FILES["sliced-derep"] if config["silva_reference_mode"] == "amplicon" else SILVA_FILES["derep"]), # Uses simulated amplicons in amplicon mode and full-length sequences otherwise.
            taxonomy=(SILVA_FILES["sliced-taxonomy"] if config["silva_reference_mode"] == "amplicon" else SILVA_FILES["derep-taxonomy"]) # Uses the taxonomy corresponding to the selected reference sequences.
        output:
            temp(SILVA_DB["staged"]) # Writes the classifier to its project-local staging path and marks it as temporary.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing q2-feature-classifier.
        resources:
            mem_mb=110000, # Declares 110,000 MB for resource-aware scheduling; does not increase Slurm's memory allocation.
            runtime=720 # Declares an estimated runtime of 720 minutes; does not extend or enforce Slurm's time limit.
        log:
            "logs/SILVA_classification_db_prep/08-train-classifier.log" # Records messages and errors from SILVA classifier training.
        shell:
            "qiime feature-classifier fit-classifier-naive-bayes " # Starts q2-feature-classifier’s naïve Bayes training action.
            "--i-reference-reads {input.sequences:q} " # Supplies the selected shell-quoted reference sequences.
            "--i-reference-taxonomy {input.taxonomy:q} " # Supplies the corresponding shell-quoted taxonomy.
            "--o-classifier {output:q} > {log:q} 2>&1" # Writes the trained classifier and redirects both output streams to the log.


#This block defines the files and rules needed only when the PR2 classifier must be trained locally.
if PR2_DB["mode"] == "train": # Includes the following PR2 preparation rules only in local-training mode.
    PR2_CLEAN_QZA = os.path.join(PR2_BUILD_ROOT, "sequences.qza") # Defines the path for the downloaded PR2 sequences.
    PR2_RAW_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "raw-taxonomy.qza") # Defines the path for the original RESCRIPt-formatted PR2 taxonomy.
    PR2_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "taxonomy.qza") # Defines the path for the normalized PR2 taxonomy.
    PR2_CULLED_QZA = os.path.join(PR2_BUILD_ROOT, "culled.qza") # Defines the path for sequences that pass quality screening.
    PR2_DEREP_QZA = os.path.join(PR2_BUILD_ROOT, "derep.qza") # Defines the path for the dereplicated full-length sequences.
    PR2_DEREP_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "derep-taxonomy.qza") # Defines the path for the full-length dereplicated taxonomy.
    PR2_SLICED_QZA = os.path.join(PR2_BUILD_ROOT, "sliced.qza") # Defines the path for the primer-extracted sequences.
    PR2_EXTRACTION_STATS_QZA = os.path.join(PR2_BUILD_ROOT, "extraction-stats.qza") # Defines the path for primer-extraction statistics.
    PR2_SLICED_DEREP_QZA = os.path.join(PR2_BUILD_ROOT, "sliced-derep.qza") # Defines the path for dereplicated simulated amplicons.
    PR2_SLICED_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "sliced-taxonomy.qza") # Defines the path for the simulated amplicons’ taxonomy.
    PR2_STAGED_CLASSIFIER = PR2_DB["staged"] # Uses the resolved project-local staging path for the trained PR2 classifier.

    #RESCRIPt: This rule downloads the selected PR2 sequences and nine-rank taxonomy for local classifier training.
    rule download_pr2_reference: # Defines the Snakemake rule that obtains the PR2 reference data.
        input:
            rules.check_pr2_reference_capability.output # Requires successful confirmation that RESCRIPt supports the selected PR2 release.
        output:
            sequences=temp(PR2_CLEAN_QZA), # Stores the downloaded PR2 sequences and marks the artifact as temporary.
            taxonomy=temp(PR2_RAW_TAXONOMY_QZA) # Stores the original PR2 taxonomy and marks the artifact as temporary.
        params:
            version=PR2_VERSION # Passes the configured PR2 release number.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        log:
            "logs/PR2_classification_db_prep/00-download-reference.log" # Records messages and errors from the PR2 download.
        shell:
            "qiime rescript get-pr2-data --p-version {params.version:q} " # Runs RESCRIPt’s PR2 downloader with the shell-quoted release number.
            "--p-ranks domain supergroup division subdivision class order family genus species " # Requests PR2’s complete ordered nine-rank taxonomy.
            "--o-pr2-sequences {output.sequences:q} " # Writes the downloaded PR2 sequences.
            "--o-pr2-taxonomy {output.taxonomy:q} > {log:q} 2>&1" # Writes the PR2 taxonomy and redirects both output streams to the log.

    #RESCRIPt: This rule removes RESCRIPt’s PR2 rank prefixes while preserving all nine ranks and :plas labels.
    rule normalize_pr2_reference_taxonomy: # Defines the Snakemake rule that normalizes PR2 taxonomy labels.
        input:
            rules.download_pr2_reference.output.taxonomy # Uses the original taxonomy produced by the PR2 download rule.
        output:
            temp(PR2_TAXONOMY_QZA) # Stores the normalized taxonomy and marks the artifact as temporary.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        log:
            "logs/PR2_classification_db_prep/01-normalize-taxonomy.log" # Records messages and errors from taxonomy normalization.
        shell:
            #The raw multiline command runs RESCRIPt’s taxonomy editor.
            #The regular expression matches PR2 rank prefixes at the start of a label or after a semicolon.
            #The replacement retains the matched start or semicolon while removing only the rank prefix.
            #The final command writes the normalized taxonomy and redirects both output streams to the log.
            r"""
            qiime rescript edit-taxonomy --i-taxonomy {input:q} \
              --p-search-strings '(^|;)(?:d__|sgr__|dv__|dvs__|c__|o__|f__|g__|s__)' \
              --p-replacement-strings '\1' --p-use-regex \
              --o-edited-taxonomy {output:q} > {log:q} 2>&1
            """ #Completes the raw multiline shell command.

    #RESCRIPt: This rule removes PR2 sequences containing excessive ambiguous bases or long homopolymer runs.
    rule cull_pr2_seqs: #Defines the Snakemake rule that screens PR2 sequences for quality problems.
        input:
            PR2_CLEAN_QZA # Uses the downloaded PR2 sequence artifact.
        output:
            temp(PR2_CULLED_QZA) # Stores sequences that pass screening and marks the artifact as temporary.
        log:
            "logs/PR2_classification_db_prep/05-cull-sequences.log" # Records messages and errors from PR2 sequence screening.
        priority: 46 #Assigns this rule Snakemake scheduling priority 46.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        shell:
            "qiime rescript cull-seqs --i-sequences {input:q} " # Runs RESCRIPt’s sequence-culling action on the PR2 reference.
            "--p-num-degenerates 5 --p-homopolymer-length 8 " # Removes sequences with at least five degenerate bases or a homopolymer of at least eight bases.
            "--o-clean-sequences {output:q} > {log:q} 2>&1" # Writes retained sequences and redirects both output streams to the log.

    #RESCRIPt: This rule dereplicates the full-length PR2 sequences while preserving PR2’s existing nine-rank taxonomy format.
    rule derep_seqs_taxonomy: #Defines the Snakemake rule that dereplicates the cleaned PR2 reference.
        input:
            culled=rules.cull_pr2_seqs.output, # Uses the PR2 sequences that passed quality screening.
            taxonomy=PR2_TAXONOMY_QZA # Uses the normalized nine-rank PR2 taxonomy.
        output:
            derepseqs=temp(PR2_DEREP_QZA), # Stores the dereplicated full-length sequences and marks the artifact as temporary.
            dereptaxa=temp(PR2_DEREP_TAXONOMY_QZA) # Stores their corresponding taxonomy and marks the artifact as temporary.
        log:
            "logs/PR2_classification_db_prep/06-dereplicate-full-reference.log" # Records messages and errors from full-length PR2 dereplication.
        priority: 45 #Assigns this rule Snakemake scheduling priority 45.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        shell:
            "qiime rescript dereplicate --i-sequences {input.culled:q} " # Starts RESCRIPt dereplication using the cleaned PR2 sequences.
            "--i-taxa {input.taxonomy:q} --p-mode uniq --p-rank-handles disable " # Retains unique taxonomic affiliations without applying SILVA-style rank backfilling.
            "--o-dereplicated-sequences {output.derepseqs:q} " # Writes the dereplicated full-length sequences.
            "--o-dereplicated-taxa {output.dereptaxa:q} > {log:q} 2>&1" # Writes their taxonomy and redirects both output streams to the log.

    #This rule extracts the configured primer-bounded region from the dereplicated PR2 sequences.
    rule extract_primers_pr2: #Defines the Snakemake rule that creates simulated PR2 amplicons.
        input:
            rules.derep_seqs_taxonomy.output.derepseqs # Uses the dereplicated full-length PR2 sequences.
        output:
            slicedDNA=temp(PR2_SLICED_QZA), # Stores the simulated PR2 amplicons and marks the artifact as temporary.
            extractionStats=temp(PR2_EXTRACTION_STATS_QZA) # Stores primer-extraction statistics and marks them as temporary.
        params:
            fwdPrimer=config["fwdPrimer"], #Uses the configured forward-primer sequence.
            revPrimer=config["revPrimer"] #Uses the configured reverse-primer sequence.
        log:
            "logs/PR2_classification_db_prep/07-extract-primer-region.log" # Records messages and errors from PR2 primer-region extraction.
        priority: 44 #Assigns this rule Snakemake scheduling priority 44.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing q2-feature-classifier.
        threads: 8 # Requests up to eight workers; Snakemake can reduce this to match available cores.
        shell:
            "qiime feature-classifier extract-reads --i-sequences {input:q} " # Runs q2-feature-classifier’s read-extraction action on the PR2 reference.
            "--p-f-primer {params.fwdPrimer:q} --p-r-primer {params.revPrimer:q} " # Supplies the shell-quoted forward and reverse primers.
            "--p-n-jobs {threads} --p-read-orientation forward " # Uses the declared workers and searches for primers in the forward orientation.
            "--o-reads {output.slicedDNA:q} " # Writes the extracted simulated amplicons.
            "--o-read-extraction-stats {output.extractionStats:q} > {log:q} 2>&1" # Writes extraction statistics and redirects both output streams to the log.

    #RESCRIPt: This rule dereplicates the extracted PR2 amplicons while retaining distinct taxonomic affiliations.
    rule dereplicate_extracted_pr2_reads: #Defines the Snakemake rule that dereplicates the PR2 primer-region reference.
        input:
            slicedDNA=rules.extract_primers_pr2.output.slicedDNA, #Uses the simulated PR2 amplicons produced by primer extraction.
            dereplicatedTaxa=rules.derep_seqs_taxonomy.output.dereptaxa # Uses the taxonomy associated with the full-length dereplicated reference.
        output:
            slicedDNAdereplicated=temp(PR2_SLICED_DEREP_QZA), # Stores the dereplicated simulated amplicons and marks the artifact as temporary.
            dereplicatedTaxaSliced=temp(PR2_SLICED_TAXONOMY_QZA) # Stores their corresponding taxonomy and marks the artifact as temporary.
        log:
            "logs/PR2_classification_db_prep/08-dereplicate-primer-region.log" # Records messages and errors from PR2 amplicon dereplication.
        priority: 43 #Assigns this rule Snakemake scheduling priority 43.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing RESCRIPt.
        shell:
            "qiime rescript dereplicate --i-sequences {input.slicedDNA:q} " # Starts RESCRIPt dereplication using the extracted PR2 sequences.
            "--i-taxa {input.dereplicatedTaxa:q} --p-mode uniq --p-rank-handles disable " # Retains unique affiliations without applying SILVA-style rank backfilling.
            "--o-dereplicated-sequences {output.slicedDNAdereplicated:q} " # Writes the dereplicated simulated amplicons.
            "--o-dereplicated-taxa {output.dereplicatedTaxaSliced:q} > {log:q} 2>&1" # Writes their taxonomy and redirects both output streams to the log.

    #This rule trains the PR2 naïve Bayes classifier using either the amplicon-specific or full-length reference.
    rule train_classifier_pr2: #Defines the Snakemake rule that trains the PR2 classifier.
        input:
            slicedDNAdereplicated=(PR2_SLICED_DEREP_QZA if config["pr2_reference_mode"] == "amplicon" else PR2_DEREP_QZA), # Uses simulated amplicons in amplicon mode and full-length sequences otherwise.
            dereplicatedTaxaSliced=(PR2_SLICED_TAXONOMY_QZA if config["pr2_reference_mode"] == "amplicon" else PR2_DEREP_TAXONOMY_QZA) # Uses the taxonomy corresponding to the selected sequences.
        output:
            temp(PR2_STAGED_CLASSIFIER) # Writes the classifier to its project-local staging path and marks it as temporary.
        log:
            "logs/PR2_classification_db_prep/09-train-classifier.log" # Records messages and errors from PR2 classifier training.
        priority: 42 # Assigns this rule Snakemake scheduling priority 42.
        conda:
            config["qiime2version"] #Uses the configured QIIME 2 environment containing q2-feature-classifier.
        threads: 1 # Declares one thread for classifier training.
        resources:
            mem_mb=32000, # Declares 32,000 MB for resource-aware scheduling; does not increase Slurm's memory allocation.
            runtime=360, # Declares an estimated runtime of 360 minutes; does not extend or enforce Slurm's time limit.
        shell:
            "qiime feature-classifier fit-classifier-naive-bayes " # Starts q2-feature-classifier’s naïve Bayes training action.
            "--i-reference-reads {input.slicedDNAdereplicated:q} " # Supplies the selected shell-quoted reference sequences.
            "--i-reference-taxonomy {input.dereplicatedTaxaSliced:q} " # Supplies the corresponding shell-quoted taxonomy.
            "--o-classifier {output:q} > {log:q} 2>&1" # Writes the trained classifier and redirects both output streams to the log.


#This rule validates and safely saves or reuses the resolved SILVA classifier.
rule ensure_silva_classifier: # Defines the readiness rule required by downstream SILVA classification.
    input:
        classifier=lambda wc: classifier_source(SILVA_DB) # Selects the existing classifier or the staged download or training result.
    output:
        marker=SILVA_READY_MARKER # Creates the marker used to confirm that the final SILVA classifier is ready.
    params:
        classifier=SILVA_CLASSIFIER, # Provides the final SILVA classifier path.
        lock=SILVA_DB["lock"], # Coordinates competing saves to the shared classifier destination.
        read_only=SILVA_DB["mode"] == "existing", # Prevents modification of a user-supplied existing classifier.
        md5=SILVA_DB["checksum"], # Provides the expected checksum when one is available.
        label="SILVA" # Provides the database label used in validation messages.
    conda:
        config["qiime2version"] #Uses the configured QIIME 2 environment for artifact validation.
    log:
        "logs/SILVA_classification_db_prep/validate-and-save.log" # Records classifier validation, reuse and saving messages.
    script:
        "../scripts/save_classifier.py" # Validates the classifier and safely saves it using locking and atomic replacement.


#This rule validates and safely saves or reuses the resolved PR2 classifier.
rule ensure_pr2_classifier: # Defines the readiness rule required by downstream PR2 classification.
    input:
        classifier=lambda wc: classifier_source(PR2_DB) # Selects the existing classifier or the staged download or training result.
    output:
        marker=PR2_READY_MARKER # Creates the marker used to confirm that the final PR2 classifier is ready.
    params:
        classifier=PR2_CLASSIFIER, # Provides the final PR2 classifier path.
        lock=PR2_BUILD_LOCK, # Coordinates competing saves to the shared classifier destination.
        read_only=PR2_DB["mode"] == "existing", # Prevents modification of a user-supplied existing classifier.
        md5=PR2_DB["checksum"], # Provides the expected checksum when one is available.
        label="PR2" # Provides the database label used in validation messages.
    conda:
        config["qiime2version"] #Uses the configured QIIME 2 environment for artifact validation.
    log:
        "logs/PR2_classification_db_prep/10-validate-and-save.log" # Records classifier validation, reuse and saving messages.
    resources:
        mem_mb=2000, # Declares 2,000 MB for resource-aware scheduling; does not increase Slurm's memory allocation.
        runtime=30 # Declares an estimated runtime of 30 minutes; does not extend or enforce Slurm's time limit.
    script:
        "../scripts/save_classifier.py" # Validates the classifier and safely saves it using locking and atomic replacement.
