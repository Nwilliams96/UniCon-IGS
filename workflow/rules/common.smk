#common.smk is the pipeline’s shared setup and configuration file. It translates the user’s configuration into settings, paths and required outputs that the other rules use consistently.
#Its main responsibilities are:
# Resolve settings: read configuration choices and apply defaults for QIIME, database versions, DADA2 and chloroplast classification.
# Set up database references: determine classifier locations, reuse/download/training routes, readiness markers and supporting reference paths.
# Configure internal standards: read their IDs, check naming, generate combinations and define expected output filenames. The actual correction calculations happen elsewhere.
# Read sample metadata: load samples.tsv, preserve sample IDs as text, sort samples and check required ISD columns.
# Define completion targets: assemble the files and directories required for processing, reports, supplementary outputs and Results-Export.
# Validate inputs: check the loaded configuration and sample sheet against their schemas.
#It does not itself run trimming, DADA2, classification or correction. It prepares the shared information those processing rules need.


#import basic packages
import os # handles filesystem paths and operating system functions
import re # Searches and validates text using regular-expression patterns.
from itertools import combinations # Imports combinations from itertools to enumerate internal-standard pairs for the expected workflow outputs; it does not install anything.
import pandas as pd
from snakemake.exceptions import WorkflowError # Imports workflow errors allowing for easy error explanation when something doesn't work.
from snakemake.utils import validate

#Reads the BBsplit reuse setting, defaulting to False when absent. False permits database preparation when needed; completed files may be reused. SILVA and PR2 independently select existing, download or training modes through database_config.
USE_PREEXISTING_BBSPLIT_DATABASE = config.get(
    "use_preexisting_bbsplit_database", False
)
#Enables the internal-standard module when True; when False, its processing rules are skipped.
USE_INTERNAL_STANDARDS = config.get("use_internal_standards", False)

#Stores the default YAML path in a variable; this assignment does not select or install an environment.
DEFAULT_QIIME2_ENV = "../envs/rachis-qiime2-linux-64-2026.7.yml"
configured_qiime2_env = str(config.get("qiime2version", "")).strip() # str() converts the value to text. .strip() removes surrounding spaces and tabs from that text.
if not configured_qiime2_env:
    config["qiime2version"] = DEFAULT_QIIME2_ENV # if not checks whether the resulting text is empty. This is true when the setting is missing, empty, or has only whitespace. If it is empty, it assigns the default YAML path to the running configuration. Otherwise, the users explicit setting is left unchanged. This won't re-write the saved config.yml or install the environment. Snakemake handles environment creation or re-use later when run with --use-conda.

#Reads the version of each database from the config.  Writes the resolved database release back to the running configuration, so that rules and report downstream use the same values.
SILVA_VERSION = str(config.get("SILVAversion", "144"))
PR2_VERSION = str(config.get("PR2version", "5.1.1"))
config["SILVAversion"] = SILVA_VERSION
config["PR2version"] = PR2_VERSION
for database in ("silva", "pr2"): # runs the following statement x2, first with database = "silva", then with database = "pr2"
    config[database + "_reference_mode"] = config.get(database + "_reference_mode", "amplicon") # builds setting names + preserves config value, default is amplicon if missing. This is for database training. "amplicon" will extract the region between configured fwd primer and rev primer before training. "full_length" will skip the primer region extraction and train on prepared full length reference sequences. Reference cleaning and dereplication still occur.

#Sets the minimum PR2 confidence for plastid reconciliation. With the default threshold, lineages containing :plas and confidence >= 0.7 can replace SILVA assignments.
CHLOROPLAST_PR2_MIN_CONFIDENCE = float( # converts the supplied value into a decimal number. Reads the threshold from the configuration. If the setting is missing it will use 0.7.
    config.get("chloroplast_pr2_min_confidence", 0.7))
if not 0 <= CHLOROPLAST_PR2_MIN_CONFIDENCE <= 1: # checks whether the threshold falls outside of 0-1, including neither values below zero or above 1 as valid. 
    raise WorkflowError("chloroplast_pr2_min_confidence must be between 0 and 1") # Will stop workflow with a clear error with a value outside of the threshold.

#Reads the DADA2 configuration dictionary and its prokaryote/16S and eukaryote/18S subsections. Missing sections default to empty dictionaries.
DADA2_CONFIG = config.get("dada2", {})
DADA2_PROK_CONFIG = DADA2_CONFIG.get("prokaryotes", {})
DADA2_EUK_CONFIG = DADA2_CONFIG.get("eukaryotes", {})

#These are backups, in case a user makes their own config and misses one of these settings, but snakemake expects it.
DADA2_PROK_DEFAULTS = {
    "trunc_len_f": 220,
    "trunc_len_r": 180,
    "max_ee_f": 2.0,
    "max_ee_r": 2.0,
    "trunc_q": 2,
    "min_overlap": 12,
    "max_merge_mismatch": 0,
    "trim_overhang": False,
    "pooling_method": "independent",
    "chimera_method": "consensus",
    "min_fold_parent_over_abundance": 1.0,
    "n_reads_learn": 1000000,
}
#These are backups, in case a user makes their own config and misses one of these settings, but snakemake expects it.
DADA2_EUK_DEFAULTS = {
    "max_ee": 2.0,
    "trunc_q": 0,
    "pooling_method": "independent",
    "chimera_method": "consensus",
    "min_fold_parent_over_abundance": 1.0,
    "n_reads_learn": 1000000,
}

DADA2_PROK = {**DADA2_PROK_DEFAULTS, **DADA2_PROK_CONFIG} # Creates effective 16S settings by combining two dictionaries. ** unpacks their contents. It starts with defaults and then applies users config afterwards (sort of like an over-write where available).
DADA2_EUK = {**DADA2_EUK_DEFAULTS, **DADA2_EUK_CONFIG} # Creates effective 18S settings by combining two dictionaries. ** unpacks their contents. It starts with defaults and then applies users config afterwards (sort of like an over-write where available).
# Internal standards use the ordered YAML list generated by the configuration page.

INTERNAL_STANDARD_IDS = config.get("intstds", []) # Reads the internal standards list from the config. If "intstds" is missing, it uses an empty list []. Later validation rejects an empty list if the internal standard correction is enabled.

#In summary, this chunk, when ISD correction is enabled, checks that at least 1 standard ID is supplied, every ID contains only permitted characters, & no IDs are duplicated. If any checks fail, workflow stops here.
if USE_INTERNAL_STANDARDS: # Run checks when ISD correction is enabled
    if not INTERNAL_STANDARD_IDS: # Checks if the ISD list may be empty
        raise WorkflowError( # Stops the workflow if correction is enabled, but no standards were specified.
            "use_internal_standards is true, but intstds does not contain any IDs"
        )
    invalid_internal_standard_ids = [ # Starts building a list of IDs and include it in the new list only if the next condition is true.
        standard_id
        for standard_id in INTERNAL_STANDARD_IDS
        if not re.fullmatch(r"[A-Za-z0-9._-]+", standard_id) # Checks whether the entire ID fails to match the permitted characters: ASCII letters, digits, periods, underscores and hyphens.
    ]
    if invalid_internal_standard_ids: # Checks if list contains entries. 
        raise WorkflowError( # Stops the workflow when one or more supplied internal-standard IDs contain invalid characters.
            "Internal-standard IDs may contain only letters, numbers, periods, "
            "underscores, and hyphens: "
            + ", ".join(invalid_internal_standard_ids) # separates the invalid names with commas.
        )
    if len(set(INTERNAL_STANDARD_IDS)) != len(INTERNAL_STANDARD_IDS):
        raise WorkflowError("Internal-standard IDs in intstds must be unique") # stops the workflow when duplicate IDs exist. This prevents two standards from sharing the same ID and consequently the same associated columns and output paths.

#This block prepares the internal standard labels, combinations and output file names that snakemake needs. It supports any number of standards and puts their actual IDs into correction method names. It doesn't calculate correction factors - R script does that.
INTERNAL_STANDARD_SLOTS = [
    f"isd_{index}" for index in range(1, len(INTERNAL_STANDARD_IDS) + 1) # creates one label per configured standard. - len counts the standards, -range generates numbers from 1 through to that count. - f"isd inserts each number into the label.
]
INTERNAL_STANDARD_ID_PAIRS = list(combinations(INTERNAL_STANDARD_IDS, 2)) # Creates every distinct pair of standards for three IDs [("BP", "DR")] etc. These are so that snakemake will expect all of these .tsv outputs at the end.
ALL_INTERNAL_STANDARD_IDS_STEM = "_".join(INTERNAL_STANDARD_IDS) # Joins IDs with underscores, preserving their configured order.
ISD_CORRECTED_LONG_TABLE = ( # starts constructing path to internal_standard_corrected_long
    "results/05-internal-std-corrected/" # sets output directory.
    + config["projectName"] # Adds project name.
    + "." # separates with a period.
    + ALL_INTERNAL_STANDARD_IDS_STEM # IDs adds configured standard IDs with underscores.
    + ".ISD_corrected_long_data.tsv" # adds descriptive file name ending and tsv extension.
)
ISD_FILTERED_LONG_TABLE = ( # starts constructing a separate output path for the long table.
    "results/05-internal-std-corrected/" # Sets the directory for the ISD-removed long table, separate from the ISD-corrected long table.
    + config["projectName"] # Adds project name.
    + ".ISD_removed_long_data.tsv" # adds descriptive file name ending and tsv extension.
)
INTERNAL_STANDARD_METHOD_STEMS = ( # starts building a combined list of corrected method names.
    [f"{standard_id}_recovery_ratio" for standard_id in INTERNAL_STANDARD_IDS] # creates one name per individual standard.
    + [
        f"mean_{ALL_INTERNAL_STANDARD_IDS_STEM}_recovery_ratio", # names method using mean recovery ratio across all configured standards.
        f"median_{ALL_INTERNAL_STANDARD_IDS_STEM}_recovery_ratio", # names the corresponding median method.
    ]
    + [
        f"mean_{first}_and_{second}_recovery_ratio" # constructs a mean name from two standard IDs.
        for first, second in INTERNAL_STANDARD_ID_PAIRS # repeats the naming operation for each pair mentioned earlier.
    ]
)

#Later the merging rules use this variable to identify the table that has the 16S + 18S amplicon molarities for correction.
AMPLICON_CONCENTRATIONS_FILE = ( # Starts to assign a file path
    "config/prok_and_euk_SSU_amplicon_molarities.tsv" # path
)

# Keep shared supporting resources separate from explicit, read-only classifiers.
DATABASE_DIR = os.path.normpath(os.path.expanduser(config["database_dir"])) # os.path.normpath simplifies separators and trailing slashes. os.path.expanduser expands a leading ~ or supported ~user form to a home directory. config[] retrieves the required setting.
DATABASE_PREFIX = DATABASE_DIR + os.sep # Adds operating systems directory separator -/ on Linux and macOS - so rules that append text to this prefix have a separating /
BBSPLIT_DB_DIR = os.path.join(DATABASE_DIR, "bbsplit-db", "EUK-PROK-bbsplit-db") # Constructs BBSplit index directory path safely.

#This block loads the local configuration helper, resolves the selected QIIME environment, and determines how SILVA and PR2 will be obtained and where their files belong. It also converts configuration errors into readable snakemake errors.
import importlib.util # loads python utilities for importing a module directly from a file path.
database_spec = importlib.util.spec_from_file_location( # Starts creating an import specification: instructions telling python which module to load and where it is. 
    "easv_database_config", str(workflow.source_path("../scripts/database_config.py"))
)
database_config = importlib.util.module_from_spec(database_spec) # creates module object
database_spec.loader.exec_module(database_config) # Executes the helper file's top-level code, making its functions available through database_config.

#This block is working out which QIIME environment and database resources the workflow should use. It inspects files and envs but does not install, download or train.
try: # Python attempts the indented statements below. If a ValueError is triggered it will jump to the WorkflowError at the bottom.
    qiime_specification = config["qiime2version"] # Read configured environment. Can be YAML path or the name of an already installed env.
    # Resolve relative YAML paths here, in the rule file's source context.
    qiime_definition = (
        workflow.source_path(qiime_specification) # resolves the path relative to the workflow source location.
        if qiime_specification.endswith((".yaml", ".yml")) else None # Checks whether the value looks like a yaml.
    )
    qiime_environment = database_config.resolve_qiime_environment( # Calls project helper :resolve_qiime_environment() - for a yaml it reads the file and calculates the fingerprint of its contents, and for an installed environment name, it finds the corresponding Conda environment and inspects the installed package inventory.
        qiime_specification, lambda specification: qiime_definition # When the helper asks for the YAML’s location, it returns the path already resolved above.
    )
    config["qiime2version"] = qiime_environment["specification"] # Updates the running configuration with the resolved specification. - For a named environment, it becomes the exact installed environment directory that was inspected. This does not rewrite the saved config.yml.
    config["qiime2_environment"] = qiime_environment # Stores the complete information returned by the helper, including its type and fingerprint. Later code can access this information without resolving it again.
    database_settings = dict(config, _database_environment_id=qiime_environment["identity"]) # Makes a shallow copy of config, with an additional internal setting called _database_environment_id. That value is the environment fingerprint: YAML input: derived from the YAML’s contents, Named environment: derived from its installed package inventory. The database helper uses this identifier in locally trained classifier filenames to distinguish builds associated with different environment specifications.
    DATABASE_RESOURCES = { # This starts a "dictionary" titled DATABASE_RESOURCES - each entry has a database name and its key resolved to that value
    database: database_config.resolve_database(database_settings, database) #creates dictionary entry.
    for database in ("SILVA", "PR2") #tells helper script what is being resolved.
}
except ValueError as error: # Catch configuration-validation errors
    raise WorkflowError(str(error)) from error # Report the problem to Snakemake

#Takes database information already resolved by the helper script and gives it a shorter name for use in the snakemake rules. It also enables the optional SILVA 144 species-matching step and makes the database info available for the final report.
SILVA_DB = DATABASE_RESOURCES["SILVA"] # Turns DATABASE_RESOURCES["SILVA"] into "SILVA_DB". Likewise for the rest.
PR2_DB = DATABASE_RESOURCES["PR2"]
SILVA_CLASSIFIER = SILVA_DB["classifier"] # Stores the path to the final SILVA classifier .qza file. Classification rules use this path.
PR2_CLASSIFIER = PR2_DB["classifier"] # Stores the equivalent final classifier path for PR2.
SILVA_READY_MARKER = SILVA_DB["ready"] # Stores the path of a small project-local marker file.
PR2_READY_MARKER = PR2_DB["ready"] # These markers let Snakemake track readiness without declaring the shared classifier itself as a project-owned output subject to cleanup.
PR2_BUILD_LOCK = PR2_DB["lock"]
SILVA_BUILD_ROOT = SILVA_DB["build_root"] # Get SILVAs working directory.
PR2_BUILD_ROOT = PR2_DB["build_root"] # Get PR2’s working directory.
RUN_SILVA_SPECIES_MATCHING = SILVA_VERSION == "144" # Decide whether to enable SILVA 144 species matching.
config["database_resources"] = DATABASE_RESOURCES # Make the database information available to the report.

#Creates a list of capability-check file paths for databases that still need to be trained. Databases being reused or downloaded—and trained classifiers already considered reusable—are excluded. These markers also precede primer removal, so unsupported releases fail before FASTQ work.
DATABASE_SETUP_CHECKS = [ # Starts a list called DATABASE_SETUP_CHECKS.
    os.path.join(record["build_root"], "capability.json") # Constructs the path to capability.json inside that database’s project build directory.
    for record in DATABASE_RESOURCES.values() # Loops through the resolved database settings—one record for SILVA and another for PR2.
    if record["mode"] == "train" # Only includes databases configured for local training. Existing-classifier and download modes are excluded.
    and not database_config.reusable_classifier(record["classifier"]) # Also requires that the helper does not consider the destination classifier reusable. If a reusable classifier already exists, its capability-check path is unnecessary and excluded.
]

# Defines the filename, download address, checksum and local location of the SILVA 144 reference used for exact species matching. It is different to the trained taxonomy classifier. These lines define and don't download anything.
SILVA_SPECIES_REFERENCE_FILENAME = "silva_v144_assignSpecies.fa.gz" # Stores the reference filename. .fa.gz is a gzip-compressed FASTA file.
SILVA_SPECIES_REFERENCE_URL = ( # Begins defining the download URL.
    "https://www.arb-silva.de/fileadmin/silva_databases/release_144/" # Provides the SILVA server address and release-specific directory.
    "DADA2/1.36.0/SSU/" + SILVA_SPECIES_REFERENCE_FILENAME # Python automatically joins the adjacent URL strings. Then + appends the filename, producing the complete download address. DADA2/1.36.0/SSU/ identifies the reference’s location on the server; it does not select or install a DADA2 version.
)
SILVA_SPECIES_REFERENCE_MD5 = "444de7c0cce0b66addda7a3f8b38e012" # Stores the expected MD5 checksum. The downloader can compare the downloaded file against this value to detect corruption or unexpected content.
SILVA_SPECIES_REFERENCE = os.path.join( # Begins constructing the local destination path. 
    DATABASE_DIR, "classification", "SILVA", SILVA_SPECIES_REFERENCE_FILENAME # Joins the configured database root, subdirectories and filename safely.
)

#This will stop the workflow if a local BBSplit database has been selected, but does not exist.
if USE_PREEXISTING_BBSPLIT_DATABASE and not os.path.isdir(BBSPLIT_DB_DIR): # checks the config says to reuse an existing database, and that the configured path is not an existing directory. Both must be true to trigger the error.
    raise WorkflowError("Pre-existing BBsplit database was not found: " + BBSPLIT_DB_DIR) # Stop right there, thank you very much - Stops Snakemake with an explanatory message. + BBSPLIT_DB_DIR appends the exact directory it could not find, helping the user correct the path.

#This block validates the project name, defines the project’s Results-Export folder, and lists the files to include. It selects the ISD-removed ecological table when internal standards are enabled, always includes the separate unmerged-16S supplement, and adds the ISD-corrected table when applicable. It defines paths and lists only. The export rule later uses them to create the folder, copy files and generate summaries.
#Project name and main table
PROJECT_NAME = config["projectName"] # Reads the required project name from the configuration and gives it the shorter variable name PROJECT_NAME.
if not re.fullmatch(r"[A-Za-z0-9._-]+", PROJECT_NAME): # Checks that the entire name contains only letters, numbers, periods, underscores or hyphens. + requires at least one character. 
    raise WorkflowError(
        "projectName may contain only letters, numbers, periods, underscores, "
        "and hyphens"
    ) #Starts an error that stops Snakemake if the name is invalid.
RESULTS_EXPORT_DIR = PROJECT_NAME + "-Results-Export" # Constructs the export-folder name. For example, AMT29 becomes AMT29-Results-Export.
RESULTS_LONG_DATA = ( # Begins choosing which main long-format table to export.
    ISD_FILTERED_LONG_TABLE # Selects the table with detected internal-standard ASVs removed when internal-standard processing is enabled (see if statement below).
    if USE_INTERNAL_STANDARDS
    else "results/04-formatted/" + config["projectName"] + ".long_data.tsv" # If no ISD, selects the ordinary formatted long table, such as results/04-formatted/AMT29.long_data.tsv
)

#Main export inputs and destinations
RESULTS_EXPORT_INPUTS = [ # Starts the list of source files required for the export.
    RESULTS_LONG_DATA, # Includes the main long table selected above.
    "results/04-formatted/" + config["projectName"] + ".asv_sequences.tsv", # Includes the table containing ASV sequences.
    "results/07-report/" + config["projectName"] + ".pipeline-report.html", # Includes the final HTML report.
    "results/02-proks/09-subsetting/tax-merged/" # Starts the source path for the PR2 plastid-routing audit. 
    + config["projectName"] # Adds the project-name prefix.
    + ".PR2-plastid-routing-audit.tsv", # Adds the audit filename ending. The comma finishes this list item.
    "results/02-proks/09-subsetting/tax-merged/" # Starts the path for the corresponding routing summary.
    + config["projectName"] # Adds the project-name prefix.
    + ".PR2-plastid-routing-summary.tsv", # Adds the summary filename ending.
]
RESULTS_EXPORT_OUTPUTS = [ # Starts the list of destination paths. Its order matches the input list: first source goes to first destination, and yada yada.
    RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".long_data.tsv", # Defines the exported main table’s destination.
    RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".asv_sequences.tsv", # Defines the exported sequence table’s destination.
    RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".pipeline-report.html", # Defines the exported report’s destination.
    RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".PR2-plastid-routing-audit.tsv", # Defines the exported plastid-routing audit’s destination.
    RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".PR2-plastid-routing-summary.tsv", # Defines the exported plastid-routing summary’s destination.
]
#Generated taxonomy summaries
RESULTS_EXPORT_SUMMARIES = [ # Starts a separate list of summary tables that the export helper generates, rather than simply copies. 
    RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".phylum_summary.tsv", # Defines the destination for the phylum-level summary.
    RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".order_summary.tsv", # Defines the destination for the order-level summary.
]
#Separate unmerged-16S supplement. Always export linked/unmerged features separately from canonical merged ASVs.
UNMERGED_16S_TABLE = ( # Begins constructing the supplementary table’s source path.
    "results/02-proks/11-unmerged-16S-rescue/" # Sets its source directory.
    + config["projectName"] # Adds the project-name prefix.
    + ".unmerged_linked_16S.tsv" # Adds the supplementary table’s filename ending.
)
UNMERGED_16S_SUMMARY = ( # Begins constructing the supplementary summary’s source path.
    "results/02-proks/11-unmerged-16S-rescue/" # Uses the same source directory.
    + config["projectName"] # Adds the project-name prefix.
    + ".unmerged_linked_16S_summary.tsv" # Adds the supplementary summary’s filename ending.
)
RESULTS_EXPORT_INPUTS.extend([UNMERGED_16S_TABLE, UNMERGED_16S_SUMMARY]) # Adds both source paths individually to the existing input list. .extend() adds multiple items.
RESULTS_EXPORT_OUTPUTS.extend( # Begins adding their corresponding destination paths.
    [
        RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".unmerged_linked_16S.tsv", # Sets the supplementary linked/unmerged-16S table's export destination.
        RESULTS_EXPORT_DIR + "/" + config["projectName"] + ".unmerged_linked_16S_summary.tsv", # Sets the supplementary linked/unmerged-16S summary's export destination.
    ]
)
#Optional internal-standard-corrected export
if USE_INTERNAL_STANDARDS: # Runs the following additions only when internal-standard processing is enabled.
    RESULTS_EXPORT_INPUTS.append(ISD_CORRECTED_LONG_TABLE) # Adds the ISD-corrected long table to the source list. .append() adds one item.
    RESULTS_EXPORT_OUTPUTS.append( # Begins adding the corrected table’s export destination.
        RESULTS_EXPORT_DIR + "/" # Starts with the export-folder path and a directory separator.
        + config["projectName"] # Adds the project-name prefix.
        + "." # Adds a period separating the project name from the standard IDs.
        + ALL_INTERNAL_STANDARD_IDS_STEM # Adds the configured standard IDs joined with underscores—for example, BP_DR_TT.
        + ".ISD_corrected_long_data.tsv" # Adds the descriptive filename ending. 
    )

# read sample sheet
samples = ( #Begins assigning the resulting pandas table—a DataFrame—to samples. The parentheses let the expression span multiple lines.
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str}) #Reads the file specified by samplesheet in the configuration
    .set_index("sample", drop=False) #Uses the sample values as row labels, allowing downstream code to look up rows by sample ID. drop=False keeps sample as an ordinary column too.
    .sort_index() #Sorts rows in ascending order by sample ID. Because IDs are text, this is text sorting, not numerical sorting—for example, 1, 10, 2.
)

#checks that the sample sheet contains an added-mass column for every configured internal standard when ISD processing is enabled. If any columns are missing, it stops the workflow and lists them. It checks column names, not the values inside them.
if USE_INTERNAL_STANDARDS: # Runs the following checks only when internal-standard processing is enabled.
    expected_internal_standard_columns = [ # Starts building a list of required column names.
        f"{standard_id}_ng" for standard_id in INTERNAL_STANDARD_IDS  # Loops through the configured standard IDs and adds _ng to each one. For example, ["BP", "DR"] becomes ["BP_ng", "DR_ng"].  
        ]
    missing_internal_standard_columns = [ # Starts building a second list containing only the missing columns.
        column # Specifies what to include in that list: the column name itself.
        for column in expected_internal_standard_columns # Checks each required column name in turn.
        if column not in samples.columns # Includes the name only if it is absent from the loaded sample sheet.
    ]
    if missing_internal_standard_columns: # Checks whether the list contains any missing columns. An empty list does not trigger the error.
        raise WorkflowError( # Starts an error that stops Snakemake.
            "The sample sheet is missing internal-standard amount column(s): " # error message
            + ", ".join(missing_internal_standard_columns) # Joins the missing names with commas and appends them to the message.
        )

#This defines a function and begins building its list of required pipeline outputs. It generates paths for compressed R1 FASTQ files for every sample, separately for the prokaryotic and eukaryotic branches. The commented-out lines are inactive alternatives.
def get_final_output():# Defines the following function.
    final_output = expand( #Calls Snakemake's expand() function and assigns the resulting list of paths to final_output.
        "results/01-split/{sample}.{organism}.R1.fastq.gz", sample=samples["sample"], organism=["prok","euk"] #Provides the filename template.
        #"results/01-split/{sample}.prok.fastq",
        #"results/00-trimmed/{sample}.{direction}.fastq",
        #sample=samples["sample"], direction=["1","2"]
    )

#This adds the required database, 16S, 18S, merged-data, report and downloadable-export products to final_output. Internal-standard-corrected data is required only when that module is enabled; the separate unmerged-16S supplement is always required.These statements list targets—they do not execute processing. Snakemake determines the execution order from rule dependencies, not the order of this list. final_output is a Python list of files and directories that Snakemake should produce before considering the pipeline complete. It is essentially the pipeline’s completion checklist.

    final_output.extend([BBSPLIT_DB_DIR, SILVA_READY_MARKER, PR2_READY_MARKER]) # Adds three targets individually: BBSPLIT_DB_DIR: the BBsplit index directory. SILVA_READY_MARKER: the marker recording successful SILVA classifier preparation or validation. PR2_READY_MARKER: the equivalent PR2 marker.
  
   #16S outputs
   #final_output.append("results/02-proks/manifest.tsv"),
    final_output.append("results/02-proks/16S.qza"), #Requires the imported 16S read artifact.
    final_output.append("results/02-proks/02-quality-plots-R1-R2/"), # Requires the directory containing forward/reverse read-quality visualizations.
    final_output.append("results/02-proks/03-DADA2d/"), # Requires the directory containing 16S DADA2 outputs.
    final_output.append("results/02-proks/04-DADA2d-plaintext-exports"), # Requires the directory containing exported DADA2 results, such as sequences, counts and statistics.
    final_output.append("results/02-proks/05-classified"), # Requires the directory containing the initial taxonomy-classification results.
    final_output.append("results/02-proks/07-SILVA-only-barplots/"), # Requires the SILVA-only taxonomy barplots, before PR2 plastid reconciliation.
    final_output.append("results/02-proks/09-subsetting/reclassified/all_16S_ASVs_PR2.classified.qza"), # Requires PR2 assignments for all canonical 16S ASVs, used to identify plastid lineages.
    final_output.append("results/02-proks/09-subsetting/tax-merged/chloroplasts-PR2-reclassified-merged-classification.qza"), # Requires the resolved taxonomy artifact combining the SILVA assignments with qualifying PR2 plastid assignments.
    final_output.append("results/02-proks/09-subsetting/tax-merged/" + config["projectName"] + ".PR2-plastid-routing-audit.tsv"), # Requires the detailed audit of plastid-routing decisions.
    final_output.append("results/02-proks/09-subsetting/tax-merged/" + config["projectName"] + ".PR2-plastid-routing-summary.tsv"), # Requires the summary of those routing decisions.
    final_output.append("results/02-proks/10-exports/" + config["projectName"] + ".taxonomy.tsv"), # Requires the final exported 16S taxonomy table.
    final_output.append("results/02-proks/10-exports/" + config["projectName"] + ".all-16S-seqs.with-tax.tsv"), # Requires the exported 16S abundance table with taxonomy.
    final_output.append("results/02-proks/sample-metadata.tsv"), #Requires the 16S sample-metadata table.
    
    #18S outputs
    final_output.append("results/02-euks/18S-viz.qza"), #Requires the imported paired-read artifact used for the initial 18S quality assessment. Despite viz in its name, .qza denotes a data artifact rather than a .qzv visualization.
    final_output.append("results/02-euks/02-quality-plots-R1-R2/"), # Requires quality visualizations for the original 18S forward/reverse reads.
    final_output.append("results/02-euks/07-quality-plots-concat"), # Requires quality visualizations for the concatenated 18S reads.
    final_output.append("results/02-euks/08-DADA2d"), # Requires the directory containing 18S DADA2 outputs.
    final_output.append("results/02-euks/18S-concat.qza"), #Requires the imported concatenated-read artifact used for single-end 18S DADA2 processing.
    final_output.append(expand("results/02-euks/04-concatenated/{sample}.euk.concatenated.fastq", sample=samples["sample"])), #expand() generates one concatenated FASTQ path for every sample in the sample sheet. .append() adds that resulting list as a nested item; Snakemake flattens nested target lists when consuming them.
    final_output.append("results/02-euks/09-DADA2d-plaintext-exports/"), # Requires the exported 18S DADA2 results.
    final_output.append("results/02-euks/10-classified/" + config["projectName"] + "_PR2.classified.qza"), # Requires the PR2 taxonomy-assignment artifact.
    final_output.append("results/02-euks/10-classified/" + config["projectName"] + "_PR2.normalized.classified.qza"), # Requires the artifact after PR2 rank-prefix normalization. 
    final_output.append("results/02-euks/sample-metadata.tsv"), #Requires the 18S sample-metadata table.
    final_output.append("results/02-euks/12-PR2-barplots/"), # Requires the PR2 taxonomy barplots.
    final_output.append("results/02-euks/15-exports/" + config["projectName"] + ".include_Metazoa_PR2_filtered_table.with-tax.tsv"), # Requires the taxonomy-annotated table containing features assigned to Metazoa.
    final_output.append("results/" + config["projectName"] + ".eukfrac-per-sample.tsv"), # Requires the per-sample 16S/18S read-splitting summary.
    final_output.append("results/" + config["projectName"] + ".eukfrac-all.tsv"), # Requires the combined read-count summary used downstream.
    final_output.append("results/03-merged/" + config["projectName"] + ".merged_uncorrected.tsv"), # Requires the combined 16S/18S table before the downstream count corrections.
    final_output.append("results/02-proks/10-exports/" + config["projectName"] + ".Synechococcales.proportal-classified.tsv"), # Requires the ProPortal classification results for the targeted cyanobacterial features.
    final_output.append("results/04-formatted/" + config["projectName"] + ".long_data.tsv"), # Requires the final formatted long table.
    final_output.append("results/07-report/" + config["projectName"] + ".pipeline-report.html") #Constructs the project-specific HTML report path and adds it as a required target.
    
    #Optional internal-standard output
    if USE_INTERNAL_STANDARDS: # Runs the following statement only when internal-standard processing is enabled.
        final_output.append(ISD_CORRECTED_LONG_TABLE) # Adds the internal-standard-corrected long table as a required target.
    final_output.extend([UNMERGED_16S_TABLE, UNMERGED_16S_SUMMARY]) # Always adds the supplementary linked/unmerged-16S table and its summary. These remain separate from canonical merged-ASV results.
    final_output.extend(RESULTS_EXPORT_OUTPUTS) # Adds every destination file previously listed for the project’s Results-Export folder.
    final_output.extend(RESULTS_EXPORT_SUMMARIES) # Adds the generated phylum- and order-level summary tables.
    final_output.append(RESULTS_EXPORT_DIR + "/README.txt") # Adds the export folder’s explanatory README.

    return final_output #Sends the completed list of required files and directories back to whatever called get_final_output(), then ends the function.
 
#Checks the loaded sample sheet and configuration against their schemas—files defining requirements such as mandatory columns, value types and permitted formats. Snakemake stops with an error if either fails validation.
validate(samples, schema="../../config/schemas/samples.schema.yml") #Checks the loaded samples table against the sample-sheet schema.
validate(config, schema="../../config/schemas/config.schema.yml") #Checks the loaded configuration dictionary against the configuration schema.
