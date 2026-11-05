#download_databases.smk prepares the shared database folders and downloads the two reference FASTA files needed to build the BBsplit database, which separates prokaryotic and eukaryotic reads.
#It contains three rules:
#1 initialize_database_directories creates the standard BBsplit, SILVA and PR2 directories beneath the configured database_dir. It then creates a marker file recording that initialization completed.
#2 download_prok_db downloads the prokaryotic reference, SILVA_132_PROK.cdhit95pc.fasta, from OSF.
#3 download_euk_db downloads the eukaryotic reference, SILVA_132_and_PR2_EUK.cdhit95pc.fasta, from OSF.

#OSF_IDs = ["eux4r", "eahds"] #References for PROK/EUK splitting only; taxonomy classifiers are configured separately.
#NAMES = ["and_PR2_EUK", "PROK"]
#database_dir=config["database_dir"]


#This rule creates the standard database directories and then creates a marker file recording successful initialization. Existing directories are left intact. It does not download databases or check their contents.
rule initialize_database_directories: # Defines a Snakemake rule named initialize_database_directories.
    output:
        marker=DATABASE_PREFIX + ".database_directories.ready" # Constructs the marker-file path and names that output marker.
    params:
        bbsplit=DATABASE_PREFIX + "bbsplit-db", # Defines the directory for BBsplit references and its index.
        silva=DATABASE_PREFIX + "classification/SILVA", # Defines the standard SILVA directory beneath the shared database root.
        pr2=DATABASE_PREFIX + "classification/PR2" # Defines the equivalent PR2 directory.
    priority: 51 # Gives this rule scheduling priority 51.
    shell:
        "mkdir -p {params.bbsplit:q} {params.silva:q} {params.pr2:q} " # Creates all three directories. mkdir -p also creates missing parent directories and accepts directories that already exist. Each {params.…} inserts the corresponding path. :q applies shell quoting so paths containing spaces or special characters are handled safely.
        "&& touch {output.marker:q}" # Python joins this string to the preceding one, forming a single command. && runs the next command only if directory creation succeeds. touch creates the marker if absent, or updates its timestamp if present. {output.marker:q} inserts its shell-quoted path.

#This rule downloads the prokaryotic FASTA reference used to build the BBsplit index. It waits for database-directory initialization, saves the reference beneath the configured database root, and uses a helper to verify its SHA-256 checksum.
rule download_prok_db: #Defines a Snakemake rule named download_prok_db.
    input:
        rules.initialize_database_directories.output.marker # Requires the marker produced by initialize_database_directories.
    output:
        DATABASE_PREFIX + "bbsplit-db/SILVA_132_PROK.cdhit95pc.fasta", # Constructs the downloaded FASTA’s destination beneath the configured database root. Kept as a persistent shared reference so Snakemake does not automatically delete it after indexing.
    log:
        "logs/bbsplit_prok_db_download.log", #Specifies where the helper records download messages and errors
    params:
        url="https://osf.io/eahds/download", #Provides the OSF download address. The helper follows redirects to retrieve the file.
        checksum="46683398d8859aeef635fff49112c4d6453d4ebd70bf8b9fee006cb5b84273b1", # Stores the expected file fingerprint.
        checksum_algorithm="sha256", # Tells the helper to calculate that fingerprint using SHA-256. SHA-256 is an algorithm that calculates a digital fingerprint of a file’s contents, called a hash or checksum. Even a tiny change to the file changes its fingerprint. This detects corruption, incomplete downloads or unexpected replacement files.
    resources:
        mem_mb=1024, #Declares a memory requirement of 1,024 MB.
        runtime=30, #Declares an expected runtime of 30 minutes under Snakemake’s runtime-resource convention.
    priority: 50 #Gives this rule scheduling priority 50.
    script:
        "../scripts/download_verified_reference.py" # Tells Snakemake to execute download_verified_reference.py which downloads reference files safely and checks that the downloaded content is what the pipeline expects.

#This rule downloads the eukaryotic reference FASTA used to build the BBsplit index. It waits for directory initialization, verifies the downloaded file using SHA-256, and saves it beneath the configured database root.
rule download_euk_db: #Defines a Snakemake rule named download_euk_db.
    input:
        rules.initialize_database_directories.output.marker # Requires the initialization rule’s marker, ensuring the database directories are prepared before downloading.
    output:
        DATABASE_PREFIX + "bbsplit-db/SILVA_132_and_PR2_EUK.cdhit95pc.fasta", # Constructs the destination path for the combined SILVA/PR2 eukaryotic reference. Kept as a persistent shared reference so Snakemake does not automatically delete it after indexing.
    log:
        "logs/bbsplit_euk_db_download.log",
    params:
        url="https://osf.io/eux4r/download", #Provides the OSF download URL. The helper follows redirects to retrieve the reference.
        checksum="c4a396035271db516a0054f1d750d3081ca5a11831b44756e91590a79ed0f109", # Stores the expected fingerprint of the reference file.
        checksum_algorithm="sha256", # Selects SHA-256 as the algorithm used to calculate that fingerprint.
    resources:
        mem_mb=1024, #Declares a memory requirement of 1,024 MB.
        runtime=30, #Declares a runtime resource of 30 minutes.
    priority: 50 #Gives the rule scheduling priority 50.
    script:
        "../scripts/download_verified_reference.py" # Tells Snakemake to execute download_verified_reference.py
