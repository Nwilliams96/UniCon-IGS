#builds the BBsplit reference index used to separate prokaryotic and eukaryotic reads. It uses the two downloaded FASTA references and stores the index in the shared database directory.
rule prepare_bbsplit_db: #Defines Snakemake rule.
    input:
        file1=DATABASE_PREFIX + "bbsplit-db/SILVA_132_and_PR2_EUK.cdhit95pc.fasta", #Defines file1 as the combined SILVA/PR2 eukaryotic reference fasta
        file2=DATABASE_PREFIX + "bbsplit-db/SILVA_132_PROK.cdhit95pc.fasta", #Defines file2 as the prokaryotic reference FASTA.
    output:
        directory(DATABASE_PREFIX + "bbsplit-db/EUK-PROK-bbsplit-db/") # Defines where the BBSplit index files are kept.
    conda:
        "../envs/bbmap.yaml" #Chooses the BBMap/BBTools environment definition relative to this rule.
    log:
        "logs/bbsplit_db_prep.log" #Identifies log location
    threads: 8 # Requests up to eight threads. Snakemake may reduce this number to fit the available cores.
    resources:
        mem_mb=64000, # Declares a memory requirement of 64,000 MB—approximately 64 GB.
        runtime=240, # Declares a runtime resource of 240 minutes, or four hours. Enforcement depends on the executor or cluster configuration.
    priority: 50 #Gives this rule scheduling priority 50.
    shell:
        "bbsplit.sh build=1 threads={threads} ref={input.file1:q},{input.file2:q} " # build=1: reference-build identifier 1. This is a BBTools identifier, not the workflow’s step number. threads={threads}: the effective thread allocation supplied by Snakemake. ref=…,…: both reference FASTAs, separated by a comma. :q: shell-quotes each input path so spaces and special characters are handled safely.
        "path={output:q} 2> {log:q}" # Python joins this string to the preceding one, completing the command: path={output:q} directs BBsplit to store its index in the declared directory. 2> {log:q} sends standard-error messages to the log, replacing its previous contents. Standard output is not redirected by this expression.
