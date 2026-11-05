#This rule collects the workflow’s final data products, summary tables and documentation in a single directory for download and sharing.
rule export_results_for_download: #Defines the final results-packaging rule.
    input:
        data_files=RESULTS_EXPORT_INPUTS #Requires all files selected for inclusion in the downloadable results package.
    output:
        data_files=RESULTS_EXPORT_OUTPUTS, #Creates packaged copies of the selected result files.
        summaries=RESULTS_EXPORT_SUMMARIES, #Creates the configured taxonomic and workflow summary tables.
        readme=RESULTS_EXPORT_DIR + "/README.txt" #Creates a README describing the exported files.
    script:
        "../scripts/export_results.py" #Copies and summarizes the selected outputs in the final download directory.
