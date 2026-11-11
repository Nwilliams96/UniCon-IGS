#Rules that help to make the final html report.

#This chunk selects optional internal-standard report inputs and generates a self-contained HTML summary of workflow configuration, quality control, read retention, taxonomic assignments and corrections.
INTERNAL_STANDARD_REPORT_FIGURES = [] #Defaults to no internal-standard figures when correction is disabled.
INTERNAL_STANDARD_REPORT_TABLE = [] #Defaults to no internal-standard abundance table when correction is disabled.
UNMERGED_16S_REPORT_SUMMARY = [UNMERGED_16S_SUMMARY] #Includes the linked/unmerged 16S recovery summary in the report.
if USE_INTERNAL_STANDARDS: #Adds internal-standard results only when that optional workflow branch is enabled.
    INTERNAL_STANDARD_REPORT_FIGURES = [ #Collects the internal-standard diagnostic figures.
        rules.intstd_correct_data.output.recovery_plot_png, #Includes the standard-recovery ratio plot.
        rules.intstd_correct_data.output.domain_plot_png, #Includes the corrected domain-composition plot.
    ]
    INTERNAL_STANDARD_REPORT_TABLE = [rules.intstd_correct_data.output.corrected] #Includes the corrected long-format abundance table.


#This rule combines configuration, quality-control, denoising, taxonomic and optional correction results into the project HTML report.
rule generate_pipeline_report: #Defines the pipeline-report generation rule.
    input:
        # Track the bundled logo so branding changes regenerate the report.
        logo=workflow.source_path("../assets/unicon-logo.png"),
        # Use the same configured sample sheet as the processing rules.
        samples=config["samplesheet"],
        split_summary="results/" + config["projectName"] + ".eukfrac-per-sample.tsv", #Provides per-sample prokaryotic and eukaryotic read counts.
        stats16s="results/02-proks/04-DADA2d-plaintext-exports/" + config["projectName"] + ".16S.latest_stats.tsv", #Provides the 16S DADA2 retention statistics.
        stats18s="results/02-euks/09-DADA2d-plaintext-exports/" + config["projectName"] + ".18S.latest_stats.tsv", #Provides the 18S DADA2 retention statistics.
        correction_factors=rules.merge_prok_euk.output.correction_factors, #Provides the applied 16S/18S abundance-correction factors.
        species_assignment=rules.assign_SILVA_species.output.summary, #Provides the SILVA species-assignment results.
        chloroplast_audit=rules.resolve_PR2_plastids.output.audit, #Provides the decision recorded for every SILVA/PR2 plastid comparison.
        chloroplast_summary=rules.resolve_PR2_plastids.output.summary, #Provides overall counts for the plastid-resolution outcomes.
        cutadapt_qc=expand( #Collects primer-trimming reports for every sample.
            "results/00-trimmed/{sample}.qc.txt", sample=samples["sample"] #Builds one Cutadapt report path per sample.
        ),
        quality_16s="results/02-proks/02-quality-plots-R1-R2/", #Provides the paired-end 16S quality visualization.
        quality_18s_paired="results/02-euks/02-quality-plots-R1-R2/", #Provides the paired-end 18S quality visualization.
        quality_18s_concatenated="results/02-euks/07-quality-plots-concat/", #Provides the quality visualization for concatenated 18S reads.
        long_data=RESULTS_LONG_DATA, #Provides the selected final long-format abundance table.
        internal_standard_figures=INTERNAL_STANDARD_REPORT_FIGURES, #Provides internal-standard figures when correction is enabled.
        internal_standard_table=INTERNAL_STANDARD_REPORT_TABLE, #Provides corrected internal-standard data when available.
        unmerged_16s_summary=UNMERGED_16S_REPORT_SUMMARY #Provides the linked/unmerged 16S recovery summary.
    output:
        html="results/07-report/" + config["projectName"] + ".pipeline-report.html" #Writes the self-contained project report.
    log:
        "logs/06-report/generate_pipeline_report.log" #Records report-generation messages and errors.
    script:
        "../scripts/generate_pipeline_report.py" #Reads the collected results and renders the final HTML report.
