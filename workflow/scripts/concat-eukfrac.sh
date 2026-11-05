#!/usr/bin/env bash

#This script appends per-sample prokaryotic/eukaryotic read counts to the combined summary and adds a header when the summary is initially empty.

if [ -s ${snakemake_output[eukfracpersample]} ] ; then #Checks whether the combined summary already exists and contains data.

    #When the summary is nonempty, append the new per-sample records without adding another header.
    cat ${snakemake_input[eukfrac]} >> ${snakemake_output[eukfracpersample]} #Appends all supplied per-sample fraction files to the existing summary.

else

    #When the summary is empty or missing, add the header before appending the records.
    printf "sample\tprok_seqs_split\teuk_seqs_split\ttotal_seqs_split\teukfrac_split\n" >> ${snakemake_output[eukfracpersample]} #Writes the combined table’s column headings.
    cat ${snakemake_input[eukfrac]} >> ${snakemake_output[eukfracpersample]} #Appends all supplied per-sample fraction records beneath the header.

fi #Ends the empty-file check.
