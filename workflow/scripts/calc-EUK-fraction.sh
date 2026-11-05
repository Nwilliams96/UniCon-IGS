#!/usr/bin/env bash

#This script totals the prokaryotic and eukaryotic reads across all samples and writes the two totals to a project-level summary file.

totalPROKseqs=`tail -n+2 ${snakemake_input[eukfracpersample]} | cut -f2 | paste -sd+ | bc` #Skips the header, extracts the prokaryotic-read column and sums its values.

totalEUKseqs=`tail -n+2 ${snakemake_input[eukfracpersample]} | cut -f3 | paste -sd+ | bc` #Skips the header, extracts the eukaryotic-read column and sums its values.

totalSeqs=$(python -c "print($totalEUKseqs + $totalPROKseqs)") #Calculates the combined number of prokaryotic and eukaryotic reads.
eukFrac=`bc <<< "scale=8; $totalEUKseqs/$totalSeqs"` #Calculates the overall eukaryotic read fraction to eight decimal places.

printf "PROK_reads\tEUK_reads\n" >> ${snakemake_output[eukfracall]} #Appends the two output-column headings.
printf "$totalPROKseqs\t$totalEUKseqs\n" >> ${snakemake_output[eukfracall]} #Appends the calculated prokaryotic and eukaryotic read totals.
