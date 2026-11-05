#!/usr/bin/env bash

#This script creates nucleotide BLAST databases for three internal-standard reference FASTAs and then creates their Snakemake completion markers.

makeblastdb -dbtype nucl -in ${snakemake_input[intstd1]} #Builds a nucleotide BLAST database from the first internal-standard FASTA.
touch ${snakemake_output[intstd1]} #Creates the first completion marker expected by Snakemake.

makeblastdb -dbtype nucl -in ${snakemake_input[intstd2]} #Builds a nucleotide BLAST database from the second internal-standard FASTA.
touch ${snakemake_output[intstd2]} #Creates the second completion marker expected by Snakemake.

makeblastdb -dbtype nucl -in ${snakemake_input[intstd3]} #Builds a nucleotide BLAST database from the third internal-standard FASTA.
touch ${snakemake_output[intstd3]} #Creates the third completion marker expected by Snakemake.
