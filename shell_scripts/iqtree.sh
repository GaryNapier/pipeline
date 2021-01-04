#! /bin/bash

# iqtree.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Generate tree of study accession samples

# Arguments to script:
# study accession, fasta location, output location

# Input:
# fasta file of samples from one study accession

# Main steps:
# run IQ tree on fasta files

# Output:
# newick file of study accession tree

# RUN
# Assume run from ~/transmission :
# ./shell_scripts/iqtree.sh <study_accession>, <fasta_dir>, <newick_output_dir>
# ./shell_scripts/iqtree.sh PRJEB7669 fasta/ newick/

# ------------------------------------------------------------------------------

# Setup

# Variables
study_accession=${1?Error: enter sudy accession number e.g. PRJEB7669}

# Directories
fasta_dir=${2?Error: enter directory of input fasta file}
newick_output_dir=${3?Error: enter where iqtree results are to be saved}

# Files
fasta_input_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa

# Parameters

# Commands

# Print arguments
printf "\n"
echo "Arguments:"
printf "\n"
echo ${study_accession}
echo ${fasta_dir}
echo ${newick_output_dir}
echo ${fasta_input_file}
# echo ${output_newick_file}
printf "\n"

# Make output dir if not exist
if [ ! -d ${newick_output_dir} ]; then
    mkdir ${newick_output_dir}
fi

# ------------------------------------------------------------------------------

# Run iqtree:

echo "Running iqtree on ${fasta_input_file}"
# Run tree
iqtree -s ${fasta_input_file} -m GTR+G -nt AUTO

# Move output files to correct place
mv ${fasta_dir}*.iqtree ${fasta_dir}*.treefile ${fasta_dir}*.log ${newick_output_dir}
