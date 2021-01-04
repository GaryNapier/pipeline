#! /bin/bash

# vcf2fasta.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Get distance matrix generated by plink

# Adapted from:


# Arguments to script:
# study accession, vcf location, output location

# Input:
# filtered vcf file of name format: <study accession>.filt.val.gt.g.vcf.gz, reference fasta file

# Main steps:
# run vcf2fasta.py on VCF and output multi-fasta file

# Output:
# multi-fasta file fro study accession

# RUN
# Assume run from ~/transmission :
# ./shell_scripts/vcf2fasta.sh <study_accession> <vcf_dir> <fasta_output_dir>  <ref_fasta>
# ./shell_scripts/vcf2fasta.sh PRJEB7669 ~/vcf/ fasta/ ~/refgenome/MTB-h37rv_asm19595v2-eg18.fa

# ------------------------------------------------------------------------------

# Setup

# Variables
study_accession=${1?Error: enter sudy accession number e.g. PRJEB7669}

# Directories
vcf_dir=${2?Error: enter directory of input vcf file}
fasta_output_dir=${3?Error: enter when multi-fasta results are to be saved}

# Files
vcf_input_file=${vcf_dir}${study_accession}.filt.val.gt.g.vcf.gz
# output_fasta_file=${fasta_output_dir}${study_accession}.fasta
ref_fasta=${4?Error: enter reference fasta file}

# Parameters

# Commands


# Print arguments
printf "\n"
echo "Arguments:"
printf "\n"
echo ${study_accession}
echo ${vcf_dir}
echo ${fasta_output_dir}
echo ${vcf_input_file}
printf "\n"

# Make output dir if not exist
if [ ! -d ${fasta_output_dir} ]; then
    mkdir ${fasta_output_dir}
fi

# ------------------------------------------------------------------------------

echo "Converting genotyped and filtered multi-VCF to multi-fasta file"
# Make fasta file of genotyped and filtered multi-vcf
vcf2fasta.py --vcf ${vcf_input_file} --ref ${ref_fasta} --threads 20 --snps

# Move results to output directory
mv ${vcf_dir}${study_accession}.filt.val.gt.g.snps.fa ${fasta_output_dir}
