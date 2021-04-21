#! /bin/bash

# vcf2fasta.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Get fasta file from multi-vcf input file

# Adapted from:
# https://github.com/pathogenseq/fastq2matrix/tree/master/scripts
# vcf2fasta.py
# And see https://github.com/pathogenseq/fastq2matrix/blob/master/fastq2matrix/vcf.py for vcf_class

# Arguments to script:
# study accession, vcf location, output location

# Input:
# filtered vcf file of name format: <study accession>.filt.val.gt.g.vcf.gz, reference fasta file

# Main steps:
# run vcf2fasta.py on VCF and output multi-fasta file

# Output:
# multi-fasta file of study accession samples

# RUN
# Assume run from ~/transmission :
# ./shell_scripts/vcf2fasta.sh <study_accession> <vcf_file>                                       <fasta_output_dir>   <ref_fasta>
# ./shell_scripts/vcf2fasta.sh PRJEB7669         ~/vcf/PRJEB7669.filt.val.gt.g.vcf.gz             fasta/               ~/refgenome/MTB-h37rv_asm19595v2-eg18.fa

# ------------------------------------------------------------------------------

# Setup

# Variables
study_accession=${1?Error: enter sudy accession number e.g. PRJEB7669}

# Files
# vcf_input_file=${vcf_dir}${study_accession}.filt.val.gt.g.vcf.gz
vcf_input_file=${2?Error: enter name of input vcf file}
# output_fasta_file=${fasta_output_dir}${study_accession}.fasta
ref_fasta=${4?Error: enter reference fasta file}

# Directories
vcf_dir=`dirname ${vcf_input_file}`
fasta_output_dir=${3?Error: enter where multi-fasta results are to be saved}

# Define output fast file name - see vcf_to_fasta() function in https://github.com/pathogenseq/fastq2matrix/blob/master/fastq2matrix/vcf.py
fasta_output_file=${vcf_dir}/${study_accession}*snps.fa

# Parameters

# Commands


# Print arguments
printf "\n"
echo "Arguments:"
printf "\n"
echo ${study_accession}
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
# Nb needs datamash - conda install -c bioconda datamash
vcf2fasta.py --vcf ${vcf_input_file} --ref ${ref_fasta} --threads 20 --snps --snps-no-filter

# Move results to output directory
mv ${fasta_output_file} ${fasta_output_dir}
