#! /bin/bash

# FILENAME
# tree_pipeline.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# INFORMATION AND DESCRIPTION
# Make tree (newick file) from list of samples file

# Arguments to script:
# project_code - string to prefix files, e.g. 'PRJ123' or 'isoniazid'
# threads - number of threads to process variant calling; default 20
# vcf_db - location of vcf files for individual samples
# samples_for_vcf_file - file of sample names - one column, no header
# refgenome_fasta - full path to ref genome fasta file; default ~/refgenome/MTB-h37rv_asm19595v2-eg18.fa
# vcf_output_dir - location to put multi-sample vcf file
# exclude_loci_bed_file - full path to bed file of excluded loci; default ~/refgenome/excluded_loci.bed
# fasta_dir - location of fasta input/output file(s)
# newick_dir - location to store tree results

# Main steps:
# variant_calling_and_concat_gvcfs.sh - input a file of samples, pull the SNPs from their individual VCF files (variant call) and put together into one mult-sample vcf file
# variant_filtering.sh - filter the variants in the multi-sample vcf file
# vcf2fasta.sh - convert multi-sample vcf to fasta file
# iqtree.sh - run tree on multi-sample fasta

# Output:
# newick file - tree of input samples

# RUN EXAMPLE
# tree_pipeline.sh 1              2        3                4                      5           6            7                 8                       9            10
# tree_pipeline.sh <project_code> <vcf_db> <vcf_output_dir> <samples_for_vcf_file> <fasta_dir> <newick_dir> <refgenome_fasta> <exclude_loci_bed_file> <vcf_remote> <threads>

# ------------------------------------------------------------------------------

# Setup

# ARG1=${1:-foo}
# ARG2=${2:-bar}
# ARG3=${3:-1}
# ARG4=${4:-$(date)}

# echo "$ARG1"
# echo "$ARG2"
# echo "$ARG3"
# echo "$ARG4"

# Variables
project_code=${1}
threads=${10:-20}

# Directories
vcf_db=${2}
vcf_output_dir=${3}
vcf_remote=${9:-/mnt/storage7/jody/tb_ena/per_sample/}
fasta_dir=${5}
newick_dir=${6}

# Files and suffixes/prefixes

# variant_calling_and_concat_gvcfs.sh
samples_for_vcf_file=${4}
gvcf_file_suffix=.g.vcf.gz 
refgenome_fasta=${7:-~/refgenome/MTB-h37rv_asm19595v2-eg18.fa}

# variant_filtering.sh
multi_samp_vcf=${vcf_output_dir}${project_code}.val.gt.g.vcf.gz
filt_multi_samp_vcf_file=${vcf_output_dir}${project_code}.filt.val.gt.g.vcf.gz
exclude_loci_bed_file=${8:-~/refgenome/excluded_loci.bed}

# vcf2fasta.sh

# iqtree.sh
fasta_file=${fasta_dir}/${project_code}.filt.val.gt.g.snps.fa

# Print all args and params/variables
printf "\n"
echo "Parameters:"
printf "\n"
echo "project code: " ${project_code}
echo "threads: " ${threads}
echo "vcf DB: " ${vcf_db}
echo "samples file: " ${samples_for_vcf_file}
echo "ref genome fasta file: " ${refgenome_fasta}
echo "vcf output directory: " ${vcf_output_dir}
echo "excluded loci bed file: " ${exclude_loci_bed_file}
echo "fasta directory: " ${fasta_dir}
echo "newick directory: " ${newick_dir}
printf "\n"

# Do joint variant calling on validated gvcfs of samples specific to study accession ready for distance matrix calculation.
variant_calling_and_concat_gvcfs.sh ${project_code} ${samples_for_vcf_file} ${vcf_db} ${gvcf_file_suffix} ${refgenome_fasta} ${threads} ${vcf_output_dir} ${vcf_remote}

# Do variant filtering (TB specific)
variant_filtering.sh ${multi_samp_vcf} ${filt_multi_samp_vcf_file} ${exclude_loci_bed_file} ${refgenome_fasta}

# Convert filtered multi-vcf to fasta, retaining just snps 
vcf2fasta.sh ${project_code} ${filt_multi_samp_vcf_file} ${fasta_dir} ${refgenome_fasta}

# Run IQtree on fasta
# iqtree.sh ${project_code} ${fasta_file} ${newick_dir}
FastTreeMP_double -gtr -nt ${fasta_file} > ${newick_dir}${project_code}.treefile

echo "---"
echo "tree_pipeline.sh DONE"
echo "---"
