#! /bin/bash

# concat_filt_vcfs.sh

set -e
set -u
set -o pipefail

# Concatenate and filter vcfs of samples specific to project ready for distance matrix calculation
# vcfs must be pre-processed up to gvcf stage

# Input:
# List of samples, project code

# Steps:
# Concat vcfs corresponding to samples > filter vcfs

# Output:
# Concatenated and filtered vcf file

# RUN
# concat_filt_vcfs.sh PRJXXX ~/metadata/PRJXXX_samples.txt .g.vcf.gz.validated

# Variables
project_code=${1}

# Files
sample_list=$(cat ${2})
input_vcf_file_suffix=${3}
output_vcf_file=${project_code}.val.concat.filt.no_indels.vcf.gz


bcftools concat ${sample_list}${input_vcf_file_suffix}
