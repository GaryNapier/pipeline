#! /bin/bash

# FILENAME
# vcf_links.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# INFORMATION AND DESCRIPTION
# Make links from remote database (DB) of vcf files and put in local directory
# Pulls in both the standard vcf file and the 'validated' file (.validated) and the index file (.tbi)

# Arguments to script:
# samples_for_vcf_file - file with list of samples
# vcf_remote - location of remote vcf DB
# vcf_db - local vcf DB/directory
# gvcf_file_suffix - suffix of vcf files; default .g.vcf.gz

# Input:
# samples_for_vcf_file
# e.g. 
# $ cat -A <file>
# ERR2864254$
# ERR760749$
# SRR5073823$
# ERR211998$
# ERR2514005$

# Main steps:
# Loop over the samples in the file
# Test if the file exists
# If not exist - make the link from the remote DB to local for the three file types

# Output:
# Links to remote vcfs in local vcf dir

# RUN EXAMPLE
# vcf_links.sh <samples_for_vcf_file> <vcf_remote> <vcf_db> <gvcf_file_suffix>


# ------------------------------------------------------------------------------

# Setup

# Variables
gvcf_file_suffix=${4:-.g.vcf.gz}

# Directories
vcf_remote=${2}
vcf_db=${3}

# Files and suffixes/prefixes
samples_for_vcf_file=${1}

# Parameters


# Print all args and params/variables
printf "\n"
echo "Parameters:"
printf "\n"
echo "samples_for_vcf_file - file with list of samples"
echo "vcf_remote - location of remote vcf DB"
echo "vcf_db - local vcf DB/directory"
echo "gvcf_file_suffix - suffix of vcf files; default .g.vcf.gz"
printf "\n"

# Make sure the vcf samples are linked from the remote vcf DB to the local DB - create symbolic links from remote DB
for samp in $(cat ${samples_for_vcf_file}); do
file=${vcf_remote}${samp}${gvcf_file_suffix}
# Test if link to local vcf DB exists
if [ ! -f ${vcf_db}${samp}${gvcf_file_suffix} ]; then
echo "MAKING LINK FROM ${file} TO ${vcf_db}${samp}${gvcf_file_suffix}"
echo ${file} | parallel -j 1 --bar ln -s {} ${vcf_db}
echo ${file}.validated | parallel -j 1 --bar ln -s {} ${vcf_db}
echo ${file}.tbi | parallel -j 1 --bar ln -s {} ${vcf_db}
fi
done