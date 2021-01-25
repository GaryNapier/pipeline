#! /bin/bash

# beast.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Description:
# Run python script to concatenate sample IDs in fasta with dates from metadata.
# Run python script for Beauti to specify parameters and generate XML file for beast.
# Run Beast and generate 3 Beast files.

# Arguments to script:
# study accession, fasta file location, metadata file name, name of date column in metadata file, xml script location, beast output location

# Input:
# fasta file of cluster from study accession, metadata file with dates for samples

# Main steps:
# Run fasta_add_annotations.py to concatenate the date from the metadata file with sample ID
# Script is from - https://github.com/pathogenseq/pathogenseq-scripts/tree/master/scripts
# Outputs new 'dated' fasta file

# Run python script for generating XML file for Beast with new fasta file from beast2-xml.py script from beast2-xml package
# See https://github.com/acorg/beast2-xml
# Outputs XML file for Beast

# Output:
# Dated fasta file
# XML file for beast
# Beast output files:
#   - <study_accession>.log
#   - <study_accession>.trees
#   - <study_accession>.xml.state

# RUN
# Assume run from transmission/
# shell_scripts/beast.sh <study accession> <fasta_file_location> <metadata_file_name>                              <date_column>   <xml_script_location> <beast_output_location>
# shell_scripts/beast.sh PRJEB7669        fasta/                 ../metadata/tb_data_collated_28_10_2020_clean.csv collection_date beast_xml/             beast_results/ 25 23 6 23

# Notes:
# See also Beauti/Beast tutorial: https://taming-the-beast.org/tutorials/Prior-selection/Prior-selection.pdf

# ------------------------------------------------------------------------------

# Setup

# Variables
study_accession=${1?Error: input study accession number e.g. PRJEB7669}

# MUTATIONRATE="${7:-$'1.0'}"
MUTATIONRATE="${MUTATIONRATE:-$'1.0'}"
# CLOCKRATE
# POPSIZELOWER
# POPSIZEUPPER
# POPSIZE
# FREQPARAMETER

# Directories
fasta_dir=${2?Error: input fasta file location}
xml_script_location=${5?Error input file location for XML output file}
beast_output_location=${6?Error input file location for Beast output files}

# Files
metadata_file=${3?Error: input metadata file name}
fasta_input_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa.filt
dated_fasta_file=${fasta_input_file}.dated.fa
xml_file=${xml_script_location}${study_accession}.xml

# Parameters
date_column=${4?Error input date column name from metadata}

# Commands

# Print arguments
printf '\n'
echo 'Arguments:'
printf '\n'
echo "Study accession: ${study_accession}"
echo "Fasta file:  ${fasta_input_file}"
echo "Dated fasta file: ${dated_fasta_file}"
echo "Metadata file: ${metadata_file}"
echo "Date column name: ${date_column}"
echo "XML file: ${xml_file}"
echo "Beast output files location: ${beast_output_location}"
echo "MUTATIONRATE: ${MUTATIONRATE}"


printf '\n'

# ------------------------------------------------------------------------------

# Run fasta_add_annotations.py to concatenate the metadata collection date with the sample ID in the fasta file

# e.g. >ERR1234 -> >ERR1234_01_01_10
# - see https://github.com/pathogenseq/pathogenseq-scripts/tree/master/scripts
# fasta_add_annotations.py --fasta fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt --csv ../metadata/tb_data_collated_28_10_2020_clean.csv --out fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt.dated.fa --annotations collection_date
fasta_add_annotations.py --fasta ${fasta_input_file} --csv ${metadata_file} --out ${dated_fasta_file} --annotations ${date_column}



# ------------------------------------------------------------------------------

# Generate XML file using beast2-xml.py (replaces Beauti GUI) - https://github.com/acorg/beast2-xml

# beast2-xml.py [-h] [--clockModel MODEL | --templateFile FILENAME]
#                      [--chainLength LENGTH] [--age ID=N [ID=N ...]]
#                      [--defaultAge N] [--dateUnit UNIT]
#                      [--dateDirection DIRECTION]
#                      [--logFileBasename BASE-FILENAME] [--traceLogEvery N]
#                      [--treeLogEvery N] [--screenLogEvery N] [--mimicBEAUti]
#                      [--sequenceIdDateRegex REGEX]
#                      [--sequenceIdDateRegexMayNotMatch] [--fastaFile FILENAME]
#                      [--readClass CLASSNAME] [--fasta | --fastq | --fasta-ss]
#
#
# beast2-xml.py --chainLength 100000000 --fastaFile ${dated_fasta_file}


# ------------------------------------------------------------------------------

# Run beast and clean up

# beast ${xml_file} && mv ${study_accession}.log ${study_accession}.trees ${study_accession}.xml.state ${beast_results_dir}


# ------------------------------------------------------------------------------
