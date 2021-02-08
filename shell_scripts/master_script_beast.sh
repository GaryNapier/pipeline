#! /bin/bash

# master_script.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Master script for transmission analysis - run Beast

# Arguments to script:
# None

# Input:
# None

# Main steps:

# Output:

# RUN
# Assume run from transmission/
# shell_scripts/master_script_beast.sh

# ------------------------------------------------------------------------------

# Setup - GLOBAL

cd ~/transmission

# Directories
metadata_dir=~/metadata/
metadata_local_dir=metadata/
fasta_dir=fasta/
xml_script_location=beast_xml/
beast_results_dir=beast_results/

# Files
metadata_file=${metadata_dir}tb_data_collated_28_10_2020_clean.csv
xml_template_file=${xml_script_location}tb_template.xml

# Variables
date_column=collection_date
study_accession_list=$(cat ${metadata_local_dir}study_accession_list.txt)

# Parameters

# Commands

# Print arguments

printf "\n"
echo "Processing study accessions:"
cat ${study_accession_list}
printf "\n"

# Beast loop

for study_accession in ${study_accession_list}; do

    # Define cluster numbers from clusters file (result of transmission_clusters.R)
    clusters_file=${metadata_dir}${study_accession}.clusters
    cluster_nums=$(cat ${clusters_file} | tail -n +2 | cut -f2 | sort -n | uniq)

    # Loop over cluster numbers
    for clust_num in ${cluster_nums}; do

        # Create XML file for Beast

        # Define files
        # In:
        dated_fasta_file=${fasta_dir}${study_accession}.clust_${clust_num}.fa
        # Out:
        xml_file=${xml_script_location}${study_accession}.clust_${clust_num}.xml


        if [ ! -f ${xml_file} ]; then

            echo "------------------------------------------------------------------------------"

            echo "Create Beast XML file from template"
            printf "\n"

            echo "Running beast_xml.R command - outputs ${xml_file}"
            set -x
            Rscript r_scripts/beast_xml.R -s ${study_accession} -t ${xml_template_file} -f ${dated_fasta_file} -o ${xml_file}
            set +x
            echo "------------------------------------------------------------------------------"
            printf "\n"
        else
            echo "------------------------------------------------------------------------------"
            echo "File ${xml_file} exists, skipping r_scripts/beast_xml.R"
            echo "------------------------------------------------------------------------------"
            printf "\n"
        fi

        # ------------------------------------------------------------------------------

        # Run Beast

        # Beast output files
        beast_log_file=${beast_results_dir}${study_accession}.clust_${clust_num}.log
        beast_trees_file=${beast_results_dir}${study_accession}.clust_${clust_num}.trees
        beast_state_file=${beast_results_dir}${study_accession}.clust_${clust_num}.xml.state

        if [ ! -f ${beast_log_file} ] || [ ! -f ${beast_trees_file} ] || [ ! -f ${beast_state_file} ];then

            echo "------------------------------------------------------------------------------"

            echo "Run Beast"
            printf "\n"

            echo "Running BEAST - outputs ${beast_log_file}, ${beast_trees_file} and ${beast_state_file}"
            set -x
            # Run Beast and clean up - put in right folder
            beast ${xml_file} && mv ${study_accession}.clust_${clust_num}.log ${study_accession}.clust_${clust_num}.trees ${study_accession}.clust_${clust_num}.xml.state ${beast_results_dir}
            set +x
            echo "------------------------------------------------------------------------------"
            printf "\n"
        else
            echo "------------------------------------------------------------------------------"
            echo "Files ${beast_log_file}, ${beast_trees_file} and ${beast_state_file} exist, skipping beast"
            echo "------------------------------------------------------------------------------"
            printf "\n"

        fi


    # End Beast loop on cluster numbers
    done
# End study accession loop
done
