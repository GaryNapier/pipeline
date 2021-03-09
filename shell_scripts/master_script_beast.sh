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
metadata_file=${metadata_dir}tb_data_18_02_2021_clean.csv
# metadata_file=thailand_test/thailand_test_metadata.csv
xml_template_file=${xml_script_location}tb_template.xml

# Variables
date_column=year
study_accession_list=$(cat ${metadata_local_dir}study_accession_list.txt)

# Parameters

# Commands

# Print arguments

printf "\n"
echo "Processing study accessions:"
echo ${study_accession_list}
printf "\n"

# Beast loop

for study_accession in ${study_accession_list}; do


    # Subset fasta - chop fasta file into the clusters outputted by r_scripts/transmission_clusters.R

    # In:
    clusters_data_file=${metadata_local_dir}${study_accession}.clusters
    fasta_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa
    # Out:
    clust_fasta_files=${fasta_dir}${study_accession}.clust_*

    # n.b. weird/ugly syntax/code for checking if multiple fasta files exist because don't know how many there will be or even if the first cluster is '1'.
    # It's likely to be '1', but the clusters are arbitrarily numbered by the transmission_clusters.R script
    if ! ls ${clust_fasta_files} 1> /dev/null 2>&1; then

        echo "------------------------------------------------------------------------------"

        echo "Subset fasta"
        printf "\n"

        echo "Running shell_scripts/subset_fasta.sh - outputs ${clust_fasta_files}"
        set -x
        # shell_scripts/subset_fasta.sh <study_accession>   <metadata_file>  <fasta_file>
        shell_scripts/subset_fasta.sh   ${study_accession}  ${metadata_file} ${fasta_file}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"
    else
        echo "------------------------------------------------------------------------------"
        echo "File(s) ${clust_fasta_files} exist(s), skipping shell_scripts/subset_fasta.sh"
        echo "------------------------------------------------------------------------------"
        printf "\n"
    fi

    # Define cluster numbers from metadata
    cluster_nums$(cat ${metadata_file} | csvtk grep -f study_accession_word -p ${study_accession} | csvtk cut -f major_lineage | tail -n +2 | sort | uniq)

    # Loop over cluster numbers
    for clust_num in ${cluster_nums}; do

        # ------------------------------------------------------------------------------

        # Create dated fasta file for Beast

        # In:
        clust_fasta_file=${fasta_dir}${study_accession}.clust_${clust_num}.fa
        # Out:
        dated_fasta_file=${fasta_dir}${study_accession}.clust_${clust_num}.dated.fa

        if [ ! -f ${dated_fasta_file} ]; then

            echo "------------------------------------------------------------------------------"

            echo "Dated fasta file"
            printf "\n"

            echo "Running fasta_add_annotations.py - outputs ${dated_fasta_file}"
            printf "\n"

            # Run fasta_add_annotations.py to concatenate the metadata collection date with the sample ID in the fasta file
            # e.g. >ERR1234 -> >ERR1234_01_01_10
            # - see https://github.com/pathogenseq/pathogenseq-scripts/tree/master/scripts
            set -x
            fasta_add_annotations.py --fasta ${clust_fasta_file} --csv ${metadata_file} --out ${dated_fasta_file} --annotations ${date_column} --id-key wgs_id
            set +x
            echo "------------------------------------------------------------------------------"
            printf "\n"
        else
            echo "------------------------------------------------------------------------------"
            echo "File ${dated_fasta_file} exists, skipping fasta_add_annotations.py"
            echo "------------------------------------------------------------------------------"
            printf "\n"
        fi


        # Create XML file for Beast

        # Define files
        # In:
        dated_fasta_file=${fasta_dir}${study_accession}.clust_${clust_num}.dated.fa
        # Out:
        # N.b. need the underscore after the study accession, not "." because Beast seems to split on the first ".", so will lose the cluster number
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

        # Beast output files - n.b. saves the .log and .trees files to the current directory and ONLY takes the study_accession
        # So files are ./<study_accession>.log and ./<study_accession>.trees
        beast_log_file=${study_accession}.log
        beast_trees_file=${study_accession}.trees
        beast_state_file=${study_accession}.clust_${clust_num}.xml.state

        if [ ! -f ${beast_results_dir}${study_accession}.clust_${clust_num}.log ] || [ ! -f ${beast_results_dir}${study_accession}.clust_${clust_num}.trees ] || [ ! -f ${beast_results_dir}${beast_state_file} ];then

            echo "------------------------------------------------------------------------------"

            echo "Run Beast"
            printf "\n"

            echo "Running BEAST - outputs ${beast_log_file}, ${beast_trees_file} and ${beast_state_file}"
            set -x
            # Run Beast and clean up - put in right folder
            beast -threads 20 ${xml_file}
            mv ${beast_log_file} ${beast_results_dir}${study_accession}.clust_${clust_num}.log
            mv ${beast_trees_file} ${beast_results_dir}${study_accession}.clust_${clust_num}.trees
            mv ${beast_state_file} ${beast_results_dir}
            set +x
            echo "------------------------------------------------------------------------------"
            printf "\n"
        else
            echo "------------------------------------------------------------------------------"
            echo "Files ${beast_results_dir}${study_accession}.clust_${clust_num}.log, ${beast_results_dir}${study_accession}.clust_${clust_num}.trees and ${beast_results_dir}${beast_state_file} exist, skipping beast"
            echo "------------------------------------------------------------------------------"
            printf "\n"

        fi

        # ------------------------------------------------------------------------------

        # Run TreeAnnotator to get MCC tree (input to TransPhylo)

        # In
        ${treeann_in}=${beast_results_dir}${study_accession}.clust_${clust_num}.trees
        # Out
        ${treeann_out}=${beast_results_dir}${study_accession}.clust_${clust_num}.mcc.tree

        if [ ! -f ${treeann_out} ]; then

            echo "------------------------------------------------------------------------------"

            echo "Running TreeAnnotator on ${treeann_in} - outputs ${treeann_out}"
            set -x
            treeannotator -heights ca -burnin 10 ${treeann_in} ${treeann_out}
            set +x
            echo "------------------------------------------------------------------------------"
            printf "\n"
        else
            echo "------------------------------------------------------------------------------"
            echo "Files ${treeann_out} exists, skipping TreeAnnotator"
            echo "------------------------------------------------------------------------------"
            printf "\n"

        fi

    # End loop on cluster numbers
    done
# End study accession loop
done
