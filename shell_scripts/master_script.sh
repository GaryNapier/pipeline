#! /bin/bash

# master_script.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Master script for transmission analysis

# Arguments to script:

# Input:

# Main steps:

# Output:

# RUN
# Assume run from

# ------------------------------------------------------------------------------



# Setup - GLOBAL

# Variables


# Directories
metadata_dir=../metadata/
fasta_dir=fasta/
xml_script_location=beast_xml/
beast_results_dir=beast_results/
date_column=collection_date
clusters_data_file_dir=metadata/
itol_annotation_dir=itol_annotations/
newick_output_dir=newick/

# Files
metadata_file=${metadata_dir}tb_data_collated_28_10_2020_clean.csv
xml_template_file=${xml_script_location}tb_template.xml

# Parameters

# Commands

# Print arguments
printf ''
echo 'Arguments:'
printf ''
echo ${}
printf ''


for study_accession in ${study_accession_list}; do

    # Setup - LOOP (i.e. anything that depends on study_accession)


    # Files
    xml_file=${xml_script_location}${study_accession}.xml
    fasta_input_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa.filt
    dated_fasta_file=${fasta_input_file}.dated.fa

    echo "Running ${study_accession}"


    # ------------------------------------------------------------------------------

    # Plink distances

    shell_scripts/plink_dist_and_pca.sh


    # ------------------------------------------------------------------------------


    # VCF to fasta

    shell_scripts/vcf2fasta.sh

    # ------------------------------------------------------------------------------

    # transmission clusters

    Rscript r_scripts/transmission_clusters.R

    # ------------------------------------------------------------------------------

    # IQ tree

    # ./shell_scripts/iqtree.sh <study_accession>   <fasta_dir>       <metadata_dir>        <newick_output_dir>
    ./shell_scripts/iqtree.sh   ${study_accession}  ${fasta_dir}      ${metadata_dir}       ${newick_output_dir}

    # ------------------------------------------------------------------------------

    # ITOL annotations
    # Rscript r_scripts/itol_annotation.R   <study_accession>   <metadata_file>     <clusters_data_file_dir>                 <itol_annotation_dir>
    Rscript r_scripts/itol_annotation.R     ${study_accession}  ${metadata_file}    ${clusters_data_file_dir}                ${itol_annotation_dir}

    # ------------------------------------------------------------------------------

    # Beast

    echo "Running BEAST"

    # Run fasta_add_annotations.py to concatenate the metadata collection date with the sample ID in the fasta file

    # e.g. >ERR1234 -> >ERR1234_01_01_10
    # - see https://github.com/pathogenseq/pathogenseq-scripts/tree/master/scripts
    # fasta_add_annotations.py --fasta fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt --csv ../metadata/tb_data_collated_28_10_2020_clean.csv --out fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt.dated.fa --annotations collection_date
    fasta_add_annotations.py --fasta ${fasta_input_file} --csv ${metadata_file} --out ${dated_fasta_file} --annotations ${date_column}


    # Run R script to fill out xml template file
    Rscript r_scripts/beast_xml.R -s ${study_accession} -t beast_xml/tb_template.xml -f ${fasta_dir} -o ${xml_script_location}


    # Run beast and clean up
    beast ${xml_file} && mv ${study_accession}.log ${study_accession}.trees ${study_accession}.xml.state ${beast_results_dir}


done
