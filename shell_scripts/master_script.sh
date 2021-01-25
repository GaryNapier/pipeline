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

cd ~/transmission

# Variables
gvcf_file_suffix=.g.vcf.gz

# Directories
metadata_dir=~/metadata/
metadata_local_dir=metadata/
fasta_dir=fasta/
vcf_dir=vcf/
ref_dir=~/refgenome/
xml_script_location=beast_xml/
beast_results_dir=beast_results/
date_column=collection_date
itol_annotation_dir=itol_annotations/
newick_output_dir=newick/
dist_and_pca_dir=dist_and_pca/
plots_dir=plots/

# Files
study_accession_list=$(cat ${metadata_local_dir}study_accession_list.txt)
metadata_file=${metadata_dir}tb_data_collated_28_10_2020_clean.txt
xml_template_file=${xml_script_location}tb_template.xml
ref_fasta_file=${ref_dir}MTB-h37rv_asm19595v2-eg18.fa
ex_loci_file=${ref_dir}excluded_loci.bed

# Parameters

# Commands

# Print arguments

printf "\n"
echo "Processing study accessions:"
echo ${study_accession_list}
printf "\n"

for study_accession in ${study_accession_list}; do

    echo "Running ${study_accession}"
    printf "\n"

    # Setup - LOOP (i.e. anything that depends on study_accession)

    # Files - define file names unless already created

    echo "Files for multiple scripts:"
    clusters_data_file=${metadata_local_dir}${study_accession}.clusters
    val_multi_vcf_file=${vcf_dir}${study_accession}.val.gt.g.vcf.gz
    filt_multi_vcf_file=${vcf_dir}${study_accession}.filt.val.gt.g.vcf.gz

    printf "${clusters_data_file} \n ${val_multi_vcf_file} \n ${filt_multi_vcf_file} \n"
    printf "\n"

    # ------------------------------------------------------------------------------

    # Concat VCFs

    echo "Concat VCF command:"

    echo "shell_scripts/variant_calling_and_concat_gvcfs.sh   ${study_accession}  ${metadata_file} ${vcf_dir} ${gvcf_file_suffix} ${ref_fasta_file}"
    printf "\n"

    echo "Running concat vcfs"

    # shell_scripts/variant_calling_and_concat_gvcfs.sh <study_accession>   <metadata_file>  <vcf_dir>  <gvcf_file_suffix>  <ref_file>
    shell_scripts/variant_calling_and_concat_gvcfs.sh   ${study_accession}  ${metadata_file} ${vcf_dir} ${gvcf_file_suffix} ${ref_fasta_file}

    # ------------------------------------------------------------------------------

    # # Variant filtering
    #
    # echo "Variant filtering command"
    # echo "shell_scripts/variant_filtering.sh   ${val_multi_vcf_file}              ${filt_multi_vcf_file} ${ex_loci_file}  ${ref_fasta_file}"
    #
    # echo "Running variant filtering"
    #
    # # shell_scripts/variant_filtering.sh <multi-sample vcf file name>       <output file name>     <bed file name>  <ref fasta>
    # shell_scripts/variant_filtering.sh   ${val_multi_vcf_file}              ${filt_multi_vcf_file} ${ex_loci_file}  ${ref_fasta_file}
    #
    # # ------------------------------------------------------------------------------
    #
    # # Plink distances
    #
    # echo "Plink distances command"
    #
    # echo "shell_scripts/plink_dist_and_pca.sh     ${study_accession}  ${filt_multi_vcf_file}  ${dist_and_pca_dir}"
    #
    # echo "Running Plink distances"
    #
    # # shell_scripts/plink_dist_and_pca.sh   <study_accession>   <vcf_input_file>        <output_dir>
    # shell_scripts/plink_dist_and_pca.sh     ${study_accession}  ${filt_multi_vcf_file}  ${dist_and_pca_dir}
    #
    #
    # # ------------------------------------------------------------------------------
    #
    # # VCF to fasta
    #
    # echo "VCF to fasta command"
    # echo "shell_scripts/vcf2fasta.sh   ${study_accession} ${filt_multi_vcf_file}  ${fasta_dir}         ${ref_fasta_file}"
    #
    # echo "Running VCF to fasta"
    #
    # # shell_scripts/vcf2fasta.sh <study_accession>  <vcf_file>              <fasta_output_dir>   <ref_fasta>
    # shell_scripts/vcf2fasta.sh   ${study_accession} ${filt_multi_vcf_file}  ${fasta_dir}         ${ref_fasta_file}
    #
    # # ------------------------------------------------------------------------------
    #
    # # Transmission clusters
    #
    # # Transmission clusters files
    # dist_file=${dist_and_pca_dir}${study_accession}.dist.dist
    # dist_id_file=${dist_and_pca_dir}${study_accession}.dist.dist.id
    #
    # echo "Transmission clusters command"
    # echo "Rscript r_scripts/transmission_clusters.R   ${study_accession}  50          ${dist_file} ${dist_id_file} ${metadata_local_dir}  ${plots_dir}"
    #
    # echo "Running transmission clusters"
    #
    # # Rscript r_scripts/transmission_clusters.R <study_accession>   <threshold> <dm_file>    <id_file>       <output dir>           <output dir for plot>
    # Rscript r_scripts/transmission_clusters.R   ${study_accession}  50          ${dist_file} ${dist_id_file} ${metadata_local_dir}  ${plots_dir}
    #
    # # ------------------------------------------------------------------------------
    #
    # # IQ tree
    #
    # # IQtree files
    # iq_fasta_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa
    #
    # echo "IQ tree command"
    # echo "shell_scripts/iqtree.sh   ${study_accession}  ${iq_fasta_file}      ${clusters_data_file}         ${newick_output_dir}"
    #
    # echo "Running IQ tree"
    #
    # # shell_scripts/iqtree.sh <study_accession>   <fasta_dir>           <cluster_file>                <newick_output_dir>
    # shell_scripts/iqtree.sh   ${study_accession}  ${iq_fasta_file}      ${clusters_data_file}         ${newick_output_dir}
    #
    #
    # # ------------------------------------------------------------------------------
    #
    # # ITOL annotations
    #
    # echo "ITOL annotations command"
    # echo "Rscript r_scripts/itol_annotation.R     ${study_accession}  ${metadata_file}    ${clusters_data_file}                ${itol_annotation_dir}"
    #
    # echo "Running ITOL annotations"
    #
    # Rscript r_scripts/itol_annotation.R     ${study_accession}  ${metadata_file}    ${clusters_data_file}                ${itol_annotation_dir}
    #
    # # ------------------------------------------------------------------------------
    #
    # # Beast
    #
    # # Beast files
    #
    # xml_file=${xml_script_location}${study_accession}.xml
    # beast_fasta_input_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa.filt
    # dated_fasta_file=${fasta_input_file}.dated.fa
    #
    # # ---
    #
    # echo "fasta_add_annotations.py command"
    # echo "fasta_add_annotations.py    --fasta ${beast_fasta_input_file} --csv ${metadata_file} --out ${dated_fasta_file} --annotations ${date_column}"
    #
    # echo "Running fasta_add_annotations.py"
    #
    # # Run fasta_add_annotations.py to concatenate the metadata collection date with the sample ID in the fasta file
    # # e.g. >ERR1234 -> >ERR1234_01_01_10
    # # - see https://github.com/pathogenseq/pathogenseq-scripts/tree/master/scripts
    # fasta_add_annotations.py    --fasta ${beast_fasta_input_file} --csv ${metadata_file} --out ${dated_fasta_file} --annotations ${date_column}
    #
    # # ---
    #
    # echo "beast_xml.R command"
    # echo "Rscript r_scripts/beast_xml.R -s ${study_accession} -t ${xml_template_file} -f ${dated_fasta_file} -o ${xml_script_location}"
    #
    # echo "Running beast_xml.R"
    #
    # # Run R script to fill out xml template file
    # Rscript r_scripts/beast_xml.R -s ${study_accession} -t ${xml_template_file} -f ${dated_fasta_file} -o ${xml_script_location}
    #
    # # ---
    #
    # echo "BEAST command"
    # "beast ${xml_file} && mv ${study_accession}.log ${study_accession}.trees ${study_accession}.xml.state ${beast_results_dir}"
    #
    # # Run beast and clean up
    #
    # echo "Running BEAST"
    # beast ${xml_file} && mv ${study_accession}.log ${study_accession}.trees ${study_accession}.xml.state ${beast_results_dir}
    #
    # # ------------------------------------------------------------------------------
    #
    # # Transphylo
    #
    # echo "Transphylo command"
    #
    # echo "Running Transphylo"


done
