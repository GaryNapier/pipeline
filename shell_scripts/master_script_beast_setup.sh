#! /bin/bash

# master_script_beast_setup.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Master script for transmission analysis up to the point of running Beast

# Arguments to script:

# Input:

# Main steps:

# Output:

# RUN
# Assume run from transmission/

# ------------------------------------------------------------------------------

# Setup - GLOBAL

cd ~/transmission

# Variables
gvcf_file_suffix=.g.vcf.gz
date_column=collection_date

# Directories
metadata_dir=~/metadata/
metadata_local_dir=metadata/
fasta_dir=fasta/
vcf_dir=~/vcf/
ref_dir=~/refgenome/
xml_script_location=beast_xml/
beast_results_dir=beast_results/
itol_annotation_dir=itol_annotations/
newick_output_dir=newick/
dist_and_pca_dir=dist_and_pca/
plots_dir=plots/

# Files
study_accession_list=$(cat ${metadata_local_dir}study_accession_list.txt)
metadata_file=${metadata_dir}tb_data_28_01_2021_clean.csv
xml_template_file=${xml_script_location}tb_template.xml
ref_fasta_file=${ref_dir}MTB-h37rv_asm19595v2-eg18.fa
ex_loci_file=${ref_dir}excluded_loci_rep_regions_dr_regions.bed

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

    if [ ! -f ${val_multi_vcf_file} ]; then

        echo "------------------------------------------------------------------------------"

        echo "Variant calling and concat VCFs of samples"
        printf "\n"
        echo "Running shell_scripts/variant_calling_and_concat_gvcfs.sh - outputs file ${val_multi_vcf_file}"
        set -x
        # shell_scripts/variant_calling_and_concat_gvcfs.sh <study_accession>   <metadata_file>  <vcf_dir>  <gvcf_file_suffix>  <ref_file>
        shell_scripts/variant_calling_and_concat_gvcfs.sh   ${study_accession}  ${metadata_file} ${vcf_dir} ${gvcf_file_suffix} ${ref_fasta_file}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"
    else
        echo "------------------------------------------------------------------------------"
        echo "File ${val_multi_vcf_file} already exists, skipping shell_scripts/variant_calling_and_concat_gvcfs.sh"
        echo "------------------------------------------------------------------------------"
        printf "\n"
    fi

    # ------------------------------------------------------------------------------

    # Variant filtering

    if [ ! -f ${filt_multi_vcf_file} ]; then

        echo "------------------------------------------------------------------------------"

        echo "Variant filtering"
        printf "\n"
        echo "Running shell_scripts/variant_filtering.sh - outputs file ${filt_multi_vcf_file}"
        printf "\n"
        set -x
        # shell_scripts/variant_filtering.sh    <multi-sample vcf file name>  <output file name>     <bed file name>  <ref fasta>
        shell_scripts/variant_filtering.sh      ${val_multi_vcf_file}         ${filt_multi_vcf_file} ${ex_loci_file}  ${ref_fasta_file}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"

    else
        echo "------------------------------------------------------------------------------"
        echo "File ${filt_multi_vcf_file} already exists, skipping shell_scripts/variant_filtering.sh"
        echo "------------------------------------------------------------------------------"
        printf "\n"
    fi

    # ------------------------------------------------------------------------------

    # Plink distances

    # Define output file for test of existence - these files are used by r_scripts/transmission_clusters.R

    dist_file=${dist_and_pca_dir}${study_accession}.dist.dist
    dist_id_file=${dist_and_pca_dir}${study_accession}.dist.dist.id

    if [ ! -f ${dist_file} ] || [ ! -f ${dist_id_file} ]; then

        echo "------------------------------------------------------------------------------"

        echo "Plink distances"
        printf "\n"

        echo "Running Plink distances command - outputs files to ${dist_and_pca_dir}"
        printf "\n"
        set -x
        # shell_scripts/plink_dist_and_pca.sh   <study_accession>   <vcf_input_file>        <output_dir>
        shell_scripts/plink_dist_and_pca.sh     ${study_accession}  ${filt_multi_vcf_file}  ${dist_and_pca_dir}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"
    else
        echo "------------------------------------------------------------------------------"
        echo "Files ${dist_file} and ${dist_id_file} already exist, skipping shell_scripts/plink_dist_and_pca.sh"
        echo "------------------------------------------------------------------------------"
        printf "\n"

    fi

    # ------------------------------------------------------------------------------

    # VCF to fasta

    # Define output file for test - files used by shell_scripts/iqtree.sh
    unfilt_fasta_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa

    if [ ! -f ${unfilt_fasta_file} ]; then

        echo "------------------------------------------------------------------------------"

        echo "VCF to fasta"
        printf "\n"

        echo "Running VCF to fasta command - outputs file ${unfilt_fasta_file}"
        set -x
        # Nb vcf2fasta.py (run by vcf2fasta.sh) needs datamash - conda install -c bioconda datamash
        # shell_scripts/vcf2fasta.sh <study_accession>  <vcf_file>              <fasta_output_dir>   <ref_fasta>
        shell_scripts/vcf2fasta.sh   ${study_accession} ${filt_multi_vcf_file}  ${fasta_dir}         ${ref_fasta_file}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"
    else
        echo "------------------------------------------------------------------------------"
        echo "File ${unfilt_fasta_file} exists, skipping shell_scripts/vcf2fasta.sh"
        echo "------------------------------------------------------------------------------"
        printf "\n"
    fi

    # ------------------------------------------------------------------------------

    # Transmission clusters

    if [ ! -f ${clusters_data_file} ]; then

        echo "------------------------------------------------------------------------------"

        echo "Transmission clusters"
        printf "\n"

        echo "Running transmission clusters - outputs file ${clusters_data_file} and plot in ${plots_dir}"
        set -x
        # Rscript r_scripts/transmission_clusters.R <study_accession>   <threshold> <dm_file>    <id_file>       <output dir>           <output dir for plot>
        Rscript r_scripts/transmission_clusters.R   ${study_accession}  100         ${dist_file} ${dist_id_file} ${metadata_local_dir}  ${plots_dir}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"
    else
        echo "------------------------------------------------------------------------------"
        echo "File ${clusters_data_file} exists, skipping r_scripts/transmission_clusters.R"
        echo "------------------------------------------------------------------------------"
        printf "\n"

    fi


    # ------------------------------------------------------------------------------

    # IQ tree

    # Define output files for test
    iqtree_file=${newick_output_dir}${study_accession}*iqtree
    treefile=${newick_output_dir}${study_accession}*treefile
    filt_fasta_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa.filt

    if [ ! -f ${iqtree_file} ] || [ ! -f ${treefile} ] || [ ! -f ${filt_fasta_file} ]; then

        echo "------------------------------------------------------------------------------"

        echo "IQ tree"
        printf "\n"

        echo "Running shell_scripts/iqtree.sh - outputs files ${iqtree_file}, ${treefile} and fasta/<study_accession>.filt.val.gt.g.snps.fa.filt"
        set -x
        # shell_scripts/iqtree.sh <study_accession>   <fasta_dir>           <cluster_file>                <newick_output_dir>
        shell_scripts/iqtree.sh   ${study_accession}  ${unfilt_fasta_file}  ${clusters_data_file}         ${newick_output_dir}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"
    else
        echo "------------------------------------------------------------------------------"
        echo "Files ${iqtree_file} and ${treefile} exist, skipping shell_scripts/iqtree.sh"
        echo "------------------------------------------------------------------------------"
        printf "\n"
    fi


    # ------------------------------------------------------------------------------

    # Create dated fasta file for Beast

    # Define out file
    dated_fasta_file=${study_accession}.dated.fa

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
        fasta_add_annotations.py --fasta ${filt_fasta_file} --csv ${metadata_file} --out ${dated_fasta_file} --annotations ${date_column}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"
    else
        echo "------------------------------------------------------------------------------"
        echo "File ${dated_fasta_file} exists, skipping fasta_add_annotations.py"
        echo "------------------------------------------------------------------------------"
        printf "\n"
    fi

    # ------------------------------------------------------------------------------

    # Subset fasta - chop dated_fasta_file into clusters outputted by r_scripts/transmission_clusters

    clust_fasta_dir=`dirname ${filt_fasta_file}`/
    clust_fasta_files=${clust_fasta_dir}${study_accession}.clust_*

    # n.b. weird/ugly syntax/code for checking if multiple fasta files exist because don't know how many there will be or even if the first cluster is '1'.
    # It's likely to be '1', but the clusters are arbitrarily numbered by the transmission_clusters.R script
    if ! ls ${clust_fasta_files} 1> /dev/null 2>&1; then

        echo "------------------------------------------------------------------------------"

        echo "Subset fasta"
        printf "\n"

        echo "Running shell_scripts/subset_fasta.sh - outputs ${clust_fasta_files}"
        set -x
        shell_scripts/subset_fasta.sh ${study_accession} ${clusters_data_file} ${dated_fasta_file}
        set +x
        echo "------------------------------------------------------------------------------"
        printf "\n"
    else
        echo "------------------------------------------------------------------------------"
        echo "File(s) ${clust_fasta_files} exist(s), skipping shell_scripts/subset_fasta.sh"
        echo "------------------------------------------------------------------------------"
        printf "\n"
    fi


# End for-loop on study accession numbers from study_accession_list
done


# NOW RUN shell_scripts/master_script_beast.sh





# Reset

# # Concat vcfs
# rm ${val_multi_vcf_file}
# rm ~/vcf/PRJEB7669.val.gt.g.vcf.gz
# rm ~/vcf/PRJEB7669*
# rm -r genomicsDB
# rm tmp/*
# rm logs/*
#
# # Variant filtering
# # rm ${filt_multi_vcf_file}
# rm ~/vcf/PRJEB7669.filt.val.gt.g.vcf.gz
#
# # Plink distances
# # rm ${dist_and_pca_dir}${study_accession}*
# rm dist_and_pca/PRJEB7669.*
#
# # VCF to fasta
# # rm ${unfilt_fasta_file}
# rm fasta/PRJEB7669.filt.val.gt.g.snps.fa
#
# # Transmission clusters
# # rm ${clusters_data_file}
# rm metadata/PRJEB7669.clusters
# rm plots/PRJEB7669*
#
# # IQ tree
# # rm ${iqtree_file} ${treefile} ${filt_fasta_file}
# rm newick/PRJEB7669*
# rm fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt
# rm fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt.ckp.gz
#
# # ITOL annotations
# # rm ${clust_itol_file} ${dr_itol_file} ${lin_itol_file}
# rm itol_annotations/PRJEB7669.*
#
# # Beast
#
# # Dated fasta file
# # rm ${dated_fasta_file}
# rm fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt.dated.fa
#
# # Create XML file for Beast
# # rm ${xml_file}
# rm beast_xml/PRJEB7669*
#
# # Run Beast
# # rm ${beast_log_file} ${beast_trees_file} ${beast_state_file}
# rm beast_results/PRJEB7669.*
