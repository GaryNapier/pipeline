#! /bin/bash

# variant_calling_and_concat_gvcfs.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Do joint variant calling on validated gvcfs of samples specific to study accession ready for distance matrix calculation.
# vcfs must be pre-processed up to gvcf stage

# Arguments to script:
# study accession number (from metadata file column), metadata file, gvcf directory location, gvcf file suffix from ready gvcfs, reference fasta

# Input:
# study accession number (from metadata file column), metadata file, gvcf directory location, gvcf file suffix from ready gvcfs, reference fasta, gvcf files

# Main steps:
# Check if validated gvcf files exist for each sample and add to map if exist > GenomicsDBImport (import gvcfs to database) >
# GenotypeGVCFs (joint variant calling) > bcftool concat (concatenate (rbind) the split genome into one vcf) > index merged file

# Output:
# Concatenated and filtered vcf file

# RUN
# Assume run from ~/transmission :
# shell_scripts/variant_calling_and_concat_gvcfs.sh <study_accession> <metadata_file>                                  <vcf_dir>  <gvcf_file_suffix>  <ref_file>
# shell_scripts/variant_calling_and_concat_gvcfs.sh PRJEB7669         ~/metadata/tb_data_collated_28_10_2020_clean.csv ~/vcf/     .g.vcf.gz           ~/refgenome/MTB-h37rv_asm19595v2-eg18.fa

# ------------------------------------------------------------------------------

# Setup

# Variables
nl="\n"
study_accession=${1}

# Directories
vcf_dir=${3}
tmp_dir=tmp/
logs_dir=logs/
genomicsDB_dir=genomicsDB/

# Files and suffixes/prefixes
metadata_file=${2}
gvcf_file_suffix=${4}
val_gvcf_file_suffix=${gvcf_file_suffix}.validated
ref_file=${5}
ref_index=${ref_file}.fai
sample_list_file=tmp/sample_list
vcf_map_file=${tmp_dir}vcf_map_file
failed_samples_file=${logs_dir}variant_calling_and_concat_gvcfs_failed_samples
windows_file=${tmp_dir}windows_file
genotyped_vcf_suffix=.val.gt${gvcf_file_suffix}
output_vcf_file=${vcf_dir}${study_accession}${genotyped_vcf_suffix}

# Parameters
threads=4
num_genome_chunks=20

printf "\n"
echo "Parameters:"
printf "\n"
echo ${study_accession}
echo ${vcf_dir}
echo ${tmp_dir}
echo ${logs_dir}
echo ${genomicsDB_dir}
echo ${metadata_file}
echo ${gvcf_file_suffix}
echo ${val_gvcf_file_suffix}
echo ${ref_file}
echo ${ref_index}
echo ${sample_list_file}
echo ${vcf_map_file}
echo ${failed_samples_file}
echo ${windows_file}
echo ${genotyped_vcf_suffix}
echo ${output_vcf_file}
printf "\n"

# Make temp dir (delete at end)
if [ ! -d ${tmp_dir} ]; then
    mkdir ${tmp_dir}
fi

# Make logs file if none
if [ ! -d ${logs_dir} ]; then
    mkdir ${logs_dir}
fi

# Make genomics DB directory if does not exist
if [ ! -d ${genomicsDB_dir} ]; then
    mkdir ${genomicsDB_dir}
fi

# Get samples from metadata based on study_accession number, store as file
grep ${study_accession} ${metadata_file} | cut -d, -f1 > ${sample_list_file}

echo "Samples:"
cat ${sample_list_file}
printf "\n"

# Make a list of the validated vcf file names
vcf_files=`cat ${sample_list_file} | sed "s|.*|${vcf_dir}&${val_gvcf_file_suffix}|"`

echo "vcf file names:"
echo ${vcf_files}
printf "\n"

# Commands
# see https://github.com/pathogenseq/fastq2matrix/blob/08480865bd2248ddb35a6d9e7321fff66fdd7a0c/scripts/merge_vcfs.py for use of gatk GenomicsDBImport 
bedtools_cmd="bedtools makewindows -n ${num_genome_chunks} -g ${ref_index}"
import_cmd="gatk GenomicsDBImport --genomicsdb-workspace-path ${genomicsDB_dir}{2}_genomics_db -L {1} --sample-name-map ${vcf_map_file} --reader-threads ${threads} --batch-size 500"
genotype_cmd="gatk GenotypeGVCFs -R ${ref_file} -V gendb://${genomicsDB_dir}{2}_genomics_db -O ${tmp_dir}{2}${genotyped_vcf_suffix}"

# Make windows for GenomicsDBImport and GenotypeGVCFs commands, store as file
${bedtools_cmd} | awk '{print $1":"$2+1"-"$3" "$1"_"$2+1"_"$3}' > ${windows_file}

# Get the name of the first windowed vcf file generated by GenotypeGVCFs to check if exists (i.e. "tmp/Chromosome_1_220576.val.gt.g.vcf.gz")
# first_gt_vcf_file=`cat x | head -1 | cut -f2 -d" " | sed "s|.*|${tmp_dir}&${genotyped_vcf_suffix}|"`
first_gt_vcf_file=`cat ${windows_file} | head -1 | cut -f2 -d" " | sed "s|.*|${tmp_dir}&${genotyped_vcf_suffix}|"`

# ------------------------------------------------------------------------------

# Samples map file and log file

# Generate tab-sep file of sample names and locations of GVCFs for GenomicsDBImport
if [ ! -s ${vcf_map_file} ] || [ ! -f ${vcf_map_file} ]; then
    touch ${vcf_map_file}
fi

# Generate failed samples log
if [ ! -s ${failed_samples_file} ] || [ ! -f ${failed_samples_file} ]; then
    touch ${failed_samples_file}
fi

# Check if validated gvcf files exist for each sample and add to map if exist
for samp in $(cat ${sample_list_file}); do
    if [ -f ${vcf_dir}${samp}${val_gvcf_file_suffix} ]; then
        # sed "s|.*|&\t${fastq_dir}&.gvcf.gz|" ${metadata_dir}${samples_list_file} >> ${metadata_dir}${study_accession}_sample_map.txt
        echo -e "${samp}\t${vcf_dir}${samp}${gvcf_file_suffix}" >> ${vcf_map_file}
    else
        printf "\n"
        echo -e "${samp}" >> ${failed_samples_file}
        printf "\n"
        echo "Failed sample!: ${samp}"
        printf "\n"
    fi
done

# Take uniq of sample map file otherwise gatk commands fail
sort ${vcf_map_file} | uniq > ${vcf_map_file}.uniq && mv ${vcf_map_file}.uniq ${vcf_map_file}

# Take uniq of failed samples file
sort ${failed_samples_file} | uniq > ${failed_samples_file}.uniq && mv ${failed_samples_file}.uniq ${failed_samples_file}

# ---------------------------------------------------------------------------------------------
# Run GenomicsDBImport - "Import single-sample GVCFs into GenomicsDB before joint genotyping"
# Run GenotypeGVCFs - Joint variant calling
# ---------------------------------------------------------------------------------------------

printf "\n"
# Run if genomics DB database is empty
if [ ! "$(ls -A ${genomicsDB_dir})" ]; then

    echo "Running GenomicsDBImport on ${vcf_map_file}"
    printf "\n"
    # Run in parallel
    cat ${windows_file} | parallel --bar -j ${threads} --col-sep " " ${import_cmd}
fi

printf "\n"

# Run if genomics DB database is NOT empty and the first genotyped (windowed) file does not exist
if [ "$(ls -A ${genomicsDB_dir})" ] && [ ! -f "${first_gt_vcf_file}" ]; then
    # Run GenotypeGVCFs to do joint variant calling
    # "This step runs quite quickly and can be rerun at any point when samples are added to the cohort" - see gatk link above
    echo "Running GenotypeGVCFs on ${genomicsDB_dir}"
    printf "\n"
    cat ${windows_file} | parallel -j ${threads} --col-sep " " ${genotype_cmd}
fi

# ------------------------------------------------------------
# Concatenate files (rbind) together with bcftools
# ------------------------------------------------------------

printf "\n"

if [ ! -f "${output_vcf_file}" ]; then
    echo "Concatenating VCFs, outputting multi-sample vcf ${output_vcf_file}"
    printf "\n"
    # ${bedtools_cmd} | awk '{print $1":"$2+1"-"$3" "$1"_"$2+1"_"$3}' | awk '{print $2".gn.genotyped.vcf.gz"}' > ${metadata_dir}genotyped_files.txt
    # bcftools concat -Oz -o ${final_concat_vcf_file} `${bedtools_cmd} | awk '{print $1":"$2+1"-"$3" "$1"_"$2+1"_"$3}' | awk '{print $2".gn.genotyped.vcf.gz"}'`
    bcftools concat -Oz -o ${output_vcf_file} `cat ${windows_file} | cut -f2 -d" " | sed "s|.*|${tmp_dir}&${genotyped_vcf_suffix}|"`
fi

# Index merged file
if [ ! -f "${output_vcf_file}.tbi" ]; then
    bcftools index -t ${output_vcf_file}
fi

# ------------------------------------------------------------------------------

# Clean up

rm -r ${tmp_dir}
rm -r ${genomicsDB_dir}
