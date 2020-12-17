#! /bin/bash

# concat_filt_vcfs.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Concatenate and filter vcfs of samples specific to project ready for distance matrix calculation
# vcfs must be pre-processed up to gvcf stage

# Arguments to script:
# study accession number (from metadata file column), metadata file, gvcf directory location, gvcf file suffix from ready gvcfs, reference fasta

# Input:
# study accession number, metadata file, gvcf directory location, gvcf file suffix from ready gvcfs, reference fasta, gvcf files

# Main steps:
# Check if validated gvcf files exist for each sample and add to map if exist > GenomicsDBImport (import gvcfs to database) >
# GenotypeGVCFs (joint variant calling) > bcftool concat (concatenate (rbind) the split genome into one vcf) > index merged file

# Output:
# Concatenated and filtered vcf file

# RUN
# concat_filt_vcfs.sh <study_accession> <metadata_file> <vcf_dir> <gvcf_file_suffix> <ref_file>
# concat_filt_vcfs.sh PRJEB7669 ../metadata/tb_data_collated_28_10_2020_clean.csv ~/vcf/ .g.vcf.gz ~/refgenome/MTB-h37rv_asm19595v2-eg18.fa

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
val_gvcf_file_suffix=${val_gvcf_file_suffix}.validated
ref_file=${5}
ref_index=${ref_file}.fai
vcf_map_file=${tmp_dir}vcf_map_file
failed_samples_file=${logs_dir}failed_samples
genotyped_vcf_suffix=val.gt.${gvcf_file_suffix}
output_vcf_file=${vcf_dir}${project_code}.${genotyped_vcf_suffix}

# Parameters
threads=4
num_genome_chunks=20

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

# Get samples from metadata based on study_accession number
sample_list=`grep ${study_accession} ${metadata_file} | cut -d, -f1`

# Make a list of the validated vcf files and put into temp file
vcf_files=`${samples_list} | sed "s|.*|${vcf_dir}&${val_gvcf_file_suffix}|"`

# Commands
bedtools_cmd="bedtools makewindows -n ${num_genome_chunks} -g ${ref_index}"
import_cmd="gatk GenomicsDBImport --genomicsdb-workspace-path ${genomicsDB_dir}{2}_genomics_db -L {1} --sample-name-map ${vcf_map_file} --reader-threads ${threads} --batch-size 500"
genotype_cmd="gatk GenotypeGVCFs -R ${ref_fasta} -V gendb://${genomicsDB_dir}{2}_genomics_db -O ${vcf_dir}{2}${genotyped_vcf_suffix}"


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
for samp in ${sample_list}; do
    if [ -f ${vcf_dir}${samp}${val_gvcf_file_suffix} ]; then
        # sed "s|.*|&\t${fastq_dir}&.gvcf.gz|" ${metadata_dir}${samples_list_file} >> ${metadata_dir}${project_code}_sample_map.txt
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
uniq ${vcf_map_file} > ${vcf_map_file}.uniq && mv ${vcf_map_file}.uniq ${vcf_map_file}

# Take uniq of failed samples file
uniq ${failed_samples_file} > ${failed_samples_file}.uniq && mv ${failed_samples_file}.uniq ${failed_samples_file}


# ---------------------------------------------------------------------------------------------
# Run GenomicsDBImport - "Import single-sample GVCFs into GenomicsDB before joint genotyping"
# Run GenotypeGVCFs - Joint variant calling
# ---------------------------------------------------------------------------------------------

printf "\n"
# Run if database is empty
if [ ! "$(ls -A ${genomicsDB_dir})" ]; then

    echo "Running GenomicsDBImport on ${vcf_map_file}"
    # Run in parallel
    ${bedtools_cmd} | awk '{print $1":"$2+1"-"$3" "$1"_"$2+1"_"$3}' | parallel --bar -j ${threads} --col-sep " " ${import_cmd}

    # Run GenotypeGVCFs to do joint variant calling
    # "This step runs quite quickly and can be rerun at any point when samples are added to the cohort" - see gatk link above
    printf "\n"
    echo "Running GenotypeGVCFs on ${genomicsDB_dir}"
    ${bedtools_cmd} | awk '{print $1":"$2+1"-"$3" "$1"_"$2+1"_"$3}' | parallel -j ${threads} --col-sep " " ${genotype_cmd}

fi

# ------------------------------------------------------------
# Concatenate files together with bcftools
# ------------------------------------------------------------

# bcftools concat -Oz -o genotyped.vcf.gz `%(window_cmd)s | awk '{print \"%(prefix)s.\"$2\".genotyped.vcf.gz\"}'
if [ ! -f "${output_vcf_file}" ]; then
    echo "Concatenating VCFs"
    printf "\n"
    # ${bedtools_cmd} | awk '{print $1":"$2+1"-"$3" "$1"_"$2+1"_"$3}' | awk '{print $2".gn.genotyped.vcf.gz"}' > ${metadata_dir}genotyped_files.txt
    # bcftools concat -Oz -o ${final_concat_vcf_file} `${bedtools_cmd} | awk '{print $1":"$2+1"-"$3" "$1"_"$2+1"_"$3}' | awk '{print $2".gn.genotyped.vcf.gz"}'`
    bcftools concat -Oz -o ${output_vcf_file} `${bedtools_cmd} | awk '{print $1":"$2+1"-"$3" "$1"_"$2+1"_"$3}' | cut -f2 -d" " | sed "s|.*|${vcf_dir}&${genotyped_vcf_suffix}|"`
fi

# Index merged file
if [ ! -f "${output_vcf_file}.tbi" ]; then
    bcftools index -t ${output_vcf_file}
fi

# ------------------------------------------------------------------------------

# Clean up

rm -r tmp
