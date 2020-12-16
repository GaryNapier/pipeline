#! /bin/bash

# concat_filt_vcfs.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Concatenate and filter vcfs of samples specific to project ready for distance matrix calculation
# vcfs must be pre-processed up to gvcf stage

# Arguments to script:
# Project code, vcf directory location, list of samples, input gvcf suffix

# Input:
# Project code, vcf directory location, list of samples, input gvcf suffix, gvcf files

# Steps:
# Concat vcfs corresponding to samples > filter vcfs > output vcf file

# Output:
# Concatenated and filtered vcf file

# RUN
# concat_filt_vcfs.sh <project_code> <vcf_dir> <sample_list> <input_vcf_file_suffix>
# concat_filt_vcfs.sh PRJEB7669 ~/vcf/ ../metadata/PRJEB7669_ena_samps.txt .g.vcf.gz.validated

# ------------------------------------------------------------------------------

# Variables
nl="\n"
project_code=${1}

# Directories
vcf_dir=${2}

# Files
sample_list=$(cat ${3})
input_vcf_file_suffix=${4}
output_vcf_file=${vcf_dir}${project_code}.val.concat.filt.no_indels.vcf.gz
files=$(cat ${3} | sed "s|.*|${vcf_dir}&${input_vcf_file_suffix}|")

# ------------------------------------------------------------------------------

# Make temp dir (delete at end)
if [ ! -d tmp ]; then
    mkdir tmp
fi

# ------------------------------------------------------------------------------

# Check vcf files exist and are not empty and store as a temp file if so

if [ ! -f tmp/${project_code}_samples.txt ]; then
    touch tmp/${project_code}_samples.txt
fi

for file in ${files}; do
    if [ -f ${file} ] && [ ! -s ${file} ]; then
        echo ${file} >> tmp/${project_code}_samples.txt
    fi
done


# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Run ValidateVariants - "validate the adherence of a file to VCF format" - https://gatk.broadinstitute.org/hc/en-us/articles/360042914291-ValidateVariants.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------

for samp in $(cat tmp/${project_code}_samples.txt); do
    if [ ! -f "${vcf_dir}${samp}.gvcf.gz.validated" ]; then
        printf "\n"
        echo "Running ValidateVariants on ${samp}"
        gatk ValidateVariants -R ${ref_dir}${ref_fasta} -V ${vcf_dir}${samp}.gvcf.gz -gvcf && touch ${vcf_dir}${samp}${validated_vcf_suffix}
        printf "\n"
    fi
done

printf "\n"

# ------------------------------------------------------------------------------
#
# # Exclude indels
# bcftools view -V indels ${output_vcf_file} | \
# # Exclude regions (repetitive regions. Will include some but not all PPE)
# bcftools view -T ^${exclude_regions_bed_file} | \
# # 'Clean' VCF file - convert all '|' to '/'. Convert e.g. 1/0:10,90 to 0/0 (score of 90 for 1 overrides score of 10 for 1)
# setGT.py | \
# # -c 1: Make sure all SNPs have at least one instance (sample) (for each row (SNP), need at least one sample).
# # Some rows (SNPs) may now be all 0/0 after previous step
# # -a: Trim off alternate alleles (in the Alt/Ref cols) which have gone to 0 in setGT
# # e.g. G    A,T - the T is now 0/0, so goes to G    A
# # See notes below
# # Output as uncompressed bcf (for speed. Otherwise next command has to unzip/uncompress)
# bcftools view -c 1 -a  | \
# # Change heterozygous calls to "." (unknown)
# # Command first says to exclude (-e), but second command (-S) says to convert
# bcftools filter -e 'GT="het"' -S . | \
# # Include ratio of AN to AC (F_PASS) - has to be at least 90% non-missing (./.)
# # i.e. if 10% samples are ./. then remove
# bcftools view -i 'F_PASS(GT!="mis")>0.9' | \
# # Re-apply 'at least one sample' rule
# bcftools view -c 1 | \
# # Re-calculate 'metadata' cols. e.g. Allele Counts (AC)
# bcftools +fill-tags | \
# # Exclude monomorphic - all samples have the SNP or all do not (AF = allele freq)
# bcftools view -e 'AF==1 || AF==0' | \
# # Trim the ref/alt if any deletions have dropped out:
# # See notes below
# bcftools norm -f ${ref_fasta} -Ob -o ${output_vcf_file}
#
# # ------------------------------------------------------------------------------
#
# # Clean up
# rm -r tmp
