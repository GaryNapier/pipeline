#! /bin/bash

# variant_filtering.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Do variant filtering (TB specific)
# Takes multi-sample vcf (output from variant_calling_and_concat_gvcfs.sh) and removes variants according to criteria (see Steps below)
# NOTE: Needs setGT.py script installed as executable (see pathogenseq link below)

# Adapted from:
# https://github.com/pathogenseq/fastq2matrix/tree/master/scripts
# filter_tb_vcf.py
# And see https://github.com/pathogenseq/fastq2matrix/blob/master/fastq2matrix/vcf.py for vcf_class

# Arguments to script:
# multi-sample vcf file name, output (filtered multi-sample vcf) file name,  bed file name excluding regions, reference fasta file name

# Input:
# multi-sample vcf file, bed file, reference fasta

# Main steps:
# Remove indels > Exclude regions in bed file > Clean vcf formatting and convert mixed calls according to threshold (setGT) >
# filter out SNPs with no samples and trim alt alleles with no samples > change heterozygous calls to unknown (./.) >
# filter out SNPs where 10% of samples are unknown (./.) > filter out SNPs with no samples again >
# recalculate 'metadata' columns e.g. Allele Count (AC) > remove monomorphic > trim the Ref/Alt cols if deletions have dropped out

# Output:
# Filtered multi-sample vcf

# RUN
# Assume run from ~/transmission :
# shell_scripts/variant_filtering.sh <multi-sample vcf file name>       <output file name>                      <bed file name>               <ref fasta>
# shell_scripts/variant_filtering.sh ~/vcf/PRJEB7669.val.gt.g.vcf.gz    ~/vcf/PRJEB7669.filt.val.gt.g.vcf.gz    ~/refgenome/excluded_loci.bed ~/refgenome/MTB-h37rv_asm19595v2-eg18.fa

# ------------------------------------------------------------------------------

# Setup

# Variables

# Directories

# Files
input_concat_vcf_file=${1}
output_vcf_file=${2}
exclude_regions_bed_file=${3}
ref_fasta=${4}

# Parameters

# Commands

printf "\n"
echo "Parameters:"
printf "\n"
echo ${input_concat_vcf_file}
echo ${output_vcf_file}
echo ${exclude_regions_bed_file}
echo ${ref_fasta}
printf "\n"

# ------------------------------------------------------------------------------

# -------------
# Do filtering
# -------------

printf "\n"
echo "Filtering vcf"
printf "\n"

# Exclude indels
bcftools view -V indels ${input_concat_vcf_file} | \
# Exclude regions (repetitive regions. Will include some but not all PPE)
bcftools view -T ^${exclude_regions_bed_file} | \
# 'Clean' VCF file - convert all '|' to '/'. Convert e.g. 1/0:10,90 to 0/0 (score of 90 for 1 overrides score of 10 for 1)
setGT.py | \
# -c 1: Make sure all SNPs have at least one instance (sample) (for each row (SNP), need at least one sample).
# Some rows (SNPs) may now be all 0/0 after previous step
# -a: Trim off alternate alleles (in the Alt/Ref cols) which have gone to 0 in setGT
# e.g. G    A,T - the T is now 0/0, so goes to G    A
# See notes below
# Output as uncompressed bcf (for speed. Otherwise next command has to unzip/uncompress)
bcftools view -c 1 -a  | \
# Change heterozygous calls to "." (unknown)
# Command first says to exclude (-e), but second command (-S) says to convert
bcftools filter -e 'GT="het"' -S . | \
# Include ratio of AN to AC (F_PASS) - has to be at least 90% non-missing (./.)
# i.e. if 10% samples are ./. then remove
bcftools view -i 'F_PASS(GT!="mis")>0.9' | \
# Re-apply 'at least one sample' rule
bcftools view -c 1 | \
# Re-calculate 'metadata' cols. e.g. Allele Counts (AC)
bcftools +fill-tags | \
# Exclude monomorphic - all samples have the SNP or all do not (AF = allele freq)
bcftools view -e 'AF==1 || AF==0' | \
# Trim the ref/alt if any deletions have dropped out:
# See notes below
bcftools norm -f ${ref_fasta} -Oz -o ${output_vcf_file}

# Index filtered file
if [ ! -f "${output_vcf_file}.tbi" ]; then
    bcftools index -t ${output_vcf_file}
fi


# NOTES:

# Explanation of x/x SNP data:

# Ref -  .....TG.....
# read1  .....T......
# read2  .....T......
# read3  .....T......
# read4  .....T......
# read5  .....T-.....
# read6  .....T......
# read7  .....T......
# read8  .....T......
# read9  .....T......
# read10 .....T......

# Represented as:

# Ref     Alt   Sample1
# TG      T     0/1:90,10
# - So 90% of reads say it's the reference, and 1 read says there's a deletion where the G is.

# In this case, setGT.py would re-set the SNP call to:
# Ref     Alt   Sample1
# TG      T     0/0
# because it's below the threshold specified in setGT

# Note that the Ref and Alt cols have not changed. These need to be changed with further commands (see below)


# Example with deletion and SNP in two diff samples
# Full stacks of reads are not shown (as in the first example) - see scores in the vcf representation for the proportions of reads agreeing with ref/alt

# Step 1. setGT

# From:

# Ref - .....TG.....
# S1  - .....T-.....
# S2  - .....AG.....

# Represented in vcf as:
#
# Ref     Alt     S1            S2
# TG      T, AG   0/1:90,10,0   2/2:0,0,10

# In Sample 1, 10% of reads are saying the genotype is "T-", that is, T followed by a deletion (as we saw above)
# In sample 2, 10% of reads are saying there is a SNP at the T position (i.e. T -> A), but the G remains in place
# It's coded as 2/2 because it's the 3rd candidate genotype

# In setGT, this converts to:

# Ref     Alt     S1     S2
# TG      T, AG   0/0    2/2

# Note the Ref and Alt haven't changed

# Step 2. bcftools view -a
# Represented in vcf as:
#
# Ref     Alt       S1      S2
# TG      T,AG      0/0     2/2

# bcftools view takes away the alternative allele in the Alt col because it is no longer there in any sample (has gone to 0/0):
# Ref     Alt     S1      S2
# TG      AG      0/0     1/1

# Note also the 2/2 is now coded as 1/1

# Step 3. bcftools norm #
# Ref - .....TG.....
#
# S2  - .....AG.....
#
# Note Sample 1 has gone

# Now bcftools norm just removes the redundant G from the Ref and Alt cols.
# Ref     Alt     S1      S2
# T       A       0/0     1/1
