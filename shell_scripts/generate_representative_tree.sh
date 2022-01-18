#!/bin/bash

# Get representative numbers of samples per lineage (10, say) and use these samples to generate tree

# Input arguments

# Enter number of samples per lineage
n_samps=${1?Error: enter number of samples per lineage}

# Meta file with all lineages per sample - both training and test - 35972 samples
# id,lineage,type
# ERR553088,1.2,training
meta_file=${2?Error: enter lineages csv file} # ~/metadata/meta_file_clean.csv

# Location of vcfs - input to merge_vcfs
vcf_dir=${3?Error: enter location of vcf files}

# Reference
ref_fasta=~/refgenome/MTB-h37rv_asm19595v2-eg18.fa

# Make temp directory for data
if [ ! -d ~/tmp/ ]; then
    mkdir ~/tmp/
fi

# Remove created files to reset so not appending
rm -f ~/tmp/unique_lins.txt ~/metadata/representative_samples_${n_samps}_per_lineage.csv representative_samples_${n_samps}_per_lineage

# Get all unique lineages, cleaning out mixed samples (";") blank lines and first line
# unique_lins=`cut -d, -f2 ${meta_file} | sort | uniq | grep -v ";"  | grep -v "lineage" | grep -v -e '^$'`  # Gets 86 lineages. Missing: 1.3, 2, 8, 9
cut -d, -f2 ${meta_file} | sort | uniq | grep -v ";"  | grep -v "lineage" | grep -v -e '^$' > ~/tmp/unique_lins.txt # Gets 86 lineages. Missing: 1.3, 2, 8, 9

# # Append 8 and 9 to unique_lins.txt. Don't need 1.3 and 2 as they are comprised of lower-level samples
# echo 8 >> ~/tmp/unique_lins.txt
# echo 9 >> ~/tmp/unique_lins.txt


echo "Randomly sampling ${n_samps} samples per lineage"
# Get X number of samples from each lineage from meta_file_clean.csv

# Set up loop
unique_lins=$(cat ~/tmp/unique_lins.txt)

for lineage in ${unique_lins}
do
    # echo ${lineage}
    echo ${lineage} | cgrep ${meta_file} 2 --sep "," | shuf -n ${n_samps} >> ~/metadata/representative_samples_${n_samps}_per_lineage.csv
done

# Cut the first collumn to get the samples only to pass to merge_vcfs.sh
cut -d, -f1 ~/metadata/representative_samples_${n_samps}_per_lineage.csv > representative_samples_${n_samps}_per_lineage

printf "\n"

echo "Merging VCFs for selected samples"
# Merge vcfs for those samples selected
~/S7_shell_scripts/merge_vcfs.sh representative_samples_${n_samps}_per_lineage ${vcf_dir}

printf "\n"

echo "Converting genotyped and filtered multi VCF to multi-fasta file"
# Make fasta file of genotyped and filtered multi-vcf
vcf2fasta.py --vcf ${vcf_dir}/representative_samples_${n_samps}_per_lineage.*.genotyped.filtered.vcf.gz --ref ${ref_fasta} --threads 20 --snps

# Output:
# ${vcf_dir}/representative_samples_${n_samps}_per_lineage.[date].genotyped.filtered.snps.fa
# ${vcf_dir}/representative_samples_${n_samps}_per_lineage.[date].genotyped.filtered.vcf.gz.csi

mv ${vcf_dir}/representative_samples_${n_samps}_per_lineage.*.genotyped.filtered.snps.fa ~/fasta

rm ${vcf_dir}/representative_samples_${n_samps}_per_lineage.*.genotyped.filtered.vcf.gz.csi

printf "\n"

echo "Running tree with ~/fasta/representative_samples_${n_samps}_per_lineage.*.genotyped.filtered.snps.fa"
# Run tree
iqtree -s ~/fasta/representative_samples_${n_samps}_per_lineage.*.genotyped.filtered.snps.fa -m GTR+G -nt AUTO
# nb.: Had to remove +ASC because: ERROR: Invalid use of +ASC because of 129 invariant sites in the alignment

mv ~/fasta/*.iqtree ~/fasta/*.treefile ~/fasta/*.log ~/newick

# Remove files
rm representative_samples_${n_samps}_per_lineage ${vcf_dir}/representative_samples_${n_samps}_per_lineage.*.genotyped.filtered.vcf.gz

# FOLDERS
# S7_shell_scripts
# ~/vcf/
# ~/tmp/
# ~/fasta
# ~/metadata
# ~/newick

# TESTING
# ./generate_representative_tree.sh 10 ~/metadata/TWO_SAMPLES.csv ~/vcf

# RUN
# ./generate_representative_tree.sh 10 ~/metadata/meta_file_clean_samples_with_vcfs.csv ~/vcf
