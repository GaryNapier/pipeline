#!/bin/bash

# Merge individual VCFs into one massive VCF

# See https://lshtm.gitbook.io/host-pathogen-tb/pathogen/untitled#merging-vcf-files

# Arguments

samples_file=${1?Error: enter input file} # Input samples file as simple list of samples <sample1> \n <sample2> \n etc
vcf_dir=${2?Error: enter VCF directory}
ref_fasta=~/refgenome/MTB-h37rv_asm19595v2-eg18.fa
prefix=${vcf_dir}/${samples_file}

# echo "Getting samples only from samples/lineage file"
# # Make samples only file
# cut -d, -f1 ${samples_file} > ${samples_file}.txt

printf "\n"

echo "Merging sample VCFs into one VCF and genotyping from ${vcf_dir}"

printf "\n"

# Import the vcf files into a genomics database (gatk)
merge_vcfs.py all --sample-file ${samples_file} --ref ${ref_fasta} --prefix ${prefix} --vcf-dir ${vcf_dir} --threads 20

# Outputs:
# ${prefix}.[date].genotyped.vcf.gz
# ${prefix}.dbconf.json
# ${prefix}.failed_samples.log
# ${prefix}.map
# ${prefix}*_genomics_db

# # Run the joint genotyping - don't need to run if above argument is "all" ("python fastq2matrix/scripts/merge_vcfs.py all")
# python fastq2matrix/scripts/merge_vcfs.py genotype --ref ${ref_fasta} --prefix ${prefix}

printf "\n"

echo "Filtering genotyped multi-VCF"
# Filtering on the file to remove false positive calls. Do the filtering with the following command which will use recommended parameters for Mtb
filter_tb_vcf.py --ref ${ref_fasta} --vcf ${prefix}.*.genotyped.vcf.gz --threads 20

# Outputs:
# ${prefix}.[date].genotyped.filtered.vcf.gz

# Clean up
rm -r ${prefix}*_genomics_db ${prefix}*.genotyped.vcf.gz ${prefix}*.genotyped.vcf.gz.csi ${prefix}*.dbconf.json ${prefix}*.map

# Leaves:
# ${prefix}.failed_samples.log
# ${prefix}.[date].genotyped.filtered.vcf.gz

printf "\n"
if [ -s ${prefix}*.failed_samples.log ]; then
    echo "FILE ${prefix}*.failed_samples.log FAILED SAMPLES:"
    cat ${prefix}*.failed_samples.log
    printf "\n"
else
    echo "FILE ${prefix}*.failed_samples.log IS EMPTY - REMOVING"
    rm ${prefix}*.failed_samples.log
    printf "\n"
fi

# If empty, leaves ${prefix}.genotyped.filtered.vcf.gz

printf "\n"
