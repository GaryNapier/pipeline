#! /bin/bash

# subset_fasta.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Description:
# Subset fasta file based on group from cluster analysis done in transmission_clusters.R
# transmission_clusters.R output file format:

# id	cluster
# ERR757145	1
# ERR757146	1
# ERR757147	1
# ERR757148	1
# ...
# ERR845304	2
# ERR845305	2
# ERR845306	2
# ERR845307	2
# ...etc

# i.e. one cluster group number assigned to each sample

# Split the fasta file based on the cluster number

# n.b. the fasta file will have already been 'filtered' in the transmission pipeline
# - i.e. samples that have not clustered (more than 50 SNPs away from any other sample) have been removed.

# Arguments to script:
# study accession, clusters file, fasta file

# Input:
# clusters file, fasta file

# Main steps:
# Filter the fasta file with seqtk and output as seperate fasta files

# Output:
# As many fasta files as there are unique cluster numbers in the clusters file
# Outputs to same folder as fasta intput file

# RUN
# Assume run from transmission/
# shell_scripts/subset_fasta.sh <study_accession> <metadata_file>                           <fasta_file>
# shell_scripts/subset_fasta.sh PRJEB7669         ~/metadata/tb_data_28_10_2021_clean.csv   fasta/PRJEB7669.filt.val.gt.g.snps.fa


# ------------------------------------------------------------------------------

# Setup

# Variables
study_accession=${1?Error: Enter study accession number}

# Input files
# clusters_file=${2?Error: Enter clusters file name}
metadata_file=${2?Error: Enter metadata file name}
fasta_file=${3?Error: Enter fasta file name}

# Directories
fasta_dir=`dirname ${fasta_file}`/

# Parameters
id_clust_tmp=id_clust_tmp

# Commands

# Print arguments
printf '\n'
echo 'Arguments to subset_fasta.sh:'
printf '\n'
echo "Study accession: ${study_accession}"
# echo "Clusters file: ${clusters_file}"
echo "Fasta file: ${fasta_file}"
printf '\n'

# Testing
# tail -n +2 ${clusters_file} | head -n 200 > test_clust_1
# tail -n +202 ${clusters_file} | sed 's/\t1/\t2/' > test_clust_2
# head -n 1 ${clusters_file}  | cat - test_clust_1 test_clust_2 > test_clusts

# Take unique clusters (1, 2, 3 etc) and store
# clusts=$(tail -n +2 ${clusters_file} | cut -f2 | sort -n | uniq)

cat ${metadata_file} | csvtk grep -f study_accession_word -p ${study_accession} | csvtk cut -f wgs_id,major_lineage > ${id_clust_tmp}

# Get list of clusters
clusts=$(cat ${id_clust_tmp} | cut -d',' -f2 | sort | uniq)

# Loop through the unique clusters and subset on cluster
for clust in ${clusts}; do
    # awk -v c=${clust} '$2 == c {print}' ${clusters_file} | cut -f1 | seqtk subseq ${fasta_file} - > ${fasta_dir}${study_accession}.clust_${clust}.fa
    cat ${id_clust_tmp} | \
    csvtk grep -f major_lineage -p ${clust} | \
    csvtk cut -f wgs_id | tail -n +2 | \
    seqtk subseq ${fasta_file} - > ${fasta_dir}${study_accession}.clust_${clust}.fa
    printf '\n'
    echo "---"
    echo "Outputting fasta file ${fasta_dir}${study_accession}.clust_${clust}.fa"
done

rm ${id_clust_tmp}

# Reset
# rm fasta/*.clust_*
