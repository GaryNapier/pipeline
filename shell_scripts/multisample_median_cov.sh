#! /bin/bash

# multisample_median_cov.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------



# Arguments to script:

# - metadata for list of samples
# - location of bam/cram files
# - output location

# Input:

study_accession=${1}
metadata_file=${2}
bam_dir=${3}

# Output
output_dir=${4}
median_cov_file=${output_dir}${study_accession}.median_cov.txt


# Main steps:

# - pull sample IDs from metadata file
# - put sample IDs through parallel and run bedtools genomecov on each position (-d) - pipe to calculate median coverage per genome pos
# - store in file
# - calculate median on this median

# Output:

# - single number: median coverage for given samples

# RUN
# multisample_median_cov.sh <study_accession> <metadata_file> <bam_dir> <output_dir>


parallel_cmd="bedtools genomecov -d -ibam ${bam_dir}{}.bqsr.cram | datamash median 3"
cat ${metadata_file} | csvtk grep -f study_accession_word -p ${study_accession} | csvtk cut -f wgs_id | tail -n +2 | parallel -j 4 --bar ${parallel_cmd} > ${median_cov_file}
printf "multisample_median_cov.sh RESULTS \n"
printf "median coverage: \n"
cat ${median_cov_file} | datamash median 1
printf "min cov: \n"
cat ${median_cov_file} | datamash min 1
printf "max cov: \n"
cat ${median_cov_file} | datamash max 1
