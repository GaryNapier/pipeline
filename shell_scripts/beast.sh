#! /bin/bash
# beast.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# Description:
# Run python script to concatenate sample IDs in fasta with dates from metadata.
# Run python script for Beauti to specify parameters and generate XML file for beast.
# Run Beast and generate newick file ready for itol.

# Arguments to script:
# study accession, fasta file location, metadata file name, name of date column in metadata file, xml script location, beast output location

# Input:
# fasta file of cluster from study accession, metadata file with dates for samples

# Main steps:
# Run fasta_add_annotations.py to concatenate the date from the metadata file with sample ID
# Script is from - https://github.com/pathogenseq/pathogenseq-scripts/tree/master/scripts
# Outputs new 'dated' fasta file

# Run python script for generating XML file for Beast with new fasta file from beast2-xml.py script from beast2-xml package
# See https://github.com/acorg/beast2-xml
# Outputs XML file for Beast

# Output:
# Dated fasta file
# XML file for beast
# Beast output files:
#   - <study_accession>.log
#   - <study_accession>.trees
#   - <study_accession>.xml.state

# RUN
# Assume run from transmission/
# beast.sh <study accession> <fasta_file_location> <metadata_file_name>                              <date_column>   <xml_script_location> <beast_output_location>
# beast.sh PRJEB7669        fasta/                 ../metadata/tb_data_collated_28_10_2020_clean.csv collection_date beast_xml             beast_results

# Notes:
# See also Beauti/Beast tutorial: https://taming-the-beast.org/tutorials/Prior-selection/Prior-selection.pdf

# ------------------------------------------------------------------------------

# Setup

# Variables
study_accession=${1?Error: input study accession number e.g. PRJEB7669}

# Directories
fasta_dir=${2?Error: input fasta file location}
xml_script_location=${5?Error input file location for XML output file}
beast_output_location=${6?Error input file location for Beast output files}

# Files
metadata_file=${3?Error: input metadata file name}
fasta_input_file=${fasta_dir}${study_accession}.filt.val.gt.g.snps.fa.filt
dated_fasta_file=${fasta_input_file}.dated.fa
xml_file=${xml_script_location}${study_accession}.xml

# Parameters
date_column=${4?Error input date column name from metadata}

# Commands

# Print arguments
printf '\n'
echo 'Arguments:'
printf '\n'
echo "Study accession: ${study_accession}"
echo "Fasta file:  ${fasta_input_file}"
echo "Dated fasta file: ${dated_fasta_file}"
echo "Metadata file: ${metadata_file}"
echo "Date column name: ${date_column}"
echo "XML file: ${xml_file}"
echo "Beast output files location: ${beast_output_location}"


printf '\n'

# ------------------------------------------------------------------------------

# Run fasta_add_annotations.py to concatenate the metadata collection date with the sample ID in the fasta file

# e.g. >ERR1234 -> >ERR1234_01_01_10
# - see https://github.com/pathogenseq/pathogenseq-scripts/tree/master/scripts
# fasta_add_annotations.py --fasta fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt --csv ../metadata/tb_data_collated_28_10_2020_clean.csv --out fasta/PRJEB7669.filt.val.gt.g.snps.fa.filt.dated.fa --annotations collection_date
fasta_add_annotations.py --fasta ${fasta_input_file} --csv ${metadata_file} --out ${dated_fasta_file} --annotations ${date_column}


# ------------------------------------------------------------------------------

# Generate XML file using beast2-xml.py (replaces Beauti GUI) - https://github.com/acorg/beast2-xml

# beast2-xml.py [-h] [--clockModel MODEL | --templateFile FILENAME]
#                      [--chainLength LENGTH] [--age ID=N [ID=N ...]]
#                      [--defaultAge N] [--dateUnit UNIT]
#                      [--dateDirection DIRECTION]
#                      [--logFileBasename BASE-FILENAME] [--traceLogEvery N]
#                      [--treeLogEvery N] [--screenLogEvery N] [--mimicBEAUti]
#                      [--sequenceIdDateRegex REGEX]
#                      [--sequenceIdDateRegexMayNotMatch] [--fastaFile FILENAME]
#                      [--readClass CLASSNAME] [--fasta | --fastq | --fasta-ss]
#
#
# beast2-xml.py --chainLength 100000000 --fastaFile ${dated_fasta_file}


# ------------------------------------------------------------------------------

# Run beast

beast ${xml_file}

# ------------------------------------------------------------------------------

# Clean up

mv ${study_accession}.log ${study_accession}.trees ${study_accession}.xml.state ${beast_results_dir}

# ------------------------------------------------------------------------------

# NOTES

# Beast parameters

# See Transmission analysis of a large tuberculosis outbreak in London: a mathematical modelling study using genomic data Open Access - Yuanwei Xu - 11 November 2020 https://doi.org/10.1099/mgen.0.000450

# The phylogenetic tree-building software beast2 (version 2.6.1) [12] was used to build timed phylogenetic trees.

# A preliminary check using TempEst [13] showed positive correlation between genetic divergence and sampling time and a moderate level of temporal signal (TempEst R^2 =0.21).

# Because of moderate temporal signal in the SNP data, we adopted a strict molecular clock, supplying the tip dates, and we used a fixed rate parameter of 1.0×10 − 7 per site per year,
# corresponding to 0.44 substitutions per genome per year [14].

# Tip Dates > 'as dates with format dd/M/yyyy'

# Clock Model > Strict Clock > Clock.rate = 1.0E-7

# We used a coalescent constant population model with a log-normal [0, 200] prior for the population size.

# Priors > Tree.t = Coalescent Constant Population
# Priors > popSize > Log Normal > initial = [Lower = 0], [Upper = 200]

# Because the K3Pu model of nucleotide substitution was not available in beast2, we used the generalized time reversible (GTR) substitution model [17],
# which had the next lowest Bayesian information criterion (BIC) score (Delta6910.964) on the basis of model testing using iq-tree [18].

# Site Model > Change JC69 to GTR

# The GTR model with prior rates having a gamma distribution with rates in [0, inf] and prior frequencies (estimated) in [0, 1] were applied, along with 0 proportion of invariant sites.

# Change nothing? - yes

# [Why is Rate CT estimated unchecked?]

#  We used the beast2 correction for ascertainment bias, specifying the number of invariant A, C, G and T sites as 758511 1449901 1444524 758336.
#  Note that this must be manually added to the xml and may not appear when the xml is loaded into the BEAUti2 (version 2.6.1) software.

# [???]

# Line in the xml file:
# <constantSiteWeights id="IntegerParameter.0" spec="parameter.IntegerParameter" dimension="4" lower="0" upper="0">758511 1449901 1444524 758336</constantSiteWeights>

#  We ran the Markov chain Monte Carlo (MCMC) method for 100000000 iterations, sampling every 10000th iteration.

#  We verified chain convergence (by confirming multiple independent chains converged to the same posterior values) as well as good mixing and an effective sample size (ESS)
#  of greater than 200 for all parameters using Tracer (version 1.7.1) [19].
#
#  A maximum clade credibility (MCC) tree was created using TreeAnnotator (version 2.6.0), with 10% of the chain discarded as burn-in, resulting in a posterior collection of 9000 trees.
#
#  Instead of trying to obtain a single optimal timed phylogenetic tree from this posterior set, we sampled a collection of 50 of them at random.
#  This ensures that we capture as much diversity as possible from the beast posterior, to achieve robust uncertainty quantification in our subsequent analysis.
