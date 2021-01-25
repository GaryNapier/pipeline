#!/usr/bin/env Rscript

# beast_xml.R


# Description
# 

# Arguments to script:

# Input:

# Main steps:

# Output:

# RUN:
  # Rscript r_scripts/beast_xml.R -s <study_acc>  -t <template_file_name>      -f <dated_fasta_file_name> -o <output_file_location>
# Rscript r_scripts/beast_xml.R -s PRJEB7669    -t beast_xml/tb_template.xml -f fasta/                    -o beast_xml/

# Setup ----

library(optparse)

heaD <- function(x,...){
  head(x, ...)
}

len_str <- function(string){
  length(unlist(strsplit(string, split = "")))
}

hs <- function(x, ...){
  print(head(x, ...))
  print("---")
  str(x, ...)
}


# Arguments ----

option_list = list(
  make_option(c("-s", "--STUDY_ACCESSION"), type="character", default=NULL, 
              help="input study accession number", metavar="character"), 
  
  make_option(c("-t", "--template_file_name"), type="character", default=NULL, 
              help="input template xml file", metavar="character"), 
  
  make_option(c("-f", "--dated_fasta_file_name"), type="character", default=NULL, 
              help="input template xml file", metavar="character"), 
  
  make_option(c("-o", "--output_file_location"), type="character", default="./", 
              help="file location to save output file to", metavar="character"), 
  
  make_option(c("-M", "--MUTATIONRATE"), type="character", default="1.0", 
              help="input mutation rate, default = 1.0", metavar="character"), 
  
  make_option(c("-C", "--CLOCKRATE"), type="character", default="1.0E-7", 
              help="input clockrate, default = 1.0E-7", metavar="character"), 
  
  make_option(c("-P", "--POPSIZE"), type="character", default="0.3", 
              help="input popsize, default = 0.3", metavar="character"), 
  
  make_option(c("-L", "--POPSIZELOWER"), type="character", default="0.0", 
              help="input popsize, lower value, default = 0.0", metavar="character"),
  
  make_option(c("-U", "--POPSIZEUPPER"), type="character", default="200.0", 
                help="input popsize, upper value, default = 200.0", metavar="character"),
  
  make_option(c("-F", "--FREQPARAMETER"), type="character", default="0.25", 
              help="input freqparameter, default = 0.25", metavar="character")
  
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)
print("---")
print(str(opt))

# Save parameter names to loop through and do gsub (except DATA)
# Note order of params - POPSIZE has to come last because is a subset of the other two. 
params <- c("STUDY_ACCESSION", "MUTATIONRATE", "CLOCKRATE", "POPSIZELOWER", "POPSIZEUPPER", "POPSIZE", "FREQPARAMETER")

# Set up DATA sub-template
DATA <- "<sequence id=\"seq_ID\" spec=\"Sequence\" taxon=\"ID\" totalcount=\"4\" value=\"SEQ\"/>"


# -------------------
# # TESTING
# STUDY_ACCESSION <- "PRJEB7669"
# template_file_name <- "beast_xml/tb_template.xml"
# dated_fasta_file_location <- "fasta/"
# output_file_location <- "beast_xml/"
# MUTATIONRATE <- "1.0"
# CLOCKRATE <- "1.0E-7"
# POPSIZE <- "0.3"
# POPSIZELOWER <- "0.0"
# POPSIZEUPPER <- "200.0"
# FREQPARAMETER <- "0.25"
# 
# opt <- c(STUDY_ACCESSION, template_file_name, dated_fasta_file_location, output_file_location, 
#               MUTATIONRATE, CLOCKRATE, POPSIZE, POPSIZELOWER, POPSIZEUPPER, FREQPARAMETER)
# 
# x <- list()
# for(i in opt){
#   x[[i]] <- i
# }
# 
# names(x) <- c("STUDY_ACCESSION", "template_file_name", "dated_fasta_file_location", "output_file_location", 
#               "MUTATIONRATE", "CLOCKRATE", "POPSIZE", "POPSIZELOWER", "POPSIZEUPPER", "FREQPARAMETER")
# 
# opt <- x
# 
# print(opt)

# -------------------
# TESTING



# Files ----

fasta_file <- paste0(opt$dated_fasta_file_location, opt$STUDY_ACCESSION, ".filt.val.gt.g.snps.fa.filt.dated.fa")
output_file <- paste0(opt$output_file_location, opt$STUDY_ACCESSION, ".xml")

print("FILES:")
print(c("dated fasta file name: ", fasta_file))
print(c("output file name: ", output_file))


# Read in data ----

template <- readChar(template_file_name, nchars = 1e6)
fasta <- readLines(fasta_file)

# Substitutions ----

# Do substitutions of parameters except DATA
for(param in params){
  template <- gsub(param, opt[[param]], template)
}


# DATA substitution using fasta file

# Get samples from fasta
samples <- grep("^>", fasta, value = T, perl = T)
# Strip ">" from samples 
samples <- gsub(">", "", samples)

# Get seqs from fasta
seqs <- fasta[!grepl("^>", fasta, perl = T)]

# Loop through and gsub in DATA sub-template
DATA_list <- vector()
for(i in seq(samples)){
  DATA_list[i] <- gsub("ID", samples[i], DATA)
  DATA_list[i] <- gsub("SEQ", seqs[i], DATA_list[i])
}

# Put together as one string separated by new lines
DATA_list <- paste(DATA_list, collapse = "\n")

# Substitute DATA in main template
template <- gsub("DATA", DATA_list, template)

# Write to file
write.table(template, file = output_file,
            row.names=F, col.names=F, quote = F)

# write.table(template, file = "beast_xml/TEST.xml",
#             row.names=F, col.names=F, quote = F)



# NOTES ----

# Beast parameters

# See Transmission analysis of a large tuberculosis outbreak in London: a mathematical modelling study using genomic data 
# Open Access - Yuanwei Xu - 11 November 2020 https://doi.org/10.1099/mgen.0.000450

# The phylogenetic tree-building software beast2 (version 2.6.1) [12] was used to build timed phylogenetic trees.

# A preliminary check using TempEst [13] showed positive correlation between genetic divergence and sampling time and a 
# moderate level of temporal signal (TempEst R^2=0.21).

# Because of moderate temporal signal in the SNP data, we adopted a strict molecular clock, supplying the tip dates, 
# and we used a fixed rate parameter of 1.0x10^7 per site per year, corresponding to 0.44 substitutions per genome per year.

# Tip Dates > 'as dates with format dd/M/yyyy'

# Clock Model > Strict Clock > Clock.rate = 1.0E-7

# We used a coalescent constant population model with a log-normal [0, 200] prior for the population size.

# Priors > Tree.t = Coalescent Constant Population
# Priors > popSize > Log Normal > initial = [Lower = 0], [Upper = 200]

# Because the K3Pu model of nucleotide substitution was not available in beast2, we used the generalized time reversible (GTR) 
# substitution model [17], which had the next lowest Bayesian information criterion (BIC) score (Delta6910.964) 
# on the basis of model testing using iq-tree [18].

# Site Model > Change JC69 to GTR

# The GTR model with prior rates having a gamma distribution with rates in [0, inf] and prior frequencies (estimated) in [0, 1]
# were applied, along with 0 proportion of invariant sites.

# Change nothing? - yes

# [Why is Rate CT estimated unchecked?]

# We used the beast2 correction for ascertainment bias, specifying the number of invariant A, C, G and T sites as 
# 758511 1449901 1444524 758336.
# Note that this must be manually added to the xml and may not appear when the xml is loaded into the BEAUti2 (version 2.6.1) software.

# [???]

# Line in the xml file:
# <constantSiteWeights id="IntegerParameter.0" spec="parameter.IntegerParameter" dimension="4" lower="0" upper="0">758511 1449901 1444524 758336</constantSiteWeights>

# We ran the Markov chain Monte Carlo (MCMC) method for 100000000 iterations, sampling every 10000th iteration.

# We verified chain convergence (by confirming multiple independent chains converged to the same posterior values) as well as good mixing 
# and an effective sample size (ESS) of greater than 200 for all parameters using Tracer (version 1.7.1) [19].
#
# A maximum clade credibility (MCC) tree was created using TreeAnnotator (version 2.6.0), with 10% of the chain discarded as burn-in,
# resulting in a posterior collection of 9000 trees.

# Instead of trying to obtain a single optimal timed phylogenetic tree from this posterior set, we sampled a collection of 
# 50 of them at random. This ensures that we capture as much diversity as possible from the beast posterior, to achieve robust
# uncertainty quantification in our subsequent analysis.









