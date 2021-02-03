#!/usr/bin/env Rscript

# merge_rep_regions_and_dr_bed.R

# Description
# There is a pre-existing bed file of repetitive regions (including PE/PPE regions) - excluded_loci.bed.

# Chromosome	33582	  33794	  Rv0031	remnant_of_A_transposase
# Chromosome	103710	104663	Rv0094c	50bp_duplicated
# Chromosome	104805	105215	Rv0095c	50bp_duplicated
# Chromosome	105324	106715	Rv0096	PPE_family_protein
# Chromosome	131382	132872	Rv0109	PE-PGRS_family_protein
# Chromosome	149533	150996	Rv0124	PE-PGRS_family_protein

# The start and end positions of genes with known drug resistance mutations need to be appended to this list (and sorted), 
# to make one big bed file of exclusion loci which is passed to shell_scripts/variant_filtering.sh

# Arguments to script:
# repetitive regions bed file, 
# gff annotations file with ALL genes and positions, 
# drug resistance genes file

# Input:
# repetitive regions bed file, 
# gff annotations file with ALL genes and positions, 
# drug resistance genes file

# Main steps:
# Parse and clean gff file
# Subset based on drug resistance genes file
# Combine with repetitive regions file 

# Output:
# One bed file of concatenated repetitive regions and drug resistance regions

# RUN:
# Rscript r_scripts/merge_rep_regions_and_dr_bed.R --rep_regions_bed_file <rep_reg_bed_file>              --gff_file <gff_file>                                 --dr_file <dr_file>                         --out_location <out_location>
# Rscript r_scripts/merge_rep_regions_and_dr_bed.R --rep_regions_bed_file ~/refgenome/excluded_loci.bed   --gff_file ~/refgenome/MTB-h37rv_asm19595v2-eg18.gff  --dr_file ~/metadata/drug_resistance_genest --out_location ~/refgenome/

# Setup ----

library(optparse)
library(rtracklayer)

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

drop_cols <- function(x, cnames){
  x[, -which(names(x) %in% cnames)]
}

# Arguments ----

option_list = list(
  # EXAMPLE:
  # make_option(c("-t", "--template_file_name"), type="character", default=NULL,
  #             help="input template xml file", metavar="character"),
  
  make_option(c("-r", "--rep_regions_bed_file"), type="str", default=NULL, 
              help="Enter name of repetitive regions bed file", metavar="character"),
  make_option(c("-g", "--gff_file"), type="str", default=NULL, 
              help="Enter name of gff annotations file", metavar="character"),
  make_option(c("-d", "--dr_file"), type="str", default=NULL, 
              help="Enter name of drug resistance genes file", metavar="character"),
  make_option(c("-o", "--out_location"), type="str", default="./", 
              help="Enter location to save output file. Default is current location ./", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)

# Files ----

rep_regions_bed_file <- opt$rep_regions_bed_file
gff_file <- opt$gff_file
dr_file <- opt$dr_file
out_location <- opt$out_location

# Read in data
rep_regions <- read.delim(rep_regions_bed_file)
gff <- data.frame(readGFF(gff_file))
dr <- read.csv(dr_file, header = T, stringsAsFactors = F)
out_file <- paste0(out_location, "rep_reg_and_dr_reg_exclude.bed")


head(rep_regions)
head(gff)
head(dr)
head(out_file)


# Parse and clean gff file
# x <- data.frame(readGFF("test"))
# 
# x <- drop_cols(x, c("Parent", "Alias"))
# 
# x <- subset(x, type == "gene")
# 
# x$Name <- ifelse(is.na(x$Name), x$gene_id, x$Name)




















