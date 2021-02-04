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
# Rscript r_scripts/merge_rep_regions_and_dr_bed.R --rep_regions_bed_file <rep_reg_bed_file>              --gff_file <gff_file>                                 --dr_file <dr_file>                             --out_location <out_location>
# Rscript r_scripts/merge_rep_regions_and_dr_bed.R --rep_regions_bed_file ~/refgenome/excluded_loci.bed   --gff_file ~/refgenome/MTB-h37rv_asm19595v2-eg18.gff  --dr_file ~/metadata/drug_resistance_genes.csv  --out_location ~/refgenome/

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
  make_option(c("-r", "--rep_regions_bed_file"), type="character", default=NULL,
              help="Enter name of repetitive regions bed file", metavar="character"),
  make_option(c("-g", "--gff_file"), type="character", default=NULL,
              help="Enter name of gff annotations file", metavar="character"),
  make_option(c("-d", "--dr_file"), type="character", default=NULL,
              help="Enter name of drug resistance genes file", metavar="character"),
  make_option(c("-o", "--out_location"), type="character", default="./",
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
rep_regions <- read.delim(rep_regions_bed_file, header = F)
gff <- data.frame(readGFF(gff_file))
dr <- read.csv(dr_file, header = T, stringsAsFactors = F)
out_file <- paste0(out_location, "excluded_loci_rep_and_dr.bed")


subset(gff, gene_id == "EBG00000313325")


print("REP REGIONS BED:")
head(rep_regions)
print("GFF")
head(gff)
print("DR")
head(dr)
print("OUT FILE")
head(out_file)

# Clean data ----

# gff data

# Remove these two cols. Not sure what they are and seem very messy
gff <- drop_cols(gff, c("Parent", "Alias"))

# Only need the genes
# gff <- subset(gff, type == "gene")

# Tidy up Name col - if name is unavailable then fill in with gene ID (to match drug resistance genes file)
gff$Name <- ifelse(is.na(gff$Name), gff$gene_id, gff$Name)

# Filter based on drug resistance genes
dr_genes <- unique(dr$Gene)
gff <- subset(gff, Name %in% dr_genes)

# Retain "gene" - anything that is actually "gene" or anything that is e.g. rRNA_gene
gff <- gff[grepl("gene", gff$type), ]

# Sort out description col
gff$description <- gsub("Probable ", "", gff$description)
gff$description <- gsub("Possible ", "", gff$description)
gff$description <- gsub(" ", "_", gff$description)

gff$description <- ifelse(is.na(gff$description), gff$biotype, gff$description)


# rep_regions bed file

# Name cols based on gff data
names(rep_regions) <- c("seqid", "start", "end", "gene_id", "description")


# Merge (rbind) data ----

# Subset cols
dr_df <- gff[, c("seqid", "start", "end", "gene_id", "description")]

# rbind
df <- rbind(rep_regions, dr_df)

# sort
df <- df[order(df$start), ]

# Remove duplicates
df <- df[!(duplicated(df[, "gene_id"]) ), ]




df[ df[, "gene_id"] %in%  df[duplicated(df[, "gene_id"]), ]$gene_id, ]


subset(gff, gene_id == "Rv0605")

subset(rep_regions, gene_id == "Rv0605")




print("DF")
head(df)
print("---")
str(df)

print("GENES IN DR LIST THAT ARE NOT IN DF LIST")
df$gene_id[!(dr_genes %in% df$gene_id)]

































