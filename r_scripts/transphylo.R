#!/usr/bin/env Rscript

# <file_name>.R


# Description
# 

# Arguments to script:

# Input:

# Main steps:

# Output:

# RUN:
# Rscript r_scripts/<file_name>.R
# Rscript r_scripts/<file_name>.R

# Setup ----

library(optparse)
library(ape)
library(TransPhylo)

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
  # EXAMPLE:
  # make_option(c("-t", "--template_file_name"), type="character", default=NULL,
  #             help="input template xml file", metavar="character"),
  
  # make_option(c("-s", "--study_accession"), type="character", default=NULL, 
  #             help="enter study accession code e.g. PAKISTAN", metavar="character"),
  make_option(c("-p", "--plots_path"), type="character", default="./", 
              help="path to save plots to; default is ./", metavar="character"),
  make_option(c("-t", "--tree_file"), type="character", default=NULL, 
              help="input (dated) tree file - i.e. MCC tree from BEAST and TreeAnnotator", metavar="character"), 
  make_option(c("-o", "--out_file"), type="character", default=NULL, 
              help="name of transphylo output file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)
print("---")
print(str(opt))


# Paths ----

plots_path <- opt$plots_path

# Files ----

study_accession <- opt$study_accession
tree_file <- opt$tree_file
out_file <- opt$out_file
mcc_plot_file <- paste0(plots_path, basename(tree_file), ".png")

print("ARGS:")
print(c("study_accession: ", study_accession))
print(c("tree_file: ", tree_file))
print(c("out_file: ", out_file))
print(c("mcc_plot_file", mcc_plot_file))


# Transphylo ----

# Load and parse MCC tree from BEAST

phy <- read.tree(tree_file)

png(mcc_plot_file, units="px", width=3000, height=2000, res = 300)
ptree<-ptreeFromPhylo(phy,dateLastSample=2007.964)
plot(ptree)
dev.off()













