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
  make_option(c("-", "--"), type="", default=NULL, 
              help="", metavar="")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)
print("---")
print(str(opt))

# Files ----

template_file <- paste0(opt$<path>, opt$<study accession>, "")

print("FILES:")
print(c(": ", ))


