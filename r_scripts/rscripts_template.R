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

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

# Arguments ----

option_list = list(
  # EXAMPLE:
  # make_option(c("-t", "--template_file_name"), type="character", default=NULL,
  #             help="input template xml file", metavar="character"),
  
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

template_file <- opt$template_file_name

print("FILES:")
print(c(": ", ))

