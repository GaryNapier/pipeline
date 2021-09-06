#!/usr/bin/env Rscript

library(ape)
library(optparse)

option_list = list(
    make_option(c("-i", "--infile"), type="character", default=NULL,
              help="name of nexus file", metavar="character"),
    make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="name of newick output file", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)
print("---")
print(str(opt))

nexus_file <- opt$infile
newick_file <- opt$outfile

# Read in nexus file 
nexus_tree <- read.nexus(nexus_file)

# Write out as newick
write.tree(nexus_tree, file = newick_file)

