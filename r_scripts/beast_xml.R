#!/usr/bin/env Rscript

# beast_xml.R


# Description
# 

# Arguments to script:

# Input:

# Main steps:

# Output:

# RUN:
# Rscript r_scripts/beast_xml.R -s <study_acc>  -t <template_file_location>  -f <dated_fasta_file_location> -o <output_file_location>
# Rscript r_scripts/beast_xml.R -s PRJEB7669    -t beast_xml/tb_template.xml -f fasta/                      -o beast_xml/

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
  
  make_option(c("-f", "--dated_fasta_file_location"), type="character", default=NULL, 
              help="input template xml file", metavar="character"), 
  
  make_option(c("-o", "--output_file_location"), type="character", default=NULL, 
              help="input file location", metavar="character"), 
  
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












