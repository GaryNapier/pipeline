

#!/usr/bin/env Rscript

# data_wrangling.R

# Description
# Wrangle data for tables and figures

# Input:

# Main steps:

# Output:


# Libraries ----

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(janitor)
library(ape)
library(scales)

# Functions ----

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

clear <- function(){rm(list=ls(envir = .GlobalEnv),envir = .GlobalEnv); }

to_table <- function(x){
  x <- x %>% adorn_totals(c("row", "col")) %>%
    adorn_percentages(c("row")) %>%
    adorn_pct_formatting(digits = 2)
  formatted_ns <- attr(x, "core") %>% # extract the tabyl's underlying Ns
    adorn_totals(c("row", "col")) %>% # to match the data.frame we're appending to
    dplyr::mutate_if(is.numeric, format, big.mark = ",")
  x %>% adorn_ns(position = "front", ns = formatted_ns)
}

# Setup ----

study_accessions <- c("PAKISTAN", "THAILAND")
id_col <- "wgs_id"
study_acc_col <- "study_accession_word"

# Paths ----

trans_path <- "~/Documents/transmission/"
metadata_path <- "~/Documents/metadata/"
fasta_path <- "fasta/"

setwd(trans_path)


# File names ----

metadata_file <- paste0(metadata_path, "tb_data_28_01_2021_clean.csv")
fasta_file <- paste0(fasta_path, "THAILAND_TEST.clust_1.dated.fa")

print("FILES:")
print(c("metadata_file: ", metadata_file))
print(c("fasta_file: ", fasta_file))


# Read in files ----

metadata <- read.csv(metadata_file, stringsAsFactors = F)
fasta <- read.FASTA(fasta_file)


# Basic stats ----

metadata <- subset(metadata, study_accession_word %in% study_accessions)


# DR summary

# Desired order of DR cols
dr_vals <- c("Sensitive", "Pre-MDR", "MDR", "Pre-XDR", "XDR", "Other", "NA") 

# Pivot table
dr <- reshape2::dcast(metadata, study_accession_word ~ dr_status, value.var = id_col, fun.aggregate = length)
# Re-order the cols
dr <- dr[, c(study_acc_col, dr_vals)]
# Re-name
names(dr) <- c("Study", dr_vals)
# Add totals and percentages
to_table(dr)




# Plot fasta alignments ----

par(mar=c(5, 10, 4, 2))
image(fasta, col = c(rainbow(4, alpha = 0.7), "black"))





























