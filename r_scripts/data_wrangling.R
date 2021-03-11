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

# Variables ----

study_accessions <- c("PAKISTAN", "THAILAND")
id_col <- "wgs_id"
study_acc_col <- "study_accession_word"
study_acc_table_name <- "Study"

# Paths ----

trans_path <- "~/Documents/transmission/"
metadata_path <- "~/Documents/metadata/"
fasta_path <- "fasta/"
plots_path <- "plots/"

setwd(trans_path)


# File names ----

metadata_file <- paste0(metadata_path, "tb_data_18_02_2021_clean.csv")
fasta_file <- paste0(fasta_path, "THAILAND_TEST.clust_1.dated.fa")
fasta_files <- paste0(fasta_path, paste0(study_accessions, ".filt.val.gt.g.snps.fa"))

print("FILES:")
print(c("metadata_file: ", metadata_file))
print(c("fasta_files: ", fasta_files))


# Read in files ----

metadata <- read.csv(metadata_file, stringsAsFactors = F)
fasta <- read.FASTA(fasta_file)
fastas <- list()
for(i in seq(fasta_files)){
  fastas[[i]] <- read.FASTA(fasta_files[i])
}


# Subset/clean ----

metadata <- subset(metadata, study_accession_word %in% study_accessions)


# DR summary ----

# Desired order of DR cols
dr_vals <- c("Sensitive", "Pre-MDR", "MDR", "Pre-XDR", "XDR", "Other", "NA") 

# Pivot table
dr <- reshape2::dcast(metadata, study_accession_word ~ dr_status, value.var = id_col, fun.aggregate = length)
# Drop dr vals if not in table
dr_vals <- dr_vals[dr_vals %in% names(dr)]
# Re-order the cols
dr <- dr[, c(study_acc_col, dr_vals)]
# Re-name
names(dr) <- c(study_acc_table_name, dr_vals)
# Add totals and percentages
to_table(dr)


# Lineages summary ---- 

lin_table <- reshape2::dcast(metadata, major_lineage ~ study_accession_word,
                             value.var = id_col, fun.aggregate = length)

# Re-name
# names(lin_table) <- c(study_acc_table_name, names(lin_table)[2:length(names(lin_table))])
names(lin_table) <- c("Major lineage", names(lin_table)[2:length(names(lin_table))])
to_table(lin_table)


# Dates ----

table(metadata$study_accession_word, metadata$year)


# Plot fasta alignments ----

par(mar=c(5, 10, 4, 2))
# image(fasta, col = c(rainbow(4, alpha = 0.7), "black"))
i <- 0
for(study in study_accessions){
  i <- i+1
  plot_file <- paste0(plots_path, study, ".fasta_block.png")
  if(length(fasta[[i]][[1]]) < 1000){
    png(plot_file, units="px", width=3000, height=2000, res = 300)
    image(fasta[[i]], col = c(rainbow(4, alpha = 0.7), "black"), cex = 0.5)
    dev.off()
  }
}


# See:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7725332/ and sup materials:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7725332/bin/mgen-6-450-s001.pdf
# for plotting/tables ideas

# Root-to-tip genetic distance of the isolates to their sampling time ---- 




# Trajectories of model parameters during the TransPhylo MCMC run

# MDS plot of (a sample from) the posterior collection of phylogenetic trees for 4 independent
# BEAST chains, through which we find evidence of a unimodal posterior.

# Table S1 shows the results of a sensitivity analysis to a number of the model priors.

# We assess the impact on several key outcomes of varying the mean generation time, the mean sampling time
# and the within-host coalescent time unit Neg. We vary each quantity by +/-50% compared to the
# main analysis, to explore the result of a significant change. We find that all measured outcomes
# are fairly robust to changes in Neg

# Figure S3 shows the estimated number of unsampled cases

# Figure S5 shows the receiver18 operator characteristic curves for the classification tasks.

# Figure S6 shows the feature importance and partial dependence on feature age for the second machine learning task: classifying
# 24 long/short generation time of hosts.


















