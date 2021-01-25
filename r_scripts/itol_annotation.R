#!/usr/bin/env Rscript

# itol_annotation.R

# Annotate trees for each study accession (project) with:
# - Drug resistance (from metadata)
# - Clusters based on min SNP distance - output of r_scripts/transmission_clusters.R -> metadata/<study_accession>.clusters)
# - Lineage - output of tbprofiler

# Arguments to script:
# Study accession, itol template files location, DR metadata file name (full path), clusters file location (dir only), lineage file name (full path), output_location

# Input:
# itol template files, DR metadata file, clusters file, lineage file

# Main steps:
# Filter DR data by study accession, clean and parse into dataframe ready for itol template
# Get cluster membership of samples and onvert to colour scheme - parse into dataframe ready for itol template
# Repeat for lineage

# Output:
# Three itol files, one for each type of annotation

# RUN:
# Rscript r_scripts/itol_annotation.R <study_accession> <metadata_file>                                            <clusters_data_file>                   <itol_location>
# Rscript r_scripts/itol_annotation.R PRJEB7669         ~/Documents/metadata/tb_data_collated_28_10_2020_clean.csv metadata/PRJEB7669.clusters            itol_annotations/ 

# Setup ----

library(rvest)
library(RColorBrewer)
library(scales)

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

# Read in args
args <- commandArgs(trailingOnly=TRUE)

# Variables
# NOTE: FOR SOME BLOODY STUPID REASON I CAN'T DO THIS LATER IN THE CODE:
# subset(metadata, study_accession == study_accession)
# SO I HAVE TO RENAME THE VARIABLE study_acc INSTEAD OF study_accession
study_acc <- args[1]
nl <- cat("\n")
alpha <- 0.7

# Directories
# clusters_data_file_dir <- args[3]
itol_location <- args[4]

# Files and suffixes/prefixes
metadata_file <- args[2]

# clusters_data_file <- paste0(clusters_data_file_dir, study_acc, ".clusters")
clusters_data_file <- args[3]

itol_dr_in_file <- paste0(itol_location, "itol.dr.txt")
itol_clusters_in_file <- paste0(itol_location, "itol.clusters.txt")
itol_lineage_in_file <- paste0(itol_location, "itol.lineages.txt")

itol_dr_out_file <- paste0(itol_location, study_acc, ".dr.txt")
itol_clusters_out_file <- paste0(itol_location, study_acc, ".clusters.txt")
itol_lineage_out_file <- paste0(itol_location, study_acc, ".lineages.txt")

nl
print("ARGUMENTS:")
nl
print(c("study_accession:", study_acc))
print("")
print(c("metadata_file:", metadata_file))
print(c("clusters_data_file:", clusters_data_file))
print("")
print(c("itol_dr_in_file:", itol_dr_in_file))
print(c("itol_clusters_in_file:", itol_clusters_in_file))
print(c("itol_lineage_in_file:", itol_lineage_in_file))
print("")
print(c("itol_dr_out_file:", itol_dr_out_file))
print(c("itol_clusters_out_file:", itol_clusters_out_file))
print(c("itol_lineage_out_file:", itol_lineage_out_file))

# Read in templates ----

# itol_dr_in_file
# itol_clusters_in_file
# itol_lineage_in_file

itol_dr_in_file <- readChar(itol_dr_in_file, nchars = 1e6)
itol_clusters_in_file <- readChar(itol_clusters_in_file, nchars = 1e6)
itol_lineage_in_file <- readChar(itol_lineage_in_file, nchars = 1e6)

# Read in sample data/metadata and clean ----

metadata <- read.csv(metadata_file, stringsAsFactors = F)
clusters_data <- read.table(clusters_data_file, header = T, stringsAsFactors = F)

# Get list of samples from clusters data
samples <- clusters_data$id

# Drug resistance ----

cols <- c("run_accession", "study_accession", "dr_status")

# Subset rows by study accession and cols by run_accession (sample ID), study_accession and drugs
dr_data <- subset(metadata, study_accession == study_acc)
dr_data <- dr_data[, cols]

# Clean
dr_data <- dr_data[!is.na(dr_data["dr_status"]), ]

# Need to convert to binary cols
sus <- ifelse(dr_data$dr_status == "Susceptible", 1, -1)
DR <- ifelse(dr_data$dr_status == "DR", 1, -1)
MDR <- ifelse(dr_data$dr_status == "MDR", 1, -1)
XDR <- ifelse(dr_data$dr_status == "XDR", 1, -1)

dr_df <- cbind(dr_data$run_accession, sus, DR, MDR, XDR)

# Write dataframe to template

# Write the template out under the new file name for the study accession
write.table(itol_dr_in_file, file = itol_dr_out_file, sep="\t",
            row.names=F, col.names=F, quote = F)

# Append the drug resistance data to the template
write.table(dr_df, file = itol_dr_out_file,
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)


# Clusters ----

# From itol template:

#Examples:
#assign a red colored strip to leaf 9606, with label 'Human' (label is displayed in the mouseover popups)
# <id> <col> <label>
#9606 #ff0000 Human

# Set colours for each cluster
n_cols_clust <- length(unique(clusters_data$cluster))
# Get colours - note: In brewer.pal minimal value for n is 3, so have to subset with [] if less than 3
# Wrap in alpha function
cluster_cols <- scales::alpha(brewer.pal(n = n_cols_clust, name = "Dark2")[1:n_cols_clust], alpha = alpha)

# Make df for unique clusters and cols 
clust_col_df <- data.frame(cluster = unique(clusters_data$cluster), col = cluster_cols)

# Merge with cluster data 
clusters_data <- merge(clusters_data, clust_col_df, by = "cluster")

# Clean - rearrange cols
clusters_data <- clusters_data[, c("id", "col", "cluster")]

# Do legend shapes and labels

#LEGEND_TITLE Dataset_legend
#LEGEND_POSITION_X 100
#LEGEND_POSITION_Y 100
#LEGEND_SHAPES 1 1 2 2
#LEGEND_COLORS #ff0000 #00ff00 rgba(0,255,0,0.5) #0000ff
#LEGEND_LABELS value1 value2 value3 value4
#LEGEND_SHAPE_SCALES 1 1 0.5 1

# LEGEND_TITLE	Clusters
# #LEGEND_POSITION_X 100
# #LEGEND_POSITION_Y 100
# LEGEND_SHAPES	1
# LEGEND_COLORS	#1B9E77E6
# LEGEND_LABELS	1
# #LEGEND_SHAPE_SCALES 1 1 0.5 1


# Legend shapes - repeat "1" for as many unique clusters and intersperse with "\t"
leg_shapes <- paste0("LEGEND_SHAPES\t", 
                     paste0(paste0(rep(1, length(unique(clusters_data$cluster))), "\t"), collapse = ""))

itol_clusters_in_file <- gsub("#LEGEND_SHAPES 1 1 2 2", leg_shapes, itol_clusters_in_file)


# Legend colours - same for labels but just the unique cols
leg_cols <- paste0("LEGEND_COLORS\t", 
                     paste0(paste0(unique(clusters_data$col), "\t"), collapse = ""))

itol_clusters_in_file <- gsub("#LEGEND_COLORS #ff0000 #00ff00 rgba\\(0,255,0,0.5\\) #0000ff", 
                              leg_cols, itol_clusters_in_file)

# Legend labels - same for labels but just the unique values
leg_labels <- paste0("LEGEND_LABELS\t", 
                     paste0(paste0(unique(clusters_data$cluster), "\t"), collapse = ""))

itol_clusters_in_file <- gsub("#LEGEND_LABELS value1 value2 value3 value4", 
                              leg_labels, itol_clusters_in_file)

# Write the template out under the new file name for the study accession
write.table(itol_clusters_in_file, file = itol_clusters_out_file, sep="\t",
            row.names=F, col.names=F, quote = F)

# Append the clusters data to the template
write.table(clusters_data, file = itol_clusters_out_file,
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)

# Lineage ----

# Example:
# <id> <"range"> <col> <lineage>
# ERR2446223	range	#FF1919E6	1
# SRR2100428	range	#FF1919E6	1

cols <- c("run_accession", "study_accession", "lineage")

# Subset data
lin_data <- subset(metadata, study_accession == study_acc)
lin_data <- lin_data[, cols]

uniq_lins <- unique(lin_data[!is.na(lin_data[, "lineage"]), "lineage"])

# Get unique cols per lineage
n_cols_lin <- length(uniq_lins)
lin_cols <- rainbow(n_cols_lin, alpha = alpha)

# Make df for unique lins and cols 
lin_col_df <- data.frame(lineage = uniq_lins, col = lin_cols)

# Merge
lin_data <- merge(lin_data, lin_col_df, by = "lineage")

# Clean - subset and order columns and add "range" column
lin_data$range <- rep("range", nrow(lin_data))
lin_data <- lin_data[, c("run_accession", "range", "col", "lineage")]

# Write the template out under the new file name for the study accession
write.table(itol_lineage_in_file, file = itol_lineage_out_file, sep="\t",
            row.names=F, col.names=F, quote = F)

# Append the lin data to the template
write.table(lin_data, file = itol_lineage_out_file,
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)













 