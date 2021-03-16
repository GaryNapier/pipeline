#! /usr/bin/env Rscript

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
# Rscript r_scripts/itol_annotation.R <metadata_file>                                       <itol_location>
# Rscript r_scripts/itol_annotation.R ~/Documents/metadata/tb_data_28_01_2021_clean.csv     itol_annotations/


# Rscript itol_annotation.R ${pakistan_metadata_file}     itol_annotations/

# itol_annotation.R metadata/pakistan_metadata.csv itol/

# Setup ----

library(rvest)
library(RColorBrewer)
library(scales) # show_col
library(colorspace) # lighten()
library(gplots) # col2hex

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

# Reset
# rm(list=ls())

# Read in args
args <- commandArgs(trailingOnly=TRUE)

# Variables
# NOTE: FOR SOME BLOODY STUPID REASON I CAN'T DO THIS LATER IN THE CODE:
# subset(metadata, study_accession == study_accession)
# SO I HAVE TO RENAME THE VARIABLE study_acc INSTEAD OF study_accession
# study_acc <- args[1]
nl <- cat("\n")
alpha <- 0.8
brewer_set <- "Set1"

# Directories
# clusters_data_file_dir <- args[3]
itol_location <- args[2]

# Files and suffixes/prefixes
metadata_file <- args[1]

# clusters_data_file <- paste0(clusters_data_file_dir, study_acc, ".clusters")
# clusters_data_file <- args[3]

itol_dr_in_file <- paste0(itol_location, "itol.dr.txt")
itol_clusters_in_file <- paste0(itol_location, "itol.clusters.txt")
itol_lineage_in_file <- paste0(itol_location, "itol.lineages.txt")
itol_major_lineage_in_file <- paste0(itol_location, "itol.major_lineages.txt")

# itol_dr_out_file <- paste0(itol_location, study_acc, ".dr.txt")
# itol_clusters_out_file <- paste0(itol_location, study_acc, ".clusters.txt")
# itol_lineage_out_file <- paste0(itol_location, study_acc, ".lineages.txt")

metadata_file_base <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(metadata_file))

itol_dr_out_file <- paste0(itol_location, metadata_file_base, ".dr.txt")
itol_clusters_out_file <- paste0(itol_location, metadata_file_base, ".clusters.txt")
itol_lineage_out_file <- paste0(itol_location, metadata_file_base, ".lineages.txt")
itol_major_lins_out_file <- paste0(itol_location, metadata_file_base, ".major_lins.txt")

print("ARGUMENTS:")
# print(c("study_accession:", study_acc))
print("")
print(c("metadata_file:", metadata_file))
# print(c("clusters_data_file:", clusters_data_file))
print("")
print(c("itol_dr_in_file:", itol_dr_in_file))
print(c("itol_clusters_in_file:", itol_clusters_in_file))
print(c("itol_lineage_in_file:", itol_lineage_in_file))
print(c("itol_major_lineage_in_file", itol_major_lineage_in_file))
print("")
print(c("itol_dr_out_file:", itol_dr_out_file))
print(c("itol_clusters_out_file:", itol_clusters_out_file))
print(c("itol_lineage_out_file:", itol_lineage_out_file))

# Read in templates ----

# itol_dr_in_file
# itol_clusters_in_file
# itol_lineage_in_file

itol_dr <- readChar(itol_dr_in_file, nchars = 1e6)
itol_clusters <- readChar(itol_clusters_in_file, nchars = 1e6)
itol_lineage <- readChar(itol_lineage_in_file, nchars = 1e6) 
itol_major_lins <- readChar(itol_major_lineage_in_file, nchars = 1e6)

# Read in sample data/metadata and clean ----

metadata <- read.csv(metadata_file, stringsAsFactors = F)
# clusters_data <- read.table(clusters_data_file, header = T, stringsAsFactors = F)

# Get list of samples from clusters data
# samples <- clusters_data$id



# Drug resistance ----

print("---")
print("DR")
print("---")

# cols <- c("run_accession", "study_accession", "dr_status")
# dr_columns <- c("wgs_id", "study_accession_word", "dr_status")

# Subset rows by study accession and cols by wgs_id, (sample ID) and drugs
# dr_data <- subset(metadata, study_accession == study_acc)
# dr_data <- subset(metadata, study_accession_word == study_acc)
# dr_data <- dr_data[, cols]
dr_data <- metadata[, c("wgs_id", "dr_status")]

# Clean
# dr_data <- dr_data[!is.na(dr_data["dr_status"]), ]

# Need to convert to binary cols
# sus <- ifelse(dr_data$dr_status == "Susceptible", 1, -1)
# DR <- ifelse(dr_data$dr_status == "DR", 1, -1)
# MDR <- ifelse(dr_data$dr_status == "MDR", 1, -1)
# XDR <- ifelse(dr_data$dr_status == "XDR", 1, -1)

dr_status <- c("Sensitive", "Pre-MDR", "MDR", "Pre-XDR", "XDR", "Other")

# sus <- ifelse(dr_data$dr_status == "Sensitive", 1, -1)
# pre_mdr <- ifelse(dr_data$dr_status == "Pre-MDR", 1, -1)
# MDR <- ifelse(dr_data$dr_status == "MDR", 1, -1)
# pre_xdr <- ifelse(dr_data$dr_status == "pre-XDR", 1, -1)
# XDR <- ifelse(dr_data$dr_status == "XDR", 1, -1)
# other <- ifelse(dr_data$dr_status == "Other", 1, -1)

# dr_df <- cbind(dr_data$run_accession, sus, DR, MDR, XDR)
# dr_df <- cbind(dr_data$wgs_id, sus, DR, MDR, XDR)
# dr_df <- cbind(dr_data$wgs_id, sus, pre_mdr, MDR, pre_xdr, XDR, other)


# Define colours for each category
dr_colours <- alpha(col2hex(c("green1", "yellow2", "orange1", "red1", "black", "grey")), alpha = alpha)

# Combine in df
dr_df <- data.frame(dr_status, dr_colours)

# Merge with dr data
dr_df <- merge(dr_data, dr_df, by = "dr_status", all.x = T, sort = F)

# Clean
dr_df <- dr_df[, c("wgs_id", "dr_colours", "dr_status")]

# DR legend

# Legend shapes - repeat "1" for as many unique clusters and intersperse with "\t"
leg_shapes_dr <- paste0("LEGEND_SHAPES\t", 
                     paste0(paste0(rep(1, length(dr_status)), "\t"), collapse = ""))

itol_dr <- gsub("#LEGEND_SHAPES 1 1 2 2", leg_shapes_dr, itol_dr)

# Legend colours - same for labels but just the unique cols
leg_cols_dr <- paste0("LEGEND_COLORS\t", 
                   paste0(paste0(dr_colours, "\t"), collapse = ""))

itol_dr <- gsub("#LEGEND_COLORS #ff0000 #00ff00 rgba\\(0,255,0,0.5\\) #0000ff", leg_cols_dr, itol_dr)

# Legend labels - same for labels but just the unique values
leg_labels_dr <- paste0("LEGEND_LABELS\t", 
                     paste0(paste0(dr_status, "\t"), collapse = ""))

itol_dr <- gsub("#LEGEND_LABELS value1 value2 value3 value4", 
                              leg_labels_dr, itol_dr)

# Write dataframe to template

# Write the template out under the new file name for the study accession
write.table(itol_dr, file = itol_dr_out_file, sep="\t",
            row.names=F, col.names=F, quote = F)

# Append the drug resistance data to the template
write.table(dr_df, file = itol_dr_out_file,
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)



# Lineage ----

print("---")
print("LINEAGE")
print("---")

# Example:
# <id>        <"range"> <col>     <lineage>
# ERR2446223	range	    #FF1919E6	1
# SRR2100428	range	    #FF1919E6	1

# Subset data
lin_data <- metadata[, c("wgs_id", "main_lineage")]

# Get unique lins and split into numeric and 'character' strains - 'character' strains should be animal ones
uniq_lins <- sort(unique(lin_data[!is.na(lin_data[, "main_lineage"]), "main_lineage"]))
uniq_lins_num <- grep("\\d", uniq_lins, value = T)
uniq_lins_char <- grep("[^0-9]", uniq_lins, value = T)

# Get unique cols per lineage, per category
n_cols_lin_num <- length(uniq_lins_num)
# lin_colours_num <- rainbow(n_cols_lin_num, alpha = alpha)
lin_colours_num <- brewer.pal(n = n_cols_lin_num, name = brewer_set)
lin_colours_char <- col2hex(c("grey10", "grey40", "grey80"))
if(length(uniq_lins_char) > 0){
  lin_colours_char <- lin_colours_char[1:length(lin_colours_char)]
  lin_colours <- c(lin_colours_num, lin_colours_char) 
} else {
  lin_colours <- lin_colours_num
}

# Make df for unique lins and cols 
lin_col_df <- data.frame(main_lineage = uniq_lins, col = lin_colours)

# Merge
lin_data <- merge(lin_data, lin_col_df, by = "main_lineage")

# Clean - subset and order columns and add "range" column
# lin_data$range <- rep("range", nrow(lin_data))
# lin_data <- lin_data[, c("wgs_id", "range", "col", "main_lineage")]
lin_data <- lin_data[, c("wgs_id", "col", "main_lineage")]

# Lineage legend

# Legend shapes - repeat "1" for as many unique clusters and intersperse with "\t"
leg_shapes_lineage <- paste0("LEGEND_SHAPES\t",
                        paste0(paste0(rep(1, length(uniq_lins)), "\t"), collapse = ""))

itol_lineage <- gsub("#LEGEND_SHAPES 1 1 2 2", leg_shapes_lineage, itol_lineage)

# Legend colours - same for labels but just the unique cols
leg_cols_lineage <- paste0("LEGEND_COLORS\t",
                      paste0(paste0(lin_colours, "\t"), collapse = ""))

itol_lineage <- gsub("#LEGEND_COLORS #ff0000 #00ff00 rgba\\(0,255,0,0.5\\) #0000ff", leg_cols_lineage, itol_lineage)

# Legend labels - same for labels but just the unique values
leg_labels_lineage <- paste0("LEGEND_LABELS\t",
                        paste0(paste0(uniq_lins, "\t"), collapse = ""))

itol_lineage <- gsub("#LEGEND_LABELS value1 value2 value3 value4",
                leg_labels_lineage, itol_lineage)

# Write the template out under the new file name for the study accession
write.table(itol_lineage, file = itol_lineage_out_file, sep="\t",
            row.names=F, col.names=F, quote = F)

# Append the lin data to the template
write.table(lin_data, file = itol_lineage_out_file,
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)



# Major lineages ----

print("---")
print("MAJOR LINEAGES")
print("---")

# Have to take major lins from sublins first (rather than adding the major lin to metadata first) 
# because of the sodding animal lineages

major_lin_data <- metadata[, c("wgs_id", "sub_lineage")]
sub_lins <- sort(unique(major_lin_data$sub_lineage))
# Remove animal strains
animal <- c("M.bovis", "M.caprae", "M.orygis")
sub_lins <- sub_lins[!(sub_lins %in% animal)]
major_lins <- unique(substr(sub_lins, 1, 3)) 
major_lins_split <- split(major_lins, substr(major_lins, 1, 1))
all_major_lins_cols <- vector()
for(i in seq(major_lins_split)){
  # all_major_lins_cols <- append(all_major_lins_cols, rainbow(length(major_lins_split[[i]]), alpha = alpha))
  if (length(major_lins_split[[i]]) == 1){
    all_major_lins_cols <- append(all_major_lins_cols, lin_colours[i])
  }else{
    n <- length(major_lins_split[[i]])
    cols <- brewer.pal(n = n, name = brewer_set)
    all_major_lins_cols <- append(all_major_lins_cols, cols[1:n])
  }
}
# Add animal strain colours to major lin colours
if(length(uniq_lins_char) > 0){
  all_major_lins_cols <- c(all_major_lins_cols, lin_colours_char)
}

# Add major lin column to metadata for merge
metadata$major_lineage <- substr(metadata$sub_lineage, 1, 3)

# Add animal strains back in
metadata$major_lineage <- ifelse(metadata$sub_lineage %in% animal, metadata$sub_lineage, metadata$major_lineage)

# Re-define major lins data now with the animal strains
major_lin_data <- metadata[, c("wgs_id", "major_lineage")]

# Re-define major_lins now with animal strains...

# major_lins <- c(major_lins, animal)
major_lins <- c(major_lins, uniq_lins_char)

# Put major lins and colours in dataframe
major_lins_df <- data.frame(major_lineage = major_lins, col = all_major_lins_cols)

# Merge
major_lins_df <- merge(major_lin_data, major_lins_df, by = "major_lineage", all.x = T, sort = F)

# Clean
major_lins_df <- major_lins_df[, c("wgs_id", "col", "major_lineage")]


# Major lineage legend


# Legend shapes - repeat "1" for as many unique clusters and intersperse with "\t"
leg_shapes_major_lineage <- paste0("LEGEND_SHAPES\t",
                             paste0(paste0(rep(1, length(major_lins)), "\t"), collapse = ""))

itol_major_lins <- gsub("#LEGEND_SHAPES 1 1 2 2", leg_shapes_major_lineage, itol_major_lins)

# Legend colours - same for labels but just the unique cols
leg_cols_major_lineage <- paste0("LEGEND_COLORS\t",
                           paste0(paste0(all_major_lins_cols, "\t"), collapse = ""))

itol_major_lins <- gsub("#LEGEND_COLORS #ff0000 #00ff00 rgba\\(0,255,0,0.5\\) #0000ff", leg_cols_major_lineage, itol_major_lins)

# Legend labels - same for labels but just the unique values
leg_labels_major_lineage <- paste0("LEGEND_LABELS\t",
                             paste0(paste0(major_lins, "\t"), collapse = ""))

itol_major_lins <- gsub("#LEGEND_LABELS value1 value2 value3 value4",
                     leg_labels_major_lineage, itol_major_lins)

# Write the template out under the new file name for the study accession
write.table(itol_major_lins, file = itol_major_lins_out_file, sep="\t",
            row.names=F, col.names=F, quote = F)

# Append the lin data to the template
write.table(major_lins_df, file = itol_major_lins_out_file,
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)





# Clusters ----

print("---")
print("SUBLINEAGES")
print("---")

# Subset data
clusters_data <- metadata[, c("wgs_id", "sub_lineage")]

# Set colours for each cluster
n_clusts_total <- length(unique(clusters_data$sub_lineage))

## Wrap in alpha function
# cluster_cols <- scales::alpha(brewer.pal(n = n_clusts, name = "Dark2")[1:n_clusts], alpha = alpha)

# Define colours for all sublins (!!!)
# Loop through the major lins and assign a new rainbow of length number of sublins in the major lin
sub_lins <- sort(unique(clusters_data$sub_lineage))
# Drop animal lineages
sub_lins <- sub_lins[!(sub_lins %in% animal)]
sublin_split <- split(sub_lins, substr(sub_lins, 1, 3))
all_sublin_cols <- vector()
for(i in seq(sublin_split)){
  all_sublin_cols <- append(all_sublin_cols, rainbow(length(sublin_split[[i]]), alpha = alpha))
}
# Append animal
all_sublin_cols <- c(all_sublin_cols, lin_colours_char)

# Make df for unique clusters and cols
clust_col_df <- data.frame(sub_lineage = sort(unique(clusters_data$sub_lineage)), col = all_sublin_cols)

# Merge with cluster data
clusters_data <- merge(clusters_data, clust_col_df, by = "sub_lineage")

# Clean - rearrange cols
# clusters_data <- clusters_data[, c("id", "col", "cluster")]
clusters_data <- clusters_data[, c("wgs_id", "col", "sub_lineage")]

# Do legend shapes and labels

# Legend shapes - repeat "1" for as many unique clusters and intersperse with "\t"
leg_shapes <- paste0("LEGEND_SHAPES\t",
                     paste0(paste0(rep(1, n_clusts_total), "\t"), collapse = ""))

itol_clusters <- gsub("#LEGEND_SHAPES 1 1 2 2", leg_shapes, itol_clusters)


# Legend colours - same for labels but just the unique cols
leg_cols <- paste0("LEGEND_COLORS\t",
                   paste0(paste0(clust_col_df$col, "\t"), collapse = ""))

itol_clusters <- gsub("#LEGEND_COLORS #ff0000 #00ff00 rgba\\(0,255,0,0.5\\) #0000ff",
                              leg_cols, itol_clusters)

# Legend labels - same for labels but just the unique values
leg_labels <- paste0("LEGEND_LABELS\t",
                     paste0(paste0(unique(clusters_data$sub_lineage), "\t"), collapse = ""))

itol_clusters <- gsub("#LEGEND_LABELS value1 value2 value3 value4",
                              leg_labels, itol_clusters)

# Write the template out under the new file name for the study accession
write.table(itol_clusters, file = itol_clusters_out_file, sep="\t",
            row.names=F, col.names=F, quote = F)

# Append the clusters data to the template
write.table(clusters_data, file = itol_clusters_out_file,
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)











