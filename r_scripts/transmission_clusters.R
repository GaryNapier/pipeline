#!/usr/bin/env Rscript

# transmission_clusters.R

# Description:
# Take the 'major lineage' - the first decimal place of the sublineage e.g. 1.2.1.1 -> 1.2 for each sample
# Clean animal strains - "M.bovis", "M.caprae", "M.orygis" - these stay the same
# Remove samples for which there are only one or two samples in that major lineage - cannot pass one or two samples to Beast

# Arguments to script:
# Study accession, metadata file, location to save file

# Input:
# Metadata

# Main steps:
# Read in data
# Define major lineage
# Clean animal
# Remove one/two sample lineages from list
# Save to file

# Output:
# Tab-sep list of samples in the study accession and their major lineage:
# wgs_id          major_lineage
# SAMEA2533789    2.2
# SAMEA2534474    2.2
# SRR5709823      2.2
# ...etc

# RUN:
# Assume run from ~/transmission:
# Rscript r_scripts/transmission_clusters.R <study_accession>   <metadata_file>                             <output_clusters_dir>
# Rscript r_scripts/transmission_clusters.R PRJEB7669           ~/metadata/tb_data_28_01_2021_clean.csv     metadata/             


# Setup ----

# Read in args
args <- commandArgs(trailingOnly=TRUE)

# Variables
study_accession <- args[1]
# threshold <- as.numeric(args[2])
nl <- cat("\n")

# Directories
output_clusters_dir <- args[3]
# output_plot_dir <- args[6]

# Files and suffixes/prefixes
# dm_file <- args[3]
# id_file <- args[4]
metadata_file <- args[2]
output_clusters_file <- paste0(output_clusters_dir, study_accession, ".clusters")
# output_plot_file <- paste0(output_plot_dir, study_accession, ".clusters.png")

nl
print("ARGUMENTS:")
nl
print(c("Study accession:", study_accession))
# print(c("Threshold:", threshold))
# print(c("DM file:", dm_file))
# print(c("ID file:", id_file))
# print(c("Output clusters directory:", output_clusters_dir))
# print(c("Output plot directory:", output_plot_dir))
print(c("Metadata file: ", metadata_file))
print(c("Cluster output file name:", output_clusters_file))
# print(c("Plot output file name:", output_plot_file))

# Read in data and clean ----


metadata <- read.csv(metadata_file, header = T, stringsAsFactors = F)

# Make "major lineage" column - first decimal place of sublineage
metadata$major_lineage <- substr(metadata$sub_lineage, 1, 3)

# Clean animal strains
animal <- c("M.bovis", "M.caprae", "M.orygis")
metadata$major_lineage <- ifelse(metadata$sub_lineage %in% animal, metadata$sub_lineage, metadata$major_lineage)

# Subset on study accession - n.b. have to change the name of the variable
study_acc <- study_accession
metadata <- subset(metadata, study_accession_word == study_acc)

# Get just the id and major linege columns
metadata <- metadata[, c("wgs_id", "major_lineage")]

# Get the samples for which there are only 1 or two sampples
single_lins <- names(which(table(metadata$major_lineage) < 3))

# Remove these 
metadata <- subset(metadata, !(major_lineage %in% single_lins))

# Write to file
write.table(metadata, file = output_clusters_file, quote = F, row.names = F, sep = "\t")




# ---------------
# ---------------
# OLD OLD OLD
# ---------------
# ---------------


# Find samples within a SNP distance of a given threshold and cluster together.

# Example:
# Given this matrix:
#    A  B  C  D  E  F  G  H  I J
# A  0  0  0  0  0  0  0  0  0 0
# B 51  0  0  0  0  0  0  0  0 0
# C  3  1  0  0  0  0  0  0  0 0
# D 71 84 15  0  0  0  0  0  0 0
# E 44 35 22 65  0  0  0  0  0 0
# F 58 27 69 48  5  0  0  0  0 0
# G 51 48 64 24 50 11  0  0  0 0
# H 56 48 49 59 52 26 62  0  0 0
# I 30  2 54 32 78 22 21 79  0 0
# J 62 32 75 13 48 13 56 97 62 0

# If the threshold distance is 12, there are two clusters of 4 and 3 samples:
# A-C-B-I, E-F-G

# See appendix at the end of this script for examples

# From Sobkowiak paper:
# https://www.microbiologyresearch.org/docserver/fulltext/mgen/6/4/mgen000361.pdf?expires=1607610831&id=id&accname=guest&checksum=02AB64D1B492C16EAF2AE14908348712
# Methods > Phylogenetic analysis

# "Initially, broad putative clusters of cases were constructed for
# application of the Bayesian transmission inference analysis.
# The clusters were constructed using R software to group cases
# with a maximum pairwise SNP distance threshold of 50 SNPs,
# i.e. a case was included in a cluster if there were ≤50 SNP
# differences from one or more cases within this cluster. The
# threshold to determine potential transmission links between
# patients has previously typically been fixed between 5 and 12
# SNP differences [6, 15, 26], although recent analysis of within
# patient variation has suggested this estimate may be too low
# [27]. We chose to define broad initial clusters for further
# analysis to allow for missing cases and potential transmission
# between more genetically distant cases."

# Arguments to script:
# Study accession, threshold for SNP dist, distance matrix & ID (Plink output) files location,
# output directory for clusters file, output directory for plot file

# Input:
# Square distance matrix generated by Plink, two-col ID list also generated by Plink

# Main steps:
# Clean data - add IDs as row and col names to distance matrix > convert upper triangle to 0 >
# do clustering with hclust(x, method = 'single') > output plot of clustering (tree) >
# pull clusters of samples below/equal to threshold > write groups

# Output:
# png of tree with clusters below hline, list of samples within a cluster

# RUN:
# Assume run from ~/transmission:
# Rscript r_scripts/transmission_clusters.R <study_accession> <threshold> <dm_file>                         <id_file>                                 <output dir clusters> <output dir plot>
# Rscript r_scripts/transmission_clusters.R PRJEB7669         50          dist_and_pca/PRJEB7669.dist.dist  dist_and_pca/PRJEB7669.dist.dist.id       metadata/             plots/

# ---

# Read in distance matrix and IDs
# dm <- read.delim(dm_file, header = F)
# ids <- read.delim(id_file, header = F)

# Only need first col
# ids <- ids[, 1]

# Add row and col names to dm
# row.names(dm) <- ids
# colnames(dm) <- ids

# Convert to 'distance matrix' (take out upper triangle)
# dm[upper.tri(dm, diag = T)] <- 0

# Need to divide matrix by 2 because Plink assumes diploid
# dm <- dm/2

# Check
# dm[1:10, 1:10]


# Cluster using 'single' method.
# Clusters all samples within distance of threshold (i.e. transitively)
# clust <- hclust(as.dist(dm), method = "single")

# Plot and save
# wh <- 2500
# png(output_plot_file, res = 300, width = wh, height = wh)
# plot(clust, cex = 0.3)
# abline(a = threshold, b = 0)
# dev.off()

# Cut tree at threshold (12 in example below)
# This numbers each group/cluster (members of the same group/cluster have the same number
# so any sample with its own number is not in a group)
# e.g.
# A C B F D E G H I J
# 1 1 2 2 3 4 4 5 6 7
# Groups/clusters are: A-C, B-F, E-G
# clusters <- sort(cutree(clust, h = threshold))

# Filter for clusters (take out the samples with their own number, i.e., those not in a group)
# clusters <- clusters[clusters %in% names(which(table(clusters) > 2))]

# dm_no_single <- dm[names(clusters), names(clusters)]
# 
# clust_no_single <- hclust(as.dist(dm_no_single), method = "single")
# 
# clusters_no_single <- sort(cutree(clust_no_single, h = threshold))
# increase <- 1
# while(any(table(clusters_no_single) == 2)){
#   increase <- increase + 1
#   clusters_no_single <- sort(cutree(clust_no_single, h = threshold+increase))
# }
# threshold+increase
# table(clusters_no_single)

# Find two-sample clusters (sample names)
# two_samp_clusts <- names(clusters[clusters %in% names(which(table(clusters) == 2))])

# Convert to dataframe
# cluster_table <- data.frame(id = names(clusters), cluster = clusters)

# Save
# write.table(cluster_table, file = output_clusters_file, quote = F, row.names = F, sep = "\t")

# ---------------------
# APPENDIX - examples
# ---------------------

# threshold <- 12
#
# set.seed(4)
#
# # Make matrix of random numbers
# x <- sample(1:100, 100, replace = T)
# m <- matrix(x, nrow = 10, ncol = 10)
#
# # Add dummy sample names
# colnames(m) <- LETTERS[1:10]
# row.names(m) <- LETTERS[1:10]
#
# # Convert to 'distance matrix' (take out upper triangle)
# m[upper.tri(m, diag = T)] <- 0

# # Check
# m
#
# # At this point we can look at the matrices for 5 seeds and manually predict the clusters (check against results):
#
# # set.seed(1)
# # A-C, B-F, E-G
#
# # set.seed(2)
# # A-D-F-H-B-E-J
#
# # set.seed(3)
# # A-C-G-B-I-J-D-F-E-H
#
# # set.seed(4)
# # A-C-B-I, E-F-G
#
# # set.seed(5)
# # A-J-D-H-G-I
#
# # Cluster using 'single' method.
# # Clusters all samples within distance of threshold (i.e. transitively)
# clust <- hclust(as.dist(m), method = "single")
#
# # Visualise
# plot(clust)
# abline(a = threshold, b = 0)
#
# # Cut tree at threshold (12 in example below)
# # This numbers each group/cluster (members of the same group/cluster have the same number
# # so any sample with its own number is not in a group)
# # e.g.
# # A C B F D E G H I J
# # 1 1 2 2 3 4 4 5 6 7
# # Groups/clusters are: A-C, B-F, E-G
# clusters <- sort(cutree(clust, h = threshold))
#
# # Filter for clusters (take out the samples with their own number, i.e., those not in a group)
# clusters <- clusters[clusters %in% which(table(clusters) > 1)]
#
# # Convert to dataframe
# cluster_table <- data.frame(id = names(clusters), cluster = clusters)
