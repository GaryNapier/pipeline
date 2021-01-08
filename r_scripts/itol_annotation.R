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
# Rscript r_scripts/itol_annotation.R <study_accession> <metadata_file> <clusters_data_file_dir> <itol_location>
# Rscript r_scripts/itol_annotation.R PRJEB7669 ~/metadata/tb_data_collated_28_10_2020_clean.csv metadata/ itol_annotations

# Setup ----

library(rvest)

heaD <- function(x,...){
  head(x, ...)
}

len_str <- function(string){
  length(unlist(strsplit(string, split = "")))
}

# Read in args
args <- commandArgs(trailingOnly=TRUE)

# Variables
study_accession <- args[1]
nl <- cat("\n")

# Directories
clusters_data_file_dir <- args[3]
itol_location <- args[5]
tmp_dir <- "tmp/"

# Files and suffixes/prefixes
metadata_file <- args[2]

clusters_data_file <- paste0(clusters_data_file_dir, study_accession, ".clusters")

itol_dr_in_file <- paste0(itol_location, study_accession, ".dr.txt")
itol_clusters_in_file <- paste0(itol_location, "itol.clusters.txt")
itol_lineage_in_file <- paste0(itol_location, study_accession, "itol.lineages.txt")

itol_dr_out_file <- paste0(itol_location, study_accession, "itol.dr.txt")
itol_clusters_out_file <- paste0(itol_location, study_accession, ".clusters.txt")
itol_lineage_out_file <- paste0(itol_location, study_accession, "lineages.txt")

nl
print("ARGUMENTS:")
nl
print(c("Study accession:", study_accession))
print(c("clusters_data_file_dir:", clusters_data_file_dir))
print(c("output_location:", output_location))
print(c("tmp_dir:", tmp_dir))
print(c("metadata_file:", metadata_file))
print(c("lineage_file:", lineage_file))
print(c("clusters_data_file:", clusters_data_file))
print(c("itol_dr_file:", itol_dr_file))
print(c("itol_clusters_file:", itol_clusters_file))
print(c("itol_lineage_file:", itol_lineage_file))



# Read in sample data/metadata and clean ----

metadata <- read.csv(metadata_file, stringsAsFactors = F)
clusters_data <- read.table(clusters_data_file, header = T, stringsAsFactors = F)

# Get list of samples from clusters data
samples <- clusters_data$id

# Drug resistance ----

# drugs <- c("rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "streptomycin", "ofloxacin", 
#            "moxifloxacin", "levofloxacin", "amikacin", "kanamycin", "capreomycin", "ciprofloxacin",
#            "prothionamide", "ethionamide", "clarithromycin", "clofazimine", "bedaquiline", "cycloserine", 
#            "linezolid", "para.aminosalicylic_acid", "rifabutin", "delamanid") 

cols <- c("run_accession", "study_accession", "dr_status")

# Subset rows by study accession and cols by run_accession (sample ID), study_accession and drugs
dr_data <- subset(metadata, study_accession == study_accession)
dr_data <- drug_resistance_data$dr_status

sus <- ifelse(dr_data == "Susceptible", 1, -1)
DR <- ifelse(dr_data == "DR", 1, -1)
MDR <- ifelse(dr_data == "MDR", 1, -1)
XDR <- ifelse(dr_data == "XDR", 1, -1)

dr_df <- cbind(metadata$run_accession, sus, DR, MDR, XDR)

# Write dataframe to template


write.table(dr_df, file = itol_dr_out_file, 
            append = T, sep="\t",
            row.names=F, col.names=F, quote = F)






















 