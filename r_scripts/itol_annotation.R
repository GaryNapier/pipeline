#!/usr/bin/env Rscript

# itol_annotation.R

# Annotate trees for each study accession (project) with:
# - Drug resistance (from metadata)
# - Clusters based on min SNP distance - output of r_scripts/transmission_clusters.R -> metadata/<study_accession>.clusters)
# - Lineage - output of tbprofiler

# Arguments to script:
# Study accession, DR metadata file name (full path), clusters file location (dir only), lineage file name (full path), output_location

# Input:
# Scraped annotation template from itol, DR metadata file, clusters file, lineage file

# Main steps:
# Scrape the annotation template from itol and clean
# Filter DR data by study accession, clean and parse into dataframe ready for itol template
# Get cluster membership of samples and onvert to colour scheme - parse into dataframe ready for itol template
# Repeat for lineage

# Output:
# Three itol files, one for each type of annotation

# RUN:
# Rscript r_scripts/itol_annotation.R <study_accession> <drug_resistance_data_file> <clusters_data_file_dir> <lineage_file> <output_location>
# Rscript r_scripts/itol_annotation.R PRJEB7669 ~/metadata/tb_data_collated_28_10_2020_clean.csv metadata/ ~/metadata/lineages_napier_barcode.txt itol_annotations

# Setup ----

# Read in args
args <- commandArgs(trailingOnly=TRUE)

# Variables
study_accession <- args[1]
nl <- cat("\n")

# Directories
clusters_data_file_dir <- args[3]
output_location <- args[5]
tmp_dir <- "tmp/"

# Files and suffixes/prefixes
drug_resistance_data_file <- args[2]
lineage_file <- args[4]
itol_template_file <- paste0(tmp_dir, "itol_template.txt")
clusters_data_file <- paste0(clusters_data_file_dir, study_accession, ".clusters")
itol_dr_file <- paste0(output_location, study_accession, ".dr.txt")
itol_clusters_file <- paste0(output_location, study_accession, ".clusters.txt")
itol_lineage_file <- paste0(output_location, study_accession, "lineages.txt")

nl
print("ARGUMENTS:")
nl
print(c("Study accession:", study_accession))
print(c("clusters_data_file_dir:", clusters_data_file_dir))
print(c("output_location:", output_location))
print(c("tmp_dir:", tmp_dir))
print(c("drug_resistance_data_file:", drug_resistance_data_file))
print(c("lineage_file:", lineage_file))
print(c("itol_template_file:", itol_template_file))
print(c("clusters_data_file:", clusters_data_file))
print(c("itol_dr_file:", itol_dr_file))
print(c("itol_clusters_file:", itol_clusters_file))
print(c("itol_lineage_file:", itol_lineage_file))


# Scrape binary template from itol and clean

# Scrape raw
binary_template_url <- "https://itol.embl.de/help/dataset_binary_template.txt"
binary_template_text <- html_text(read_html(binary_template_url))

# Uncomment SEPARATOR TAB line and comment other two lines
binary_template_text <- gsub("#SEPARATOR TAB", "SEPARATOR TAB", binary_template_text)
binary_template_text <- gsub("SEPARATOR SPACE", "#SEPARATOR SPACE", binary_template_text)
binary_template_text <- gsub("SEPARATOR COMMA", "#SEPARATOR COMMA", binary_template_text)

# Replace commas with tabs
binary_template_text <- gsub("DATASET_LABEL,label1", "DATASET_LABEL\told/new", binary_template_text)
binary_template_text <- gsub("COLOR,#ff0000", "COLOR\t#ff0000", binary_template_text)
binary_template_text <- gsub("FIELD_SHAPES,1", "FIELD_SHAPES\t5\t2", binary_template_text)
binary_template_text <- gsub("FIELD_LABELS,f1", "FIELD_LABELS\told\tnew", binary_template_text)
binary_template_text <- gsub("#FIELD_COLORS,#ff0000", "FIELD_COLORS\t#ff0000\t#000000", binary_template_text)
# binary_template_text <- gsub("", "", binary_template_text)
# binary_template_text <- gsub("", "", binary_template_text)


write.table(binary_template_text, file = itol_template_file, 
            row.names=F,col.names=F, quote = F)

# Read in sample data/metadata and clean ----

drug_resistance_data <- read.csv(drug_resistance_data_file, stringsAsFactors = F)
clusters_data <- read.table(clusters_data_file, header = T, stringsAsFactors = F)
lineage_data <- read.delim(lineage_file, stringsAsFactors = F)

# Get list of samples from clusters data
samples <- clusters_data$id

# Drug resistance clean/subset:

drugs <- c("rifampicin", "isoniazid", "ethambutol", "pyrazinamide", "streptomycin", "ofloxacin", 
           "moxifloxacin", "levofloxacin", "amikacin", "kanamycin", "capreomycin", "ciprofloxacin",
           "prothionamide", "ethionamide", "clarithromycin", "clofazimine", "bedaquiline", "cycloserine", 
           "linezolid", "para.aminosalicylic_acid", "rifabutin", "delamanid") 

cols <- c("run_accession", "study_accession", drugs)

# Subset rows by study accession and cols by run_accession (sample ID), study_accession and drugs
drug_resistance_data <- subset(drug_resistance_data, study_accession == study_accession)
drug_resistance_data <- drug_resistance_data[, cols]





















 