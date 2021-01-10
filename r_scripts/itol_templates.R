#!/usr/bin/env Rscript

# itol_templates.R

# Create templates for itol annotation for:
# - Drug resistance 
# - Clusters (SNP distance thresholds)
# - Lineage


# Arguments to script:
# Location to save files to

# Input:
# Webscrape from itol URLs

# Main steps:
# Scrape the annotation template from itol and clean for each category

# Output:
# Three itol files, one for each type of annotation

# RUN
# Rscript r_scripts/itol_templates.R <itol_templates_location>
# Rscript r_scripts/itol_templates.R itol_annotations/

# Setup ----

library(rvest)

heaD <- function(x,...){
  head(x, ...)
}

len_str <- function(string){
  length(unlist(strsplit(string, split = "")))
}

# Read in args

# Variables

# Directories
itol_templates_location <- "itol_annotations/"

# Files and suffixes/prefixes
itol_dr_file <- paste0(itol_templates_location, "itol.dr.txt")
itol_clusters_file <- paste0(itol_templates_location, "itol.clusters.txt")
itol_lineage_file <- paste0(itol_templates_location, "itol.lineages.txt")

print("ARGUMENTS:")
print("")
print(c("itol_templates_location:", itol_templates_location))
print(c("itol_dr_file:", itol_dr_file))
print(c("itol_clusters_file:", itol_clusters_file))
print(c("itol_lineage_file:", itol_lineage_file))


# Scrape binary template from itol and clean ----

# Scrape raw
binary_template_url <- "https://itol.embl.de/help/dataset_binary_template.txt"
binary_template_text <- html_text(read_html(binary_template_url))

# Uncomment SEPARATOR TAB line and comment other two lines
binary_template_text <- gsub("#SEPARATOR TAB", "SEPARATOR TAB", binary_template_text)
binary_template_text <- gsub("SEPARATOR SPACE", "#SEPARATOR SPACE", binary_template_text)
binary_template_text <- gsub("SEPARATOR COMMA", "#SEPARATOR COMMA", binary_template_text)
# Replace comma with tab
binary_template_text <- gsub("COLOR,#ff0000", "COLOR\t#ff0000", binary_template_text)


# Scrape colourstrip template from itol

strip_template_url <- "https://itol.embl.de/help/dataset_color_strip_template.txt"
strip_template_text <- html_text(read_html(strip_template_url))

strip_template_text <- gsub("#SEPARATOR TAB", "SEPARATOR TAB", strip_template_text)
strip_template_text <- gsub("SEPARATOR SPACE", "#SEPARATOR SPACE", strip_template_text)


strip_template_text <- gsub("COLOR #ff0000", "COLOR\t#ff0000", strip_template_text)
strip_template_text <- gsub("COLOR_BRANCHES 0", "COLOR_BRANCHES\t0", strip_template_text)


# Scrape range template from itol

range_template_url <- "https://itol.embl.de/help/colors_styles_template.txt"
range_template_text <- html_text(read_html(range_template_url))

range_template_text <- gsub("#SEPARATOR TAB", "SEPARATOR TAB", range_template_text)
range_template_text <- gsub("SEPARATOR SPACE", "#SEPARATOR SPACE", range_template_text)



# Drug resistance ----

legend_text <- "LEGEND_TITLE\tDrug susceptibility\nLEGEND_SHAPES\t1\t2\t3\t4\nLEGEND_COLORS\t#e8837d\t#ace87d\t#7de2e8\t#b97de8\nLEGEND_LABELS\tSusceptible\tDR\tMDR\tXDR"

# Replace commas with tabs
dr_template_text <- gsub("DATASET_LABEL,label1", "DATASET_LABEL\tDrug_resistance", binary_template_text)
dr_template_text <- gsub("FIELD_SHAPES,1", "FIELD_SHAPES\t1\t2\t3\t4", dr_template_text)
dr_template_text <- gsub("FIELD_LABELS,f1", "FIELD_LABELS\tSusceptible\tDR\tMDR\tXDR", dr_template_text)
dr_template_text <- gsub("#FIELD_COLORS,#ff0000,#00ff00,#ffff00,#0000ff", "", dr_template_text)
dr_template_text <- gsub("#FIELD_COLORS,#ff0000", "FIELD_COLORS\t#e8837d\t#ace87d\t#7de2e8\t#b97de8", dr_template_text)
dr_template_text <- gsub("#LEGEND_SHAPE_SCALES,1,1,0.5", legend_text, dr_template_text)

write.table(dr_template_text, file = itol_dr_file,
            row.names=F,col.names=F, quote = F)


# Clusters ----

clusters_template_text <- gsub("DATASET_LABEL label1", "DATASET_LABEL\tClusters", strip_template_text)

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

clusters_template_text <- gsub("#LEGEND_TITLE Dataset_legend", "LEGEND_TITLE\tClusters", clusters_template_text)

write.table(clusters_template_text, file = itol_clusters_file,
            row.names=F,col.names=F, quote = F)


# Lineage ----

write.table(range_template_text, file = itol_lineage_file,
            row.names=F,col.names=F, quote = F)













